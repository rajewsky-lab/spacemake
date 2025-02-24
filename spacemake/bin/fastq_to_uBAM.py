#!/usr/bin/env python3
from spacemake.contrib import __version__, __author__, __license__, __email__
from collections import defaultdict
from time import time
import numpy as np
import spacemake.util as util
import mrfifo as mf


class SeqData(object):
    """
    Minimal Sequencing Read Data container. Holds the following attributes
        .qname
        .r1
        .q1
        .r2
        .q2
        .tags: `defaultdict(list)`

    qname, seq and qual (1 and 2) map directly onto the paired-end FASTQ input. `tags` is a
    dictionary with BAM tag names as keys and list as values to which every processing
    function may `append` strings. The value of the tag in the end is a comma-separated
    string version of the list elements.
    """

    __slots__ = ("qname", "r1", "q1", "r2", "q2", "tags")

    def __init__(self, qname, r1, q1, r2, q2):
        self.qname = qname
        self.r1 = r1
        self.q1 = q1
        self.r2 = r2
        self.q2 = q2
        self.tags = defaultdict(list)

    def render_tags(self, extra_tags=[]):
        # print(self.tags)
        bam_tags = [(n, f"{','.join(records)}") for n, records in self.tags.items()]
        tag_str = "\t".join([f"{tag}:Z:{val}" for tag, val in bam_tags] + extra_tags)
        return tag_str

    def render_SAM(self, flag=4, extra_tags=["RG:Z:A"]):
        tag_str = self.render_tags(extra_tags=extra_tags)
        # if len(self.r2) != len(self.q2):
        #     print(self.r2, self.q2, tag_str)
        return f"{self.qname}\t{flag}\t*\t0\t0\t*\t*\t0\t0\t{self.r2}\t{self.q2}\t{tag_str}\n"

    def __str__(self):
        return self.render_SAM()

    @classmethod
    def from_single_end(cls, fq2):
        fq_src2 = util.FASTQ_src(fq2)
        for fqid, seq, qual in fq_src2:
            yield cls(fqid, "NA", "##", seq, qual)

    @classmethod
    def from_paired_end(cls, fq1, fq2):
        fq_src1 = util.FASTQ_src(fq1)
        fq_src2 = util.FASTQ_src(fq2)
        for (fqid, r1, q1), (_, r2, q2) in zip(fq_src1, fq_src2):
            yield cls(fqid, r1, q1, r2, q2)

    @classmethod
    def from_BAM(cls, bam_src):
        "TODO: allow pre-processing of raw data that is already in BAM format"
        pass


## Preprocessing is done by sequentially acting on a SeqData object to modify
# it in-place. We use factory functions that return the actual function which
# does that, because it allows to have out-of-loop setup code.


def quantize_quality(bin_width=10):
    QMAP = {}
    for Q in range(41):
        C = chr(Q + 33)
        q = int(np.round(Q / float(bin_width), 0)) * bin_width
        c = chr(q + 33)
        QMAP[C] = c

    def qquant(sdata):
        sdata.q2 = "".join([QMAP[q] for q in sdata.q2])

    return qquant


def name_simplifier(n_workers=8, w=0, chunk_size=10):
    i = 0

    def simplify(sdata):
        nonlocal i
        batch = i // chunk_size
        within_chunk = i % chunk_size
        round = batch // n_workers
        done = round * n_workers * chunk_size
        n = done + w * chunk_size + within_chunk
        sdata.qname = f"{int(n):d}"
        i += 1

    return simplify


def clip(left=0, right=0):
    def c(sdata):
        if left:
            sdata.r2 = sdata.r2[left:]
            sdata.q2 = sdata.q2[left:]
        if right:
            sdata.r2 = sdata.r2[:-right]
            sdata.q2 = sdata.q2[:-right]

    return c


def nextseq_trim(cutoff=30):
    from cutadapt.qualtrim import nextseq_trim_index
    from dnaio import SequenceRecord

    def Q(sdata):
        end = len(sdata.r2)
        sr = SequenceRecord(name=sdata.qname, sequence=sdata.r2, qualities=sdata.q2)
        # we use the cutadapt function here (which implements BWA's logic).
        new_end = nextseq_trim_index(sr, cutoff=cutoff)
        n_trimmed = len(sdata.q2) - new_end
        if n_trimmed:
            sdata.q2 = sdata.q2[:new_end]
            sdata.r2 = sdata.r2[:new_end]

            sdata.tags["A3"].append("NextQ")
            sdata.tags["T3"].append(str(end - new_end))

    return Q


def qual_trim(left=0, right=25):
    from cutadapt.qualtrim import quality_trim_index

    def Q(sdata):
        end = len(sdata.r2)
        # we use the cutadapt function here (which implements BWA's logic).
        new_start, new_end = quality_trim_index(sdata.q2, left, right)
        n_trimmed = len(sdata.q2) - (new_end - new_start)
        if n_trimmed:
            sdata.q2 = sdata.q2[new_start:new_end]
            sdata.r2 = sdata.r2[new_start:new_end]
            if new_start:
                sdata.tags["A5"].append("Q")
                sdata.tags["T5"].append(str(new_start))

            if new_end != end:
                sdata.tags["A3"].append("Q")
                sdata.tags["T3"].append(str(end - new_end))

    return Q


def polyA_trim(rev_comp=False):
    from cutadapt.qualtrim import poly_a_trim_index

    def polyA(sdata):
        end = len(sdata.r2)
        index = poly_a_trim_index(sdata.r2, revcomp=rev_comp)
        # print(
        #     f"PolyA: seq='{seq}' index={index} end={end} seq[:index]={seq[:index]} seq[index:]={seq[index:]}"
        # )

        n_trimmed = end - index
        if n_trimmed:
            sdata.r2 = sdata.r2[:index]
            sdata.q2 = sdata.q2[:index]
            sdata.tags["A3"].append("polyA")
            sdata.tags["T3"].append(str(n_trimmed))

    return polyA


def adapter_trim(
    name="SMART_TSO", seq="AAGCAGTGGTATCAACGCAGAGTGAATGGG", where="right", **kw
):
    import cutadapt.adapters

    if where == "left":
        adapter = cutadapt.adapters.NonInternalFrontAdapter(seq, name=name, **kw)
    else:
        adapter = cutadapt.adapters.BackAdapter(seq, name=name, **kw)

    def adap(sdata):
        match = adapter.match_to(sdata.r2)
        if match:
            if where == "left":
                n_trimmed = match.rstop
                sdata.tags["A5"].append(name)
                sdata.tags["T5"].append(str(n_trimmed))
                sdata.r2 = sdata.r2[match.rstop :]
                sdata.q2 = sdata.q2[match.rstop :]
            else:
                n_trimmed = len(sdata.r2) - match.rstart
                sdata.tags["A3"].append(name)
                sdata.tags["T3"].append(str(n_trimmed))
                sdata.r2 = sdata.r2[: match.rstart]
                sdata.q2 = sdata.q2[: match.rstart]

    return adap


format_func_template = """
def format_func(sdata):
    {i5i7}
    r1 = sdata.r1
    q1 = sdata.q1
    r2 = sdata.r2
    r2_qual = sdata.q2

    sdata.tags['CB'] = [{cell}, ]
    sdata.tags['MI'] = [{UMI}, ]

    sdata.r2 = {seq}
    sdata.q2 = {qual}
"""
i5i7_extract_template = "i5i7 = sdata.qname.split(':N:0:')[1].replace('+', '')"


def barcode(cell="r1[:12]", UMI="r1[12:20]", seq="r2", qual="q2", disable_safety=False):
    def safety_check_eval(s, danger="();."):
        chars = set(list(s))
        if chars & set(list(danger)):
            return False
        else:
            return True

    if not disable_safety:
        assert safety_check_eval(cell)
        assert safety_check_eval(UMI)
        assert safety_check_eval(seq)
        assert safety_check_eval(qual)

    # This trickery here compiles a python function that
    # embeds the user-specified code for cell, UMI
    # we replace globals() with an empty dictionary d,
    i5i7 = i5i7_extract_template if "i5i7" in cell else ""

    d = {}
    exec(
        format_func_template.format(cell=cell, UMI=UMI, seq=seq, qual=qual, i5i7=i5i7),
        d,
    )
    # from said dict we can now retrieve the function object
    # and call it many times over
    func = d["format_func"]
    return func


class PreProcessor(object):
    type_dict = {
        "right": int,
        "left": int,
        "min_overlap": int,
        "bin_width": int,
        "cutoff": int,
        "max_errors": float,
        "rev_comp": bool,
        "disable_safety": bool,
    }
    fn_dict = {
        "quality": qual_trim,
        "nextseq_quality": nextseq_trim,
        "polyA": polyA_trim,
        "adapter": adapter_trim,
        "qquant": quantize_quality,
        "simple_name": name_simplifier,
        "barcode": barcode,
        "clip": clip,
    }

    def __init__(self, processing_str="", flavor_dict={}, **kw):  # , min_len=18):
        self.kw = kw
        # collect statistics here
        self.stats = defaultdict(int)

        # parse processing_str and assemble a pre-processing pipeline
        if flavor_dict:
            self.pipeline = self.pipeline_from_flavor(flavor_dict)
        else:
            self.pipeline = self.pipeline_from_str(processing_str)

    def pipeline_from_flavor(self, flavor_dict):
        pipeline = []
        from pprint import pprint

        for step_d in flavor_dict:
            # this should always only have one item
            for fn_name, kw_d in step_d.items():
                kw = self.kw.get(fn_name, {}).copy()
                if kw_d is not None:
                    kw.update(kw_d)

                fn = self.fn_dict[fn_name]
                pipeline.append(fn(**kw))

        return pipeline

    def pipeline_from_str(self, processing_str):
        pipeline = []

        for token in processing_str.split(";"):
            if ":" in token:
                fn_name, kw_str = token.split(":", 1)
            else:
                fn_name = token
                kw_str = ""

            kw = self.kw.get(fn_name, {}).copy()
            for param in kw_str.split(","):
                if not param:
                    continue

                k, v_str = param.split("=")
                kw[k] = self.type_dict.get(k, str)(v_str)

            fn = self.fn_dict[fn_name]
            pipeline.append(fn(**kw))

        return pipeline

    def process(self, sdata):
        # columns of the stats:
        # (sample) measure obsname value count
        self.stats[("reads", "N", "input")] += 1
        self.stats[("bases", "L_in", len(sdata.r2))] += 1
        for func in self.pipeline:
            func(sdata)  # modifies sdata in-place
            # print(f"after func {func}: {sdata}")

        self.record_tag_stats(sdata.tags)
        self.stats[("bases", "L_out", len(sdata.r2))] += 1

        return sdata

    def record_tag_stats(self, tags):
        if "A5" in tags:
            for a, t in zip(tags["A5"], tags["T5"]):
                self.stats[("reads", "A5", a)] += 1
                self.stats[(f"bases", a, t)] += 1
                self.stats[(f"bases", a, "A5_total")] += int(t)

        if "A3" in tags:
            for a, t in zip(tags["A3"], tags["T3"]):
                self.stats[("reads", "A3", a)] += 1
                self.stats[(f"bases", a, t)] += 1
                self.stats[(f"bases", a, "A3_total")] += int(t)


def process_fastq(fq1, fq2, sam_out, args, _extra_args={}, **kwargs):

    logger = util.setup_logging(args, "fastq_to_uBAM.worker", rename_process=False)
    logger.debug(
        f"starting up with fq1={fq1}, fq2={fq2} sam_out={sam_out} and args={args}"
    )

    ingress = SeqData.from_paired_end(fq1, fq2) if fq1 else SeqData.from_single_end(fq2)

    # TODO:
    # decide how we want to handle multiple input files
    # somehow I prefer sequential processing over complicated state-change
    # question is, how do we keep the output sam from being closed?
    # ideally, insert a serializer there and provide a way to join() some processes before starting new
    # ones... (needs thinking)
    # kw = vars(args)
    # # kw.update(par)
    # import argparse
    # args = argparse.Namespace(**kw)
    if args.flavor:
        args = util.load_config_with_fallbacks(args)
        flavor_dict = args.config["preprocessing_flavors"][args.flavor]
    else:
        flavor_dict = {}

    pre = PreProcessor(
        processing_str=args.processing,
        flavor_dict=flavor_dict,
        simple_name=dict(
            # w=_extra_args["n"], # needed for read simplification but not yet in mrfifo upstream. Disable for now
            n_workers=args.threads_work,
            chunk_size=args.chunk_size,
        ),
        # LEGACY CONFIGURATION
        barcode=dict(
            cell=args.cell,
            UMI=args.UMI,
            seq=args.seq,
            qual=args.qual,
        ),
    )

    # N = mf.util.CountDict()

    for sdata in ingress:
        # N.count("total")
        sdata = pre.process(sdata)
        sam_out.write(sdata.render_SAM(flag=4))

    return pre.stats


def min_length_filter(input, output, min_len=18):
    N = 0
    for line in input:
        if not line.startswith("@"):
            seq = line.split("\t")[9]
            if len(seq) < min_len:
                continue

        output.write(line)
        N += 1

    return N
    # if args.paired_end:
    #             # if args.paired_end[0] == 'r':
    #             #     # rev_comp R1 and reverse qual1
    #             #     q1 = q1[::-1]
    #             #     r1 = util.rev_comp(r1)

    #             # if args.paired_end[1] == 'r':
    #             #     # rev_comp R2 and reverse qual2
    #             #     q2 = q2[::-1]
    #             #     r2 = util.rev_comp(r2)

    #             rec1 = fmt.make_bam_record(
    #                 qname=fqid,
    #                 r1="THISSHOULDNEVERBEREFERENCED",
    #                 r2=r1,
    #                 R1=r1,
    #                 R2=r2,
    #                 r2_qual=q1,
    #                 r2_qname=fqid,
    #                 flag=69,  # unmapped, paired, first in pair
    #             )
    #             rec2 = fmt.make_bam_record(
    #                 qname=fqid,
    #                 r1="THISSHOULDNEVERBEREFERENCED",
    #                 r2=r2,
    #                 R1=r1,
    #                 R2=r2,
    #                 r2_qual=q2,
    #                 r2_qname=fqid,
    #                 flag=133,  # unmapped, paired, second in pair
    #             )
    #             results.append(rec1)
    #             results.append(rec2)

    #         else:


def main(args):
    input_reads1, input_reads2, input_params = get_input_params(args)
    have_read1 = set([str(r1) != "None" for r1 in input_reads1]) == set([True])

    # queues for communication between processes
    w = (
        mf.Workflow("fastq_to_uBAM", total_pipe_buffer_MB=args.pipe_buffer)
        # open reads2.fastq.gz
        .gz_reader(inputs=input_reads2, output=mf.FIFO("read2", "wb")).distribute(
            input=mf.FIFO("read2", "rt"),
            outputs=mf.FIFO("r2_{n}", "wt", n=args.threads_work),
            chunk_size=args.chunk_size * 4,
        )
    )

    if have_read1:
        # open reads1.fastq.gz
        w.gz_reader(inputs=input_reads1, output=mf.FIFO("read1", "wb"))
        w.distribute(
            input=mf.FIFO("read1", "rt"),
            outputs=mf.FIFO("r1_{n}", "wt", n=args.threads_work),
            chunk_size=args.chunk_size * 4,  # 4 FASTQ lines per record!
        )
        # process in parallel workers
        w.workers(
            func=process_fastq,
            fq1=mf.FIFO("r1_{n}", "rt"),
            fq2=mf.FIFO("r2_{n}", "rt"),
            sam_out=mf.FIFO("sam_{n}", "wt"),
            args=args,
            n=args.threads_work,
        )
    else:
        # process in parallel workers
        w.workers(
            func=process_fastq,
            fq1=None,
            fq2=mf.FIFO("r2_{n}", "rt"),
            sam_out=mf.FIFO("sam_{n}", "wt"),
            args=args,
            n=args.threads_work,
        )

    # combine output streams
    w.collect(
        inputs=mf.FIFO("sam_{n}", "rt", n=args.threads_work),
        chunk_size=args.chunk_size,
        custom_header=mf.util.make_SAM_header(
            prog_id="fastq_to_uBAM",
            prog_name="fastq_to_uBAM.py",
            prog_version=__version__,
            rg_name=args.sample,
        ),
        # output="/dev/stdout"
        # )
        output=mf.FIFO("sam_combined", "wt"),
        log_rate_every_n=1000000,
        log_rate_template="written {M_out:.1f} M BAM records ({mps:.3f} M/s, overall {MPS:.3f} M/s)",
        log_name="fastq_to_uBAM.collect",
    )
    w.funnel(
        func=min_length_filter,
        input=mf.FIFO("sam_combined", "rt"),
        output=mf.FIFO("sam_filtered", "wt"),
        min_len=args.min_len,
    )
    # compress to BAM
    fmt_opt = " ".join([f"--output-fmt-option {o}" for o in args.out_fmt_option])
    fmt = f"Sh -O {args.out_fmt} {fmt_opt}"

    w.funnel(
        func=mf.parts.bam_writer,  # mf.parts.null_writer, #
        input=mf.FIFO("sam_filtered", "rt"),
        output=args.out_file,
        _manage_fifos=False,
        fmt=fmt,
        threads=args.threads_write,
    )
    res = w.run()
    stats = defaultdict(int)
    for w, d in res.result_dict.items():
        if "worker" in w:
            for k, v in d.items():
                stats[k] += v
        elif "funnel0" in w:
            stats[("reads", "N", "output")] = d

    import pandas as pd

    data = [[args.sample] + list(key) + [value] for key, value in stats.items()]
    # data.append()
    df = pd.DataFrame(
        data, columns=["sample", "measure", "name", "value", "count"]
    ).sort_values(["name", "value"])

    if args.out_stats:
        df.to_csv(args.out_stats, sep="\t", index=None)

    return df


def get_input_params(args):
    import pandas as pd
    import csv

    if str(args.matrix) != "None":
        df = pd.read_csv(args.matrix, sep=",", index_col=None, quoting=csv.QUOTE_NONE)
        if not "R1" in df.columns:
            df["R1"] = "None"

        R1 = df["R1"].to_list()
        R2 = df["R2"].to_list()
        # TODO extract other column values into a list of dictionaries (most importantly cell=...)
        # which can override how the formatter in the workers processes the raw reads
        if "cell" in df.columns:
            # params = [{'cell': f'"{c}"'} for c in df['cell']]
            params = [{"cell": c} for c in df["cell"]]
        else:
            params = [
                {},
            ] * len(R1)

    else:
        R1 = [
            args.read1,
        ]
        R2 = [
            args.read2,
        ]
        params = [
            {},
        ]

    return R1, R2, params


def parse_args():
    import spacemake.util as util

    parser = util.make_minimal_parser(
        "fastq_to_uBAM",
        description="Convert raw reads1 and reads2 FASTQ into a single BAM file with cell barcode and UMI as BAM-tags",
    )
    ## input options
    parser.add_argument(
        "--matrix",
        default=None,
        help="sample_matrix.csv file desribing from where to get read1 and read2 (FASTQ format)",
    )
    parser.add_argument(
        "--read1",
        default=None,
        help="source from where to get read1 (FASTQ format)",
    )
    parser.add_argument(
        "--read2",
        default="/dev/stdin",
        help="source from where to get read2 (FASTQ format)",
        # required=True,
    )
    parser.add_argument(
        "--paired-end",
        default=None,
        choices=["fr", "ff", "rf", "rr"],
        help="read1 and read2 are paired end mates and store both in the BAM",
    )
    parser.add_argument(
        "--phred-base",
        default=33,
        type=int,
        help="phred quality base in the input (default=33)",
    )

    ## pre-processing options
    parser.add_argument(
        "--processing",
        help="string encoding the processing of the cDNA",
        default="quality:right=25;polyA;adapter:name=SMART,seq=AAGCAGTGGTATCAACGCAGAGTGAATGGG,max_errors=0.1,min_overlap=10;barcode",  #:cell=r1[8:20][::-1],UMI=r1[0:8]",
    )
    parser.add_argument(
        "--config",
        help="path to config.yaml",
        default="config.yaml",
    )
    parser.add_argument(
        "--flavor",
        help="read a preprocessing_flavor from config.yaml to configure the pre-processing pipeline, instead of parsing --processing",
        default="",
    )
    parser.add_argument(
        "--min-len",
        default=18,
        type=int,
        help="minimum read2 length to keep after preprocessing, which involves trimming (default=18)",
    )
    # Barcode/UMI extraction LEGACY options
    parser.add_argument("--cell", default="r1[8:20][::-1]", help="DEPRECATED")
    parser.add_argument("--UMI", default="r1[0:8]", help="DEPRECATED")
    parser.add_argument("--seq", default="r2", help="DEPRECATED")
    parser.add_argument("--qual", default="r2_qual", help="DEPRECATED")
    parser.add_argument("--disable-safety", default=False, type=bool, help="DEPRECATED")

    ## parallelization
    parser.add_argument(
        "--pipe-buffer",
        default=4,
        type=int,
        help="How many megabytes of pipe-buffer to use. kernel settings is usually 64MB per user (default=4MB)",
    )
    parser.add_argument(
        "--chunk-size",
        default=10,
        type=int,
        help="how many consecutive reads are assigned to the same worker (default=10)",
    )
    parser.add_argument(
        "--threads-work",
        help="number of worker threads for actual trimming (default=8)",
        type=int,
        default=8,
    )
    parser.add_argument(
        "--threads-write",
        help="number of threads for writing bam_out (default=8)",
        type=int,
        default=8,
    )

    ## output options
    parser.add_argument(
        "--out-file",
        default="/dev/stdout",
        help="output for unaligned BAM records (default=/dev/stdout) ",
    )
    # SAM standard now supports CB=corrected cell barcode, CR=original cell barcode, and MI=molecule identifier/UMI
    parser.add_argument(
        "--out-fmt",
        default="BAM",
        choices=[
            "SAM",
            "BAM",
            "CRAM",
        ],
        help="SAM, BAM (default), or CRAM",
    )
    parser.add_argument(
        "--out-fmt-option",
        default="",
        nargs="+",
        help="passed options to `samtools view` as --output-fmt-options",
    )
    parser.add_argument(
        "--out-stats",
        default="preprocessing_stats.txt",
        help="store statistics in this file",
    )
    parser.add_argument(
        "--out-cell-barcodes",
        default="",
        help="store (raw) cell barcode counts in this file. Numbers add up to number of total raw reads. WARNING: may use excessive RAM for Open-ST data!",
    )

    args = parser.parse_args()
    return args


def cmdline():
    from time import time

    args = parse_args()
    logger = util.setup_logging(args, name="spacemake.bin.fastq_to_uBAM")
    for k, v in sorted(vars(args).items()):
        logger.info(f"cmdline arg\t{k}={v}")

    t0 = time()
    df = main(args)
    dt = time() - t0
    N = df.query("name == 'N' and value == 'input'")["count"].iloc[0]
    N_kept = df.query("name == 'N' and value == 'output'")["count"].iloc[0]
    print(N)
    rate = N / dt / 1000
    logger.info(
        f"processed {N/1e6:.3f} M reads in {dt:.1f} seconds ({rate:.1f} k reads/sec)"
    )
    logger.info(f"kept {N_kept/1e6:.3f} M reads in output ({100 * N_kept/N:.2f} %)")
    return df


if __name__ == "__main__":
    res = cmdline()
