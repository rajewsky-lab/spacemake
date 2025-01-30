#!/usr/bin/env python3
from spacemake.contrib import __version__, __author__, __license__, __email__
import logging
import os
import sys
import numpy as np
import multiprocessing as mp
from collections import defaultdict
from time import time
import spacemake.util as util
import mrfifo as mf

# TODO:
# * load params from YAML (code from opseq)
# * save params to YAML in run-folder to document
# * refactor into multiple modules?

format_func_template = """
import re

def format_func(qname=None, r2_qname=None, r2_qual=None, r1=None, r2=None, R1=None, R2=None, **kw):
    #qparts = r2_qname.split(':N:0:')
    #if len(qparts) > 1:
    #    i5i7 = qparts[1].replace('+', '')
    #else:
    #    i5i7 = None

    cell = {cell}

    attrs = dict(
        cell=cell,
        raw=cell,
        UMI={UMI},
        seq={seq},
        qual={qual},
        fqid=r2_qname,
    )

    return attrs
"""


def safety_check_eval(s, danger="();."):
    chars = set(list(s))
    if chars & set(list(danger)):
        return False
    else:
        return True


def make_formatter_from_args(args):
    if not args.disable_safety:
        assert safety_check_eval(args.cell)
        assert safety_check_eval(args.UMI)
        assert safety_check_eval(args.seq)
        assert safety_check_eval(args.qual)

    # This trickery here compiles a python function that
    # embeds the user-specified code for cell, UMI
    # we replace globals() with an empty dictionary d,
    d = {}
    exec(
        format_func_template.format(
            cell=args.cell, UMI=args.UMI, seq=args.seq, qual=args.qual
        ),
        d,
    )
    # from said dict we can now retrieve the function object
    # and call it many times over
    func = d["format_func"]
    return func


def make_sam_record(
    fqid,
    seq,
    qual,
    flag=4,
    tags=[("CB", "{cell}"), ("MI", "{UMI}"), ("CR", "{raw}"), ("RG", "A")],
    **kw,
):
    tag_str = "\t".join([f"{tag}:Z:{tstr.format(**kw)}" for tag, tstr in tags])
    return f"{fqid}\t{flag}\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}\t{tag_str}\n"


def simplify_qname(qname, n):
    return f"{int(n):d}"


QMAP = {}
for Q in range(41):
    C = chr(Q + 33)
    q = int(np.round(Q / 10, 0)) * 10
    c = chr(q + 33)
    QMAP[C] = c


def quantize_quality(qual):
    return "".join([QMAP[q] for q in qual])


# pre-processing function factories
def qual_trim(left=0, right=25):
    from cutadapt.qualtrim import quality_trim_index

    def Q(qname, seq, qual, tags):
        end = len(seq)
        # we use the cutadapt function here (which implements BWA's logic).
        new_start, new_end = quality_trim_index(qual, left, right)
        n_trimmed = len(qual) - (new_end - new_start)
        if n_trimmed:
            qual = qual[new_start:new_end]
            seq = seq[new_start:new_end]
            if new_start:
                tags["A5"].append("Q")
                tags["T5"].append(str(new_start))

            if new_end != end:
                tags["A3"].append("Q")
                tags["T3"].append(str(end - new_end))

        return qname, seq, qual, tags

    return Q


def polyA_trim(rev_comp=True):
    from cutadapt.qualtrim import poly_a_trim_index

    def polyA(qname, seq, qual, tags):
        end = len(seq)
        index = poly_a_trim_index(seq, revcomp=rev_comp)

        n_trimmed = end - index
        if n_trimmed:
            seq = seq[:index]
            qual = qual[:index]
            tags["A3"].append("polyA")
            tags["T3"].append(str(n_trimmed))

        return qname, seq, qual, tags

    return polyA


def adapter_trim(
    name="SMART_TSO", seq="AAGCAGTGGTATCAACGCAGAGTGAATGGG", where="right", **kw
):
    import cutadapt.adapters

    if where == "left":
        adapter = cutadapt.adapters.NonInternalFrontAdapter(seq, name=name, **kw)
    else:
        adapter = cutadapt.adapters.BackAdapter(seq, name=name, **kw)

    def adap(qname, qseq, qual, tags):
        match = adapter.match_to(qseq)
        if match:
            if where == "left":
                n_trimmed = match.rstop
                tags["A5"].append(name)
                tags["T5"].append(n_trimmed)
                qseq = qseq[match.rstop :]
                qual = qual[match.rstop :]
            else:
                n_trimmed = len(seq) - match.rstart
                tags["A3"].append(name)
                tags["T3"].append(n_trimmed)
                qseq = qseq[: match.rstart]
                qual = qual[: match.rstop]

        return qname, qseq, qual, tags

    return adap


class PreProcessor(object):
    def __init__(self, processing_str):  # , min_len=18):
        # parse processing_str and assemble a processor function
        self.pipeline = self.parse(processing_str)
        # collect statistics here
        self.stats = defaultdict(lambda: defaultdict(int))
        # self.min_len = min_len

    def parse(self, processing_str):
        type_dict = {"right": int, "left": int, "min_overlap": int, "max_errors": float}
        fn_dict = {
            "Q": qual_trim,
            "polyA": polyA_trim,
            "adapter": adapter_trim,
            # TODO: add 'clip', 'barcode', 'qquant', 'simple_name'
        }

        pipeline = []

        for token in processing_str.split(";"):
            if ":" in token:
                fn_name, kw_str = token.split(":", 1)
            else:
                fn_name = token
                kw_str = ""

            kw = {}
            for param in kw_str.split(","):
                if not param:
                    continue

                k, v_str = param.split("=")
                kw[k] = type_dict.get(k, str)(v_str)

            fn = fn_dict[fn_name]
            pipeline.append(fn(**kw))

        return pipeline

    def process(self, qname, seq, qual):
        tags = defaultdict(list)
        self.stats["len_in"][len(seq)] += 1
        for func in self.pipeline:
            qname, seq, qual, tags = func(qname, seq, qual, tags)

        self.record_stats(tags)
        self.stats["len_out"][len(seq)] += 1
        bam_tags = [(n, f"{','.join(records)}") for n, records in tags.items()]
        return qname, seq, qual, bam_tags

    def record_stats(self, tags):
        for a, t in zip(tags["A5"], tags["T5"]):
            self.stats[f"freq"][a] += 1
            self.stats[f"bases_{a}"][t] += 1
            self.stats[f"bases"][a] += int(t)

        for a, t in zip(tags["A3"], tags["T3"]):
            self.stats[f"freq"][a] += 1
            self.stats[f"bases_{a}"][t] += 1
            self.stats[f"bases"][a] += int(t)


def render_to_sam(fq1, fq2, sam_out, args, _extra_args={}, **kwargs):

    w = _extra_args["n"]

    logger = util.setup_logging(args, "fastq_to_uBAM.worker", rename_process=False)
    logger.debug(
        f"starting up with fq1={fq1}, fq2={fq2} sam_out={sam_out} and args={args} _extra_args={_extra_args}"
    )

    def iter_paired(fq1, fq2):
        fq_src1 = util.FASTQ_src(fq1)
        fq_src2 = util.FASTQ_src(fq2)
        # if args.min_qual_trim:
        #     fq_src2 = quality_trim(
        #         fq_src2, min_qual=args.min_qual_trim, phred_base=args.phred_base
        #     )

        return zip(fq_src1, fq_src2)

    def iter_single(fq2):
        fq_src2 = util.FASTQ_src(fq2)
        # if args.min_qual_trim:
        #     fq_src2 = quality_trim(
        #         fq_src2, min_qual=args.min_qual_trim, phred_base=args.phred_base
        #     )

        for fqid, seq, qual in fq_src2:
            # if args.qual_quantization:
            #     qual = quantize_quality(qual)
            yield (fqid, "NA", "NA"), (fqid, seq, qual)

    if fq1:
        ingress = iter_paired(fq1, fq2)
    else:
        ingress = iter_single(fq2)

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

    pre = PreProcessor(args.processing)

    N = mf.util.CountDict()
    c = args.chunk_size
    p = args.parallel

    fmt = make_formatter_from_args(args)  # , **params

    for (fqid, r1, q1), (_, r2, q2) in ingress:
        if pre:
            fqid, r2, q2, preprocess_tags = pre.process(fqid, r2, q2)

        # if args.qual_quantization:
        #     q2 = quantize_quality(q2)

        # if args.simplify_read_id:
        #     i = N.stats["total"]
        #     n = c * (w + (i // c) * (p + w)) + i % c
        #     fqid = simplify_qname(fqid, n)

        N.count("total")
        attrs = fmt(r2_qname=fqid, r1=r1, r1_qual=q1, r2=r2, r2_qual=q2)
        sam_out.write(
            make_sam_record(
                flag=4,
                tags=preprocess_tags
                + [("CB", "{cell}"), ("MI", "{UMI}"), ("CR", "{raw}"), ("RG", "A")],
                **attrs,
            )
        )

    return N


def main(args):
    input_reads1, input_reads2, input_params = get_input_params(args)
    have_read1 = set([str(r1) != "None" for r1 in input_reads1]) == set([True])

    # queues for communication between processes
    w = (
        mf.Workflow(
            f"[{args.sample}]fastq_to_uBAM", total_pipe_buffer_MB=args.pipe_buffer
        )
        # open reads2.fastq.gz
        .gz_reader(inputs=input_reads2, output=mf.FIFO("read2", "wb")).distribute(
            input=mf.FIFO("read2", "rt"),
            outputs=mf.FIFO("r2_{n}", "wt", n=args.parallel),
            chunk_size=args.chunk_size * 4,
        )
    )

    if have_read1:
        # open reads1.fastq.gz
        w.gz_reader(inputs=input_reads1, output=mf.FIFO("read1", "wb"))
        w.distribute(
            input=mf.FIFO("read1", "rt"),
            outputs=mf.FIFO("r1_{n}", "wt", n=args.parallel),
            chunk_size=args.chunk_size * 4,
        )
        # process in parallel workers
        w.workers(
            func=render_to_sam,
            fq1=mf.FIFO("r1_{n}", "rt"),
            fq2=mf.FIFO("r2_{n}", "rt"),
            sam_out=mf.FIFO("sam_{n}", "wt"),
            args=args,
            n=args.parallel,
        )
    else:
        # process in parallel workers
        w.workers(
            func=render_to_sam,
            fq1=None,
            fq2=mf.FIFO("r2_{n}", "rt"),
            sam_out=mf.FIFO("sam_{n}", "wt"),
            args=args,
            n=args.parallel,
        )

    # combine output streams
    w.collect(
        inputs=mf.FIFO("sam_{n}", "rt", n=args.parallel),
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
    # compress to BAM
    fmt_opt = " ".join([f"--output-fmt-option {o}" for o in args.output_fmt_option])
    fmt = f"Sh -O {args.output_fmt} {fmt_opt}"

    w.funnel(
        func=mf.parts.bam_writer,  # mf.parts.null_writer, #
        input=mf.FIFO("sam_combined", "rt"),
        output=args.out_bam,
        _manage_fifos=False,
        fmt=fmt,
        threads=16,
    )
    return w.run()


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
    # input options
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

    # output options
    parser.add_argument(
        "--out-bam",
        default="/dev/stdout",
        help="output for unaligned BAM records (default=/dev/stdout) ",
    )
    parser.add_argument(
        "--save-stats",
        default="preprocessing_stats.txt",
        help="store statistics in this file",
    )
    parser.add_argument(
        "--save-cell-barcodes",
        default="",
        help="store (raw) cell barcode counts in this file. Numbers add up to number of total raw reads. Warning: may use excessive RAM for Open-ST data!",
    )

    # pre-processing options
    parser.add_argument(
        "--processing",
        help="string encoding the processing of the cDNA",
        default="Q:right=25;polyA;adapter:name=SMART,seq=AAGCAGTGGTATCAACGCAGAGTGAATGGG,max_errors=0.1,min_overlap=10",  # ;barcode:cell='r1[8:20][::-1]',UMI='r1[0:8]'
    )

    # filtering options
    # parser.add_argument(
    #     "--min-qual-trim",
    #     default=0,
    #     type=int,
    #     help="clip low quality-bases from a read2 3' end (e.g. pre-detected adapter)",
    # )
    parser.add_argument(
        "--min-len",
        default=18,
        type=int,
        help="minimum read2 length to keep after preprocessing trimming (default=18)",
    )

    # Barcode/UMI extraction options
    parser.add_argument("--cell", default="r1[8:20][::-1]")
    parser.add_argument("--UMI", default="r1[0:8]")
    parser.add_argument("--seq", default="r2")
    parser.add_argument("--qual", default="r2_qual")
    parser.add_argument("--disable-safety", default=False, type=bool)
    # parser.add_argument(
    #     "--fq-qual",
    #     default="E",
    #     help="phred qual for assigned barcode bases in FASTQ output (default='E')",
    # )

    # experimental options to reduce disk footprint
    parser.add_argument(
        "--qual-quantization",
        default=False,
        # type=bool,
        action="store_true",
        help="[EXPERIMENTAL] Quantize quality scores to just 4 levels for better compression/smaller file size (default=False)",
    )
    parser.add_argument(
        "--qname-simplification",
        default=False,
        action="store_true",
        help="[EXPERIMENTAL] Truncate flow-cell tile coordinates from QNAME and replace with counter for better compression/smaller file size (default=False)",
    )

    # parallelization
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
        "--parallel", default=1, type=int, help="how many processes to spawn"
    )
    parser.add_argument(
        "--threads-read",
        help="number of threads for reading bam_in (default=2)",
        type=int,
        default=2,
    )
    parser.add_argument(
        "--threads-write",
        help="number of threads for writing bam_out (default=4)",
        type=int,
        default=4,
    )
    parser.add_argument(
        "--threads-work",
        help="number of worker threads for actual trimming (default=8)",
        type=int,
        default=8,
    )

    # SAM standard now supports CB=corrected cell barcode, CR=original cell barcode, and MI=molecule identifier/UMI
    parser.add_argument(
        "--bam-tags",
        default="CR:{cell},CB:{cell},MI:{UMI},RG:A",
        help="a template of comma-separated BAM tags to generate. Variables are replaced with extracted cell barcode, UMI etc.",
    )
    parser.add_argument(
        "--output-fmt",
        default="BAM",
        choices=[
            "SAM",
            "BAM",
            "CRAM",
        ],
        help="SAM, BAM (default), or CRAM",
    )
    parser.add_argument(
        "--output-fmt-option",
        default="",
        nargs="+",
        help="passed options to `samtools view` as --output-fmt-options",
    )

    args = parser.parse_args()
    if args.parallel < 1:
        raise ValueError(f"--parallel {args.parallel} is invalid. Must be >= 1")

    return args


def cmdline():
    args = parse_args()
    util.setup_logging(args, name="spacemake.bin.fastq_to_uBAM")

    if not args.read2:
        raise ValueError("bam output requires --read2 parameter")

    return main(args)


if __name__ == "__main__":
    res = cmdline()
