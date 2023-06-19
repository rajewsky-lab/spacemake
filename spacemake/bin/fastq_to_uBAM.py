#!/usr/bin/env python3
from spacemake.contrib import __version__, __author__, __license__, __email__
import logging
import os
import sys
import numpy as np
import pysam
import multiprocessing as mp
from collections import defaultdict
from time import time
from spacemake.parallel import (
    put_or_abort,
    chunkify,
    queue_iter,
    iter_queues_round_robin,
    join_with_empty_queues,
    ExceptionLogging,
    log_qerr,
)
import spacemake.util as util

# TODO:
# * load params from YAML (code from opseq)
# * save params to YAML in run-folder to document
# * refactor into multiple modules?


def read_to_queues(input_files, params, Qfq, args, Qerr, abort_flag):
    """
    Reads from a fastq (.gz) file, gathers the raw FASTQ data
    into chunks, and places them - round-robin - onto Queues
    for the worker processes.

    Args:
        input_files (list): List of input files to be read
        params (list): List of parameter dictionaries to supersede args for each input file
        Qfq (_type_): list of n_workers Queues to place chunks on
        args (_type_): parsed cmdline args
        Qerr (_type_): Queue connecting to main process for error message retrieval
        abort_flag (_type_): Shared variable that indicates shutdown due to an
                             unhandled exception somewhere
    """
    import gzip

    chunk = []
    n_chunk = 0
    n_reads = 0
    logger = util.setup_logging(args, "spacemake.bin.fastq_to_uBAM.read_to_queues")
    T0 = time()
    dT = 0

    def dispatch(par):
        nonlocal chunk, n_chunk, n_reads, T0
        # send this chunk to worker i
        # using round-robin
        i = n_chunk % args.parallel
        logger.debug(f"putting chunk {n_chunk} onto queue {i}")
        Qfq[i].put((n_chunk, chunk, par))
        n_chunk += 1
        n_reads += len(chunk) / 4
        chunk = []

        dT = time() - T0
        if dT > 10:
            logger.info(
                f"ingested {n_reads/1000:.1f}k reads in {dT:.1f} seconds ({n_reads/dT/1000:.2f} reads/s)."
            )
            T0 = time()
            n_reads = 0

    with ExceptionLogging(
        f"spacemake.fastq_to_uBAM.read_to_queues {args.sample}",
        Qerr=Qerr,
        exc_flag=abort_flag,
    ) as el:
        for fname, par in zip(input_files, params):
            logger.info(f"iterating over reads from '{fname}' with special params={par}")
            if fname.endswith(".gz"):
                src = gzip.open(fname, mode="rt")
            else:
                src = open(fname)

            for line in src:
                chunk.append(line)
                if len(chunk) >= args.chunk_size:
                    dispatch(par)

            if len(chunk):
                dispatch(par)

        el.logger.info("finished! Closing down.")
        for i in range(args.parallel):
            Qfq[i].put(None)  # each worker consumes exactly one None


format_func_template = """
def format_func(qname=None, r2_qname=None, r2_qual=None, r1=None, r2=None, R1=None, R2=None):
    qparts = r2_qname.split(':N:0:')
    if len(qparts) > 1:
        i5i7 = qparts[1].replace('+', '')
    else:
        i5i7 = None

    cell = {cell}

    attrs = dict(
        cell=cell,
        raw=cell,
        UMI={UMI},
        seq={seq},
        qual={qual}
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


def make_BAM_header(args):
    prog = os.path.basename(__file__)
    header = {
        "HD": {"VN": "1.6"},
        "PG": [
            {
                "ID": "fastq_to_uBAM",
                "PN": prog,
                "CL": " ".join(sys.argv[1:]),
                "VN": __version__,
            },
        ],
        "RG": [
            {"ID": "A", "SM": args.sample},
        ],
    }
    return pysam.AlignmentHeader.from_dict(header)


class Formatter:
    logger = logging.getLogger("spacemake.fastq_to_uBAM.Formatter")

    def __init__(self, args):
        # assert Output.safety_check_eval(args.cell_raw)

        self.bam_header = make_BAM_header(args)
        self.format = make_formatter_from_args(args)

    def make_bam_record(
        self,
        flag=4,
        tags=[("CB", "{cell}"), ("MI", "{UMI}"), ("CR", "{raw}"), ("RG", "A")],
        **kw,
    ):

        # pass what we already know about the read into format()
        # which contains the cmdline-specified code for cell and UMI.
        # since it returns a dictionary, lets just update kw and
        # keep using that to construct the BAM record.
        try:
            kw.update(self.format(**kw))
        except Exception as E:
            self.logger.error(
                "Unhandled exception inside user-defined read handler.\n"
                f"Input: {kw}\n"
                f"Error: {E}"
            )
            raise E

        a = pysam.AlignedSegment(self.bam_header)

        # STAR does not like spaces in read names so we have to split
        a.query_name = kw["r2_qname"].split()[0]
        a.query_sequence = kw["seq"]
        a.flag = flag
        a.query_qualities = pysam.qualitystring_to_array(kw["qual"])
        a.tags = [(name, templ.format(**kw)) for name, templ in tags]
        # if self.count_cb:
        #     self.raw_cb_counts[kw["cell"]] += 1

        return a.to_string()


def quality_trim(fq_src, min_qual=20, phred_base=33):
    for (name, seq, qual) in fq_src:
        end = len(seq)
        q = np.array(bytearray(qual.encode("ASCII"))) - phred_base
        qtrim = q >= min_qual
        new_end = end - (qtrim[::-1]).argmax()

        # TODO: yield A3,T3 adapter-trimming tags
        # TODO: convert to cutadapt/BWA qual-trim logic
        if new_end != end:
            qual = qual[:new_end]
            seq = seq[:new_end]

        yield (name, seq, qual)


def parallel_worker(Qfq1, Qfq2, Qres, args, Qerr, abort_flag, stat_lists):
    logger = util.setup_logging(args, "spacemake.bin.fastq_to_uBAM.worker")

    def iter_paired(Qfq1, Qfq2):
        src1 = queue_iter(Qfq1, abort_flag)
        src2 = queue_iter(Qfq2, abort_flag)
        for (n1, chunk1, par1), (n2, chunk2, par2) in zip(src1, src2):
            logger.debug(
                f"received chunk {n1} {n2} of paired {len(chunk1)} {len(chunk2)} raw FASTQ lines"
            )
            assert n1 == n2
            assert par1 == par2

            fq_src1 = util.FASTQ_src(chunk1)
            fq_src2 = util.FASTQ_src(chunk2)
            if args.min_qual_trim:
                fq_src2 = quality_trim(
                    fq_src2, min_qual=args.min_qual_trim, phred_base=args.phred_base
                )

            out = []
            for (_, seq1, qual1), (name2, seq2, qual2) in zip(fq_src1, fq_src2):
                out.append((name2, seq1, qual1, seq2, qual2))

            yield n1, out, par1

    def iter_single(Qfq2):
        for n, chunk, par in queue_iter(Qfq2, abort_flag):
            logger.debug(
                f"received chunk {n} of single-end {len(chunk)} raw FASTQ lines"
            )

            fq_src2 = util.FASTQ_src(chunk)
            if args.min_qual_trim:
                fq_src2 = quality_trim(
                    fq_src2, min_qual=args.min_qual_trim, phred_base=args.phred_base
                )

            out = [(name2, "NA", "NA", seq2, qual2) for name2, seq2, qual2 in fq_src2]
            yield n, out, par

    with ExceptionLogging(
        f"spacemake.bin.fastq_to_uBAM.worker {args.sample}",
        Qerr=Qerr,
        exc_flag=abort_flag,
    ) as el:
        el.logger.debug(
            f"starting up with Qfq1={Qfq1}, Qfq2={Qfq2} Qres={Qres} and args={args}"
        )

        N = defaultdict(int)

        if Qfq1:
            ingress = iter_paired(Qfq1, Qfq2)
        else:
            ingress = iter_single(Qfq2)

        for n_chunk, chunk, par in ingress:
            # TODO: also retrieve parameters such as CB=NNNN for each chunk
            # and recreate a formatter. Since chunk_size is ~100k, this should
            # be negligible but add flexibility for multiple input files etc.
            kw = vars(args)
            kw.update(par)
            import argparse
            args = argparse.Namespace(**kw)

            fmt = Formatter(args)  # , **params
            results = []
            for fqid, r1, q1, r2, q2 in chunk:
                N["total"] += 1
                if args.paired_end: 
                    # if args.paired_end[0] == 'r':
                    #     # rev_comp R1 and reverse qual1
                    #     q1 = q1[::-1]
                    #     r1 = util.rev_comp(r1)

                    # if args.paired_end[1] == 'r':
                    #     # rev_comp R2 and reverse qual2
                    #     q2 = q2[::-1]
                    #     r2 = util.rev_comp(r2)

                    rec1 = fmt.make_bam_record(
                        qname=fqid,
                        r1="THISSHOULDNEVERBEREFERENCED",
                        r2=r1,
                        R1=r1,
                        R2=r2,
                        r2_qual=q1,
                        r2_qname=fqid,
                        flag=69,  # unmapped, paired, first in pair
                    )
                    rec2 = fmt.make_bam_record(
                        qname=fqid,
                        r1="THISSHOULDNEVERBEREFERENCED",
                        r2=r2,
                        R1=r1,
                        R2=r2,
                        r2_qual=q2,
                        r2_qname=fqid,
                        flag=133,  # unmapped, paired, second in pair
                    )
                    results.append(rec1)
                    results.append(rec2)

                else:
                    rec = fmt.make_bam_record(
                        qname=fqid,
                        r1=r1,
                        r2=r2,
                        r2_qual=q2,
                        r2_qname=fqid,
                    )
                    results.append(rec)

            Qres.put((n_chunk, results))

            if Qfq1:
                Qfq1.task_done()
            if Qfq2:
                Qfq2.task_done()

        # return our counts and observations
        # via these list proxies
        el.logger.debug("synchronizing cache and counts")
        Ns = stat_lists[0]
        Ns.append(N)

        # cb_counts = stat_lists[1]
        # cb_counts.append(out.raw_cb_counts)


def collect_from_queues(Qres, args, Qerr, abort_flag):
    with ExceptionLogging(
        f"spacemake.bin.fastq_to_uBAM.collector {args.sample}",
        Qerr=Qerr,
        exc_flag=abort_flag,
    ) as el:
        header = make_BAM_header(args)
        bam = pysam.AlignmentFile(args.out_bam, "wbu", header=header, threads=4)

        for n_chunk, chunk in iter_queues_round_robin(
            Qres, abort_flag, logger=el.logger
        ):
            el.logger.debug(f"received chunk: {n_chunk}")
            for rec in chunk:
                bam.write(pysam.AlignedSegment.fromstring(rec, header))

        bam.close()


def count_dict_sum(sources):
    dst = defaultdict(float)
    for src in sources:
        for k, v in src.items():
            dst[k] += v

    return dst


def dict_merge(sources):
    dst = {}
    for src in sources:
        dst.update(src)
    return dst


def get_input_params(args):
    import pandas as pd
    import csv
    if str(args.matrix) != "None":
        df = pd.read_csv(args.matrix, sep=',', index_col=None, quoting=csv.QUOTE_NONE)
        if not 'R1' in df.columns:
            df['R1'] = 'None'

        R1 = df['R1'].to_list()
        R2 = df['R2'].to_list()
        # TODO extract other column values into a list of dictionaries (most importantly cell=...)
        # which can override how the formatter in the workers processes the raw reads
        if 'cell' in df.columns:
            # params = [{'cell': f'"{c}"'} for c in df['cell']]
            params = [{'cell': c} for c in df['cell']]
        else:
            params = [{},] * len(R1)

    else:
        R1 = [args.read1,]
        R2 = [args.read2,]
        params = [{},]

    return R1, R2, params

def main_parallel(args):
    input_reads1, input_reads2, input_params = get_input_params(args)
    have_read1 = set([str(r1) != "None" for r1 in input_reads1]) == set([True])

    # queues for communication between processes
    if have_read1:
        Qfq1 = [mp.JoinableQueue(5) for i in range(args.parallel)]
    else:
        Qfq1 = [None for i in range(args.parallel)]

    Qfq2 = [mp.JoinableQueue(5) for i in range(args.parallel)]
    Qres = [mp.Queue(5) for i in range(args.parallel)]  # BAM records as string
    Qerr = mp.Queue()  # child-processes can report errors back to the main process here

    # Proxy objects to allow workers to report statistics about the run
    manager = mp.Manager()
    abort_flag = mp.Value("b")
    abort_flag.value = False
    Ns = manager.list()
    cb_counts = manager.list()
    stat_lists = [Ns, cb_counts]
    workers = []

    res = 0
    with ExceptionLogging(
        "spacemake.bin.fastq_to_uBAM.main_parallel", exc_flag=abort_flag
    ) as el:

        # read FASTQ in chunks and put them in Qfq
        if have_read1:
            reader1 = mp.Process(
                target=read_to_queues,
                name="reader1",
                args=(input_reads1, input_params, Qfq1, args, Qerr, abort_flag),
            )
            reader1.start()

        reader2 = mp.Process(
            target=read_to_queues,
            name="reader2",
            args=(input_reads2, input_params, Qfq2, args, Qerr, abort_flag),
        )

        reader2.start()

        el.logger.debug("Started readers")

        # worker i consumes chunks of paired FASTQ from Qfq1[i] and Qfq2[i],
        # processes them, and put the results in Qres[i]
        for i in range(args.parallel):
            w = mp.Process(
                target=parallel_worker,
                name=f"worker_{i}",
                args=(Qfq1[i], Qfq2[i], Qres[i], args, Qerr, abort_flag, stat_lists),
            )
            w.start()
            workers.append(w)

        el.logger.debug("Started workers")
        writer = mp.Process(
            target=collect_from_queues,
            name="writer",
            args=(Qres, args, Qerr, abort_flag),
        )
        writer.start()
        el.logger.debug("Started writer")

        # wait until all sequences have been thrown onto Qfq
        if have_read1:
            q1 = join_with_empty_queues(reader1, Qfq1 + [Qerr], abort_flag)

        q2 = join_with_empty_queues(reader2, Qfq2 + [Qerr], abort_flag)
        el.logger.debug("The readers exited")

        # qfq = qfq1 + qfq2
        # qerr = qerr1 + qerr2
        # if qfq or qerr:
        #     el.logger.error("Some exception occurred! Tearing it all down")
        #     res = 1
        #     el.logger.warning(f"{len(qfq)} chunks were drained from Qfq upon abort.")
        #     log_qerr(qerr)
        #     reader1.terminate()
        #     writer.terminate()
        #     for w in workers:
        #         w.terminate()
        # if False:
        #     pass
        # else:
        # signal all workers to finish
        # el.logger.debug("Signalling all workers to finish")
        # for i in range(args.parallel):
        #     Qfq1[i].put(None)  # each worker consumes exactly one None
        #     Qfq2[i].put(None)  # from each queue

        for w in workers:
            # make sure all results are on Qres by waiting for
            # workers to exit. Or, empty queues if aborting.
            q = join_with_empty_queues(w, Qres + [Qerr], abort_flag)
            # if qres or qerr:
            #     res = 1
            #     el.logger.info(f"{len(qres)} chunks were drained from Qres upon abort.")
            #     log_qerr(qerr)

        el.logger.info("All worker processes have joined. Signalling writer to finish.")
        # signal the collector to stop
        [Q.put(None) for Q in Qres]

        # and wait until all output has been generated
        q = join_with_empty_queues(writer, Qres, abort_flag)
        el.logger.info("Collector has joined. Merging worker statistics.")

        N = count_dict_sum(Ns)
        if N["total"]:
            el.logger.info(f"Run completed, {N['total']} reads processed.")
        else:
            el.logger.error("No reads were processed!")

        cb_count = count_dict_sum(cb_counts)
        if args.save_cell_barcodes:
            el.logger.info(
                f"writing {len(cb_count)} barcode counts to '{args.save_cell_barcodes}'"
            )
            import gzip

            with gzip.open(util.ensure_path(args.save_cell_barcodes), "wt") as f:
                f.write("cell_bc\traw_read_count\n")
                for bc, count in sorted(cb_count.items()):
                    f.write(f"{bc}\t{count}\n")

        if args.save_stats:
            with open(util.ensure_path(args.save_stats), "w") as f:
                for k, v in sorted(N.items()):
                    f.write(f"freq\t{k}\t{v}\t{100.0 * v/max(N['total'], 1):.2f}\n")

    if el.exception:
        res = 1
        for w in workers:
            w.terminate()
            w.join()

    return res


def parse_args():
    import spacemake.util as util

    parser = util.make_minimal_parser(
        "fastq_to_uBAM.py",
        description="Convert raw reads1 and reads2 FASTQ into a single BAM file with cell barcode and UMI as BAM-tags",
    )
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
    parser.add_argument("--cell", default="r1[8:20][::-1]")
    parser.add_argument("--UMI", default="r1[0:8]")
    parser.add_argument("--seq", default="r2")
    parser.add_argument("--qual", default="r2_qual")
    parser.add_argument(
        "--paired-end",
        default=None,
        choices=['fr', 'ff', 'rf', 'rr'],
        help="read1 and read2 are paired end mates and store both in the BAM",
    )
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
        help="store (raw) cell barcode counts in this file. Numbers add up to number of total raw reads.",
    )
    parser.add_argument(
        "--fq-qual",
        default="E",
        help="phred qual for assigned barcode bases in FASTQ output (default='E')",
    )
    parser.add_argument(
        "--min-qual-trim",
        default=0,
        type=int,
        help="clip low quality-bases from a read2 3' end (e.g. pre-detected adapter)",
    )
    parser.add_argument(
        "--min-len",
        default=18,
        type=int,
        help="minimum read2 length to keep after quality trimming (default=18)",
    )
    parser.add_argument(
        "--phred-base",
        default=33,
        type=int,
        help="phred quality base (default=33)",
    )
    parser.add_argument(
        "--parallel", default=1, type=int, help="how many processes to spawn"
    )
    parser.add_argument(
        "--chunk-size",
        default=10000,
        type=int,
        help="how many reads (mate pairs) are grouped together for parallel processing (default=10000)",
    )
    # SAM standard now supports CB=corrected cell barcode, CR=original cell barcode, and MI=molecule identifier/UMI
    parser.add_argument(
        "--bam-tags",
        default="CR:{cell},CB:{cell},MI:{UMI},RG:A",
        help="a template of comma-separated BAM tags to generate. Variables are replaced with extracted cell barcode, UMI etc.",
    )
    args = parser.parse_args()
    if args.parallel < 1:
        raise ValueError(f"--parallel {args.parallel} is invalid. Must be >= 1")

    return args


def cmdline():
    args = parse_args()
    util.setup_logging(args, name="spacemake.bin.fastq_to_uBAM")

    with ExceptionLogging("spacemake.preprocess.main") as el:
        if not args.read2:
            raise ValueError("bam output format requires --read2 parameter")

        res = main_parallel(args)
    if el.exception:
        res = 1

    el.logger.debug(f"exit code={res}")
    return res


if __name__ == "__main__":
    import os, sys

    ## We catch errors in the workers, but somehow
    # multiprocessing masks our non-zero exit code.
    # I think this is a bug...
    # sys.stderr.write("running that shit\n")
    res = cmdline()
    # sys.stderr.write(f"calling sys.exit({res})\n")
    # os._exit(1)
    # sys.exit(1)
    sys.exit(res)
