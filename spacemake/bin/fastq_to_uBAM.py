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


def quality_trim(fq_src, min_qual=20, phred_base=33):
    for name, seq, qual in fq_src:
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


def render_to_sam(fq1, fq2, sam_out, args, **kwargs):
    logger = util.setup_logging(args, "spacemake.bin.fastq_to_uBAM.worker")
    logger.debug(
        f"starting up with fq1={fq1}, fq2={fq2} sam_out={sam_out} and args={args}"
    )

    def iter_paired(fq1, fq2):
        fq_src1 = util.FASTQ_src(fq1)
        fq_src2 = util.FASTQ_src(fq2)
        if args.min_qual_trim:
            fq_src2 = quality_trim(
                fq_src2, min_qual=args.min_qual_trim, phred_base=args.phred_base
            )

        return zip(fq_src1, fq_src2)

    def iter_single(fq2):
        fq_src2 = util.FASTQ_src(fq2)
        if args.min_qual_trim:
            fq_src2 = quality_trim(
                fq_src2, min_qual=args.min_qual_trim, phred_base=args.phred_base
            )

        for fqid, seq, qual in fq_src2:
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

    N = mf.util.CountDict()

    fmt = make_formatter_from_args(args)  # , **params

    for (fqid, r1, q1), (_, r2, q2) in ingress:
        N.count("total")
        attrs = fmt(r2_qname=fqid, r1=r1, r1_qual=q1, r2=r2, r2_qual=q2)
        sam_out.write(make_sam_record(flag=4, **attrs))

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
    )
    # compress to BAM
    w.funnel(
        func=mf.parts.bam_writer,  # mf.parts.null_writer, #
        input=mf.FIFO("sam_combined", "rt"),
        output=args.out_bam,
        _manage_fifos=False,
        fmt="Sbh",
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
    parser.add_argument("--disable-safety", default=False, type=bool)

    parser.add_argument(
        "--paired-end",
        default=None,
        choices=["fr", "ff", "rf", "rr"],
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
        "--pipe-buffer",
        default=4,
        type=int,
        help="How many megabytes of pipe-buffer to use. kernel settings is usually 64MB per user (default=4MB)",
    )

    parser.add_argument(
        "--parallel", default=1, type=int, help="how many processes to spawn"
    )
    parser.add_argument(
        "--chunk-size",
        default=10,
        type=int,
        help="how many consecutive reads are assigned to the same worker (default=10)",
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

    if not args.read2:
        raise ValueError("bam output requires --read2 parameter")

    return main(args)


if __name__ == "__main__":
    res = cmdline()
