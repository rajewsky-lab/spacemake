import logging
import pysam
import os
from collections import defaultdict, OrderedDict
import argparse
import pandas as pd
import numpy as np
from spacemake.contrib import __version__, __license__, __author__, __email__
import spacemake.util as util


def parse_cmdline():
    parser = util.make_minimal_parser(
        prog="cell_counter.py",
        description="count how often cell barcode values occur in the provided BAM files"
    )

    parser.add_argument(
        "bam_in",
        help="bam input (default=stdin)",
        default=["/dev/stdin"],
        nargs="+",
    )
    parser.add_argument(
        "--sample-name",
        help="sample identifier (default='NA')",
        default="NA",
    )
    parser.add_argument(
        "--skim",
        help="skim through the BAM(s) by investigating only every <skim>-th record (default=0 off)",
        default=0,
        type=int,
    )
    parser.add_argument(
        "--top",
        help="keep only the top n most observed cell barcodes (default=100000)",
        default=100000,
        type=int,
    )
    parser.add_argument(
        "--tag",
        help="bam tag to investigate (default='CB')",
        default="CB",
    )
    parser.add_argument(
        "--unique",
        help="count each read only once (by qname)",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--UMI",
        help="count UMIs as well. Use this bam tag (default='MI'). set to '' to disable.",
        default="MI",
    )
    parser.add_argument(
        "--unmapped",
        help="count also unmapped reads",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--out",
        help="filename of output-file (default={{args.sample_name}}.top_barcodes.tsv",
        default="{args.sample_name}.top_barcodes.tsv",
    )
    args = parser.parse_args()
    return args


def main(args):
    logger = util.setup_logging(args, "spacemake.cell_counter")
    read_counts = defaultdict(int)
    read_counts_by_ref = OrderedDict()
    umi_counts = defaultdict(set)
    umi_counts_by_ref = {}

    n_aln = 0
    unique = set()
    for bam_name in args.bam_in:
        reference_name = os.path.basename(bam_name).split(".")[0]
        ref_read_counts = defaultdict(int)
        ref_umi_counts = defaultdict(set)
        logger.info(f"processing {bam_name}")
        bam = util.quiet_bam_open(bam_name, check_sq=False)

        for aln in util.timed_loop(bam.fetch(until_eof=True), logger, skim=args.skim):
            n_aln += 1
            if aln.is_unmapped and not args.unmapped:
                continue

            if args.unique:
                if aln.query_name in unique:
                    continue

                unique.add(aln.query_name)

            tags = dict(aln.get_tags())
            cell = tags.get(args.tag, "NA")
            umi = tags.get(args.UMI, "NA")

            read_counts[cell] += 1
            ref_read_counts[cell] += 1

            umi_counts[cell].add(umi)
            ref_umi_counts[cell].add(umi)

        read_counts_by_ref[reference_name] = ref_read_counts
        umi_counts_by_ref[reference_name] = dict(
            [(k, len(v)) for k, v in ref_umi_counts.items()]
        )

    all_cbs = list(read_counts.keys())
    umi_counts = dict([(k, len(v)) for k, v in umi_counts.items()])

    data = dict(
        cell_bc=all_cbs,
        reads=[read_counts[cb] for cb in all_cbs],
        umis=[umi_counts[cb] for cb in all_cbs],
    )
    for ref, ref_counts in read_counts_by_ref.items():
        data[ref + "_reads"] = [ref_counts[cb] for cb in all_cbs]
        data[ref + "_umis"] = [umi_counts_by_ref[ref][cb] for cb in all_cbs]

    df = pd.DataFrame(data).set_index("cell_bc").sort_values("umis", ascending=False)
    if args.top > 0:
        df = df.iloc[: args.top]

    logger.info(
        f"counted {len(all_cbs)} barcodes in {n_aln} total alignments. Writing output"
    )
    df.to_csv(util.ensure_path(args.out.format(args=args)), sep="\t")


if __name__ == "__main__":
    args = parse_cmdline()
    main(args)
