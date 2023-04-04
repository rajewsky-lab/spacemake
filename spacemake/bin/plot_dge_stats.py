import scanpy as sc
import matplotlib as mpl
import os

mpl.use("Agg")

import argparse
import logging
import numpy as np

import spacemake.util as util
from spacemake.pl import plot_dge_stats, plot_loglog_knee

def parse_args():
    parser = util.make_minimal_parser(prog="plot_dge_stats.py")
    parser.add_argument("dge", help="path to h5ad AnnData file")
    parser.add_argument(
        "--reference",
        default="",
        help="subset analysis only to genes from this mapping reference",
    )

    parser.add_argument(
        "--out-knee",
        default="dge_stats/knee.pdf",
        help="destination path for knee plot",
    )
    parser.add_argument(
        "--out-metrics",
        default="dge_stats/metrics.pdf",
        help="destination path for metrics plots",
    )

    args = parser.parse_args()
    return args

def main(args):
    adata = sc.read_h5ad(args.dge)

    # subsetting to miRNA only counts
    if args.reference:
        mask = adata.var["reference"] == args.reference
        adata = adata[:, mask]

    fig = plot_dge_stats(adata)
    fig.savefig(util.ensure_path(args.out_metrics))


    # loglog version of knee plot
    fig = plot_loglog_knee(adata, key="n_counts")
    fig.savefig(util.ensure_path(args.out_knee))

    return adata


if __name__ == "__main__":
    args = parse_args()
    util.setup_logging(args, "spacemake.bin.plot_dge_stats")
    adata = main(args)
