import scanpy as sc
import matplotlib as mpl
import os

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

import argparse
import logging
import numpy as np

import spacemake.util as util


def cell_hist_plot(ax, df, key="n_reads", n_bins=100, xlog=True, ylog=True, **kw):
    values = df[key].values
    # print(values)
    if xlog:
        lv = np.log10(values + 1)
    else:
        lv = values

    # exclude extreme outliers
    lmin, lmax = np.percentile(lv, [0.1, 99.9])
    lv = lv[(lv >= lmin) & (lv <= lmax)]
    hist, bin_edges = np.histogram(lv, bins=n_bins)
    if xlog:
        bin_edges = 10 ** bin_edges

    w = bin_edges[1:] - bin_edges[:-1]
    ax.bar(bin_edges[:-1], hist, width=w, **kw)

    if xlog:
        ax.set_xscale("log")

    if ylog:
        ax.set_yscale("log")

    ax.set_xlabel(f"{key}")
    ax.set_ylabel("no. cell barcodes")

    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)

    # Only show ticks on the left and bottom spines
    # ax.yaxis.set_ticks_position("left")
    # ax.xaxis.set_ticks_position("bottom")

    for axis in [ax.xaxis, ax.yaxis]:
        formatter = ScalarFormatter()
        formatter.set_scientific(False)
        axis.set_major_formatter(formatter)


def loglog_knee_plot(ax, df, key="n_counts"):
    # UMI = np.sort(np.array(df[key].values, dtype=int))[::-1]
    UMI = np.array(df[key].values, dtype=int)
    if not len(UMI):
        return

    UMIsum = UMI.sum()
    # ax.plot(UMI.values.cumsum()[:n_rank], label=f"{key} (total={UMIsum/1e6:.1f} M)")
    # upper_q = np.percentile(UMI, 95)
    # if idx == "miRNA":
    # axh.hist(UMI, bins=100, label=idx, histtype="step", range=(0, upper_q))
    ax.loglog(
        len(UMI) - np.bincount(UMI).cumsum(),
        label=f"{key} (total={UMIsum/1e6:.1f} M)",
        alpha=0.75,
        lw=3,
    )

    # ax.axvline(100, color="k", ls="dashed", lw=0.1)
    ax.set_ylabel("cumulative no. cells")
    ax.set_xlabel(f"{key} > x")

    # ax.set_ylabel(f"cumulative {key}")
    # ax.set_xlabel("cell BC rank")
    # ax.set_xlim(0, n_rank)
    ax.legend()
    # fig.tight_layout()


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

    sc.pp.filter_cells(adata, min_counts=1)
    # re-generating the metrics
    adata.obs["n_exonic_reads"] = adata.layers["exonic_reads"].sum(axis=1)[:, 0]
    adata.obs["n_exonic_counts"] = adata.X.sum(axis=1)[:, 0]
    adata.obs["n_genes"] = (adata.X > 0).sum(axis=1)[:, 0]
    adata.obs["reads_per_counts"] = (
        adata.obs["n_exonic_reads"] / adata.obs["n_exonic_counts"]
    )

    # print(adata)
    # print(adata.obs)
    # print(adata.var)
    sample_name = adata.uns["sample_name"]

    # cell metrics plot
    fig, axes = plt.subplots(2, 2, figsize=(8, 5))
    axes = np.array(axes).ravel()
    cell_hist_plot(axes[0], adata.obs, key="n_exonic_reads")
    cell_hist_plot(axes[1], adata.obs, key="n_counts")
    cell_hist_plot(axes[2], adata.obs, key="n_genes")
    cell_hist_plot(axes[3], adata.obs, key="reads_per_counts", xlog=False)
    fig.tight_layout()
    fig.savefig(util.ensure_path(args.out_metrics))

    # loglog version of knee plot
    fig, ax = plt.subplots(figsize=(6, 5))
    fig.suptitle(sample_name)
    loglog_knee_plot(ax, adata.obs, key="n_counts")
    ax.grid(axis="y")
    fig.tight_layout()
    fig.savefig(util.ensure_path(args.out_knee))

    return adata


if __name__ == "__main__":
    args = parse_args()
    util.setup_logging(args, "spacemake.bin.plot_dge_stats")
    adata = main(args)
