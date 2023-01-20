import argparse
import spacemake.reporting as rep
import spacemake.util as util
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from spacemake.contrib import __license__, __version__, __author__, __email__
import pandas as pd
import numpy as np
import logging
import os
from collections import defaultdict, OrderedDict


def parse_args():
    parser = util.make_minimal_parser(prog=__file__)
    # parser.add_argument("--sample", type=str, default="a_sample", help="sample_name")
    parser.add_argument(
        "--raw",
        required=True,
        help="readlen statistics file for the raw, unaligned BAM",
    )
    parser.add_argument(
        "--trimmed",
        required=True,
        default="",
        help="stats output of the cutadapt_bam script",
    )
    parser.add_argument(
        "--map-strategy",
        required=True,
        default="",
        help="map-strategy used for this sample",
    )
    parser.add_argument(
        "--mapped",
        type=str,
        default=[],
        nargs="+",
        help="readlen statistics of a mapped BAM file",
    )
    parser.add_argument(
        "--not-mapped",
        type=str,
        default=[],
        nargs="+",
        help="readlen statistics of a not_x unmapped BAM file",
    )
    parser.add_argument(
        "--out-tsv",
        type=str,
        default="overview.tsv",
        help="aggregate table with sample statistics",
    )
    parser.add_argument(
        "--out-pdf",
        type=str,
        default="overview.pdf",
        help="aggregate PDF plot with sample statistics",
    )
    args = parser.parse_args()

    return args


def load_readlen_stats(fname):
    df = pd.read_csv(fname, sep="\t")
    return df


def load_cutadapt(fname):
    df = pd.read_csv(fname, sep="\t")
    reads = df.loc["reads"].set_index("key")
    bases = df.loc["bases"].set_index("key")
    rlens = df.loc["L_final"].set_index("key")

    def name_mangle(names):
        out = []
        for name in names:
            o = ""
            if name.startswith("N_kept"):
                o = "kept"
            elif name.startswith("N_too_short_after"):
                o = (
                    name.replace("N_too_short_after_", "")
                    .replace("left_", "")
                    .replace("right_", "")
                )

            out.append(o)

        return out

    reads["label"] = name_mangle(reads.index)
    reads = reads.query("label != ''")
    reads = reads.set_index("label")

    return reads, bases, rlens


def parse_map_strategy(map_strategy):
    tokens = map_strategy.replace("->", ",").split(",")
    targets = []

    for tok in tokens:
        ref, mapper = tok.split(":")[:2]
        targets.append(f"{mapper}.{ref}")

    return targets


def stacked_bars(ax, df, columns, total, scale=1e6):
    x = np.arange(len(columns))
    bottom = np.zeros(len(columns))
    for label in df.index:
        y = df[columns].loc[label].values
        b = ax.bar(x, y / scale, bottom=bottom, align="center", label=label)
        percent = np.array([100.0 * v / total for v in y])
        blabels = np.array([f"({p:.1f}%)" for p in percent])
        blabels = np.where(percent >= 5, blabels, "")
        ax.bar_label(b, labels=blabels, label_type="center")
        bottom += y / scale

    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)

    ax.legend()
    ax.set_xticks(x)
    ax.set_xticklabels(columns)


def length_hist(ax, df, ref_name, mapped=[False, True]):
    ax.set_title(ref_name)
    ax.set_xlabel("read length [nt]")
    # w = 1./(len(mapped) + 1)

    for i, m in enumerate(mapped):
        _df = (
            df.query(f"mapped == {m} and ref_name == '{ref_name}'")
            .reset_index()
            .sort_values("readlen")
        )
        x = _df["readlen"]
        y = _df["count"]
        # print("plotting", _df)
        if len(mapped) > 1:
            label = "aligned" if m else "not aligned"
        else:
            label = "trimmed"

        ax.bar(x, y, width=1, align="center", label=label, alpha=0.5)

    ax.legend(loc="upper left", frameon=False)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)


def main(args):
    raw_rlen = load_readlen_stats(args.raw)
    N_raw = raw_rlen["count"].sum()

    trimmed_reads, trimmed_bases, trimmed_rlens = load_cutadapt(args.trimmed)
    # print(trimmed_reads)
    # print(trimmed_bases)
    # print(trimmed_rlens)

    targets = parse_map_strategy(args.map_strategy)
    # print(targets)

    rlen_dfs = []
    data = defaultdict(list)
    # print(args.mapped)
    # print(trimmed_rlens)
    df_trim = trimmed_rlens.reset_index()
    df_trim["ref_name"] = "adapter_trimmed"
    df_trim["mapped"] = False
    df_trim["readlen"] = df_trim["key"].astype(int)
    rlen_dfs.append(df_trim[["ref_name", "mapped", "readlen", "count"]])

    for m in args.not_mapped:
        name = os.path.basename(m).replace(".readlen.tsv", "")
        # print(m, name)
        df = load_readlen_stats(m)
        ref_name = name.replace("not_", "")
        data["name"].append(ref_name)
        data["count"].append(df["count"].sum())
        data["status"].append("not_aligned")

        df["ref_name"] = ref_name
        df["mapped"] = False
        rlen_dfs.append(df[["ref_name", "mapped", "readlen", "count"]])

    for m in args.mapped:
        name = os.path.basename(m).replace(".readlen.tsv", "")
        # print(m, name)
        df = load_readlen_stats(m)
        data["name"].append(name)
        data["count"].append(df["count"].sum())
        data["status"].append("aligned")

        ref_name = name.replace("not_", "")
        df["ref_name"] = ref_name
        df["mapped"] = True
        rlen_dfs.append(df[["ref_name", "mapped", "readlen", "count"]])

    df_mapped = (
        pd.DataFrame(data).set_index("status").pivot(columns="name", values="count")
    )
    # print(df_mapped)
    df_mapped.to_csv(args.out_tsv, sep="\t")

    # mapstats = df_mapped[["count", "bamname"]].groupby("bamname").aggregate('sum')

    # aligned = []
    # not_aligned = []

    # for t in targets:
    #     aligned.append(mapstats.loc[t, 'count'])
    #     not_aligned.append(mapstats.loc[f"not_{t}", 'count'])

    # df = pd.DataFrame()

    # plot read count statistics
    n_mapped = len(args.mapped)
    n_col = 1 + n_mapped
    fig, axes = plt.subplots(
        1, 2, figsize=(3 * n_col, 5), width_ratios=[1, n_mapped], sharey=True
    )
    fig.suptitle(args.sample)
    stacked_bars(
        axes[0],
        trimmed_reads,
        columns=[
            "count",
        ],
        total=N_raw,
        scale=1e6,
    )
    axes[0].set_ylabel("raw reads [M]")
    axes[0].set_title("adapter trimming")
    stacked_bars(
        axes[1],
        df_mapped.loc[["not_aligned", "aligned"]],
        columns=targets,
        total=N_raw,
        scale=1e6,
    )
    axes[1].set_title("read mapping")
    plt.savefig(args.out_pdf)

    # plot read-length histograms
    df_readlen = pd.concat(rlen_dfs).query("readlen > 0").set_index("readlen")
    # print(df_readlen)

    df_readlen.to_csv(args.out_tsv.replace(".tsv", ".readlen.tsv"), sep="\t")
    fig, axes = plt.subplots(1, n_col, figsize=(3 * n_col, 5), sharex=True)
    fig.suptitle(args.sample)
    length_hist(axes[0], df_readlen, ref_name="adapter_trimmed", mapped=[False])
    for ref_name, ax in zip(targets, axes[1:]):
        length_hist(ax, df_readlen, ref_name=ref_name)

    fig.tight_layout()
    fig.savefig(args.out_pdf.replace(".pdf", ".readlen.pdf"))

    return df_mapped


if __name__ == "__main__":
    args = parse_args()
    util.setup_logging(args, "spacemake.bin.plot_overview")
    df_mapped = main(args)
