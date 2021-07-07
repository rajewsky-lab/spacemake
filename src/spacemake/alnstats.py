__version__ = "0.9.6"
__author__ = ["Marvin Jens"]
__license__ = "MIT"
__email__ = ["marvin.jens@mdc-berlin.de"]

import os
import sys
import logging
import pandas as pd
import numpy as np
import pysam
import argparse
from collections import defaultdict

import spacemake.reporting as rep

# TODO:
# * sampling of reads
# * graphical output of most common soft-clipped sequences
#   (MUSCLE?)


def gf_prio(gf, hierarchy=["CODING", "UTR", "INTRONIC", "INTERGENIC"]):
    parts = gf.split(",")
    if len(parts) == 1:
        # fast path if only one gf is given
        return gf

    parts = set(parts)
    for h in hierarchy:
        if h in parts:
            return h

    # could not resolve??
    return gf


def parse_args():
    parser = argparse.ArgumentParser("alnstats")
    parser.add_argument("fname", help="a gene-annotated BAM file")
    parser.add_argument(
        "--n-max",
        type=int,
        default=0,
        help="number of BAM records to process (default=0 -> all)",
    )
    parser.add_argument(
        "--out-csv", default="./", help="where to store the CSV raw stats"
    )
    parser.add_argument(
        "--out-pdf", default="./", help="where to store the PDF reports"
    )
    parser.add_argument(
        "--out-png", default="./", help="where to store the PNG reports"
    )
    args = parser.parse_args()

    return args


class Results:
    pass


def scan_bam(fname, n_max=0):
    res = Results()
    res.fname = fname
    res.n_max = n_max
    res.aln_types = defaultdict(int)
    res.cigar_types = defaultdict(int)
    res.tag_types = defaultdict(int)
    res.match_len_by_cigtype = defaultdict(lambda: defaultdict(int))
    res.match_len_by_tag = defaultdict(lambda: defaultdict(int))
    res.sc_lengths = defaultdict(lambda: defaultdict(int))
    res.sc_seq = defaultdict(lambda: defaultdict(int))

    op_dict = list("MIDNSHP=XB")
    sam = pysam.Samfile(fname, "rb")
    for r in sam.fetch(until_eof=True):
        res.aln_types["N_reads"] += 1
        if r.is_unmapped:
            res.aln_types["unmapped"] += 1
            continue

        if r.mapping_quality < 255:
            res.aln_types["multimapper"] += 1
        else:
            res.aln_types["unique"] += 1

        minus_strand = r.is_reverse
        if r.has_tag("gf"):
            tag = gf_prio(r.get_tag("gf"))
        else:
            tag = "INTERGENIC"

        res.tag_types[tag] += 1

        cigar = list(r.cigartuples)
        cigtype = "".join([op_dict[c[0]] for c in cigar])
        cigtype = cigtype.replace("MNM", "M").replace("MNM", "M")  # ignore splicing
        res.cigar_types[cigtype] += 1

        if cigar[0][0] == pysam.CSOFT_CLIP:
            n = cigar[0][1]
            seq = r.query[:n]
            if minus_strand:
                res.sc_seq["3'"][seq] += 1
                res.sc_lengths["3'"][n] += 1
            else:
                res.sc_seq["5'"][seq] += 1
                res.sc_lengths["5'"][n] += 1

        if cigar[-1][0] == pysam.CSOFT_CLIP:
            n = cigar[-1][1]
            seq = r.query[-n:]
            if minus_strand:
                res.sc_seq["5'"][seq] += 1
                res.sc_lengths["5'"][n] += 1
            else:
                res.sc_seq["3'"][seq] += 1
                res.sc_lengths["3'"][n] += 1

        n_match = 0
        for op, n in r.cigartuples:
            if op == pysam.CMATCH:
                n_match += n

        res.match_len_by_cigtype[cigtype][n_match] += 1
        res.match_len_by_cigtype["all"][n_match] += 1
        res.match_len_by_tag[tag][n_match] += 1
        res.match_len_by_tag["all"][n_match] += 1

        if (n_max > 0) and (res.aln_types["N_reads"] >= n_max):
            break

    res.N_total = res.aln_types["N_reads"]

    return res


def make_plots(res):
    import matplotlib.pyplot as plt

    # preparing the plot
    fig, axes = plt.subplots(
        3, 3, figsize=(14, 11)
    )  # gridspec_kw={'height_ratios': [1, 1, 1]})
    axes = np.array(axes).flatten()

    # preparing count dictionaries for the donut plots (pie charts)
    cig_types_d, _ = rep.count_dict_collapse_misc(res.cigar_types, total=res.N_total)
    tag_types_d, _ = rep.count_dict_collapse_misc(res.tag_types, total=res.N_total)
    aln_types_d, _ = rep.count_dict_collapse_misc(
        res.aln_types, total=res.N_total, add_up="N_reads"
    )

    # First row: donut plots
    ciglabels, cigcolors = rep.donut_plot(
        axes[0], cig_types_d, title="CIGAR types", cmap="Accent"
    )

    taglabels, tagcolors = rep.donut_plot(
        axes[1], tag_types_d, title="gf types", cmap="Dark2"
    )

    alnlabels, alncolors = rep.donut_plot(
        axes[2],
        aln_types_d,
        title="Alignment status",
        labels=["unique", "multimapper", "unmapped"],
        cmap="Set2",
    )

    # Second row: histograms of no. matching bases
    rep.len_plot(
        axes[3],
        res.match_len_by_cigtype,
        labels=ciglabels + ["all"],
        colors=cigcolors + ["k"],
        title="CIGAR types",
    )

    rep.len_plot(
        axes[4],
        res.match_len_by_tag,
        labels=taglabels + ["all"],
        colors=tagcolors + ["k"],
        title="gf (annotation) types",
    )

    rep.len_plot(
        axes[5],
        res.sc_lengths,
        title="read position",
        xlabel="soft-clipped bases",
        colors=["r", "b"],
    )

    # Third row: cumulative histograms of no. matching bases
    rep.len_plot(
        axes[6],
        res.match_len_by_cigtype,
        labels=ciglabels + ["all"],
        colors=cigcolors + ["k"],
        title="CIGAR types",
        cumulative=True,
        ylabel="cumulative fraction",
        legend=False,
    )

    rep.len_plot(
        axes[7],
        res.match_len_by_tag,
        labels=taglabels + ["all"],
        colors=tagcolors + ["k"],
        title="gf (annotation) types",
        cumulative=True,
        ylabel="cumulative fraction",
        legend=False,
    )

    rep.len_plot(
        axes[8],
        res.sc_lengths,
        title="read position",
        colors=["r", "b"],
        xlabel="soft-clipped bases",
        cumulative=True,
        ylabel="cumulative fraction",
        legend=False,
    )
    return fig


def cmdline():
    args = parse_args()
    path = args.fname.split("/")

    if len(path) >= 4 and path[-4].startswith("sts"):
        sample_name = path[-4]
    else:
        sample_name = os.path.basename(args.fname).rsplit(".", maxsplit=1)[0]

    logging.info(f"Starting analysis of sample_name='{sample_name}'")
    res = scan_bam(args.fname, n_max=args.n_max)

    import matplotlib

    matplotlib.use("Agg")

    fig = make_plots(res)
    fig.suptitle(sample_name)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(f"{args.out_pdf}/{sample_name}.pdf")
    fig.savefig(f"{args.out_png}/{sample_name}.png", dpi=300)

    # Done with the figure. Text output.
    # TODO: store as CSV instead

    rep.count_dict_out(
        res.aln_types, "alignment types", total=res.N_total, add_up="N_reads"
    )
    rep.count_dict_out(res.cigar_types, "CIGAR types", total=res.N_total)
    rep.count_dict_out(res.tag_types, "annotation (gf) types", total=res.N_total)

    # report the most frequent soft-clipped sequences
    def most_frequent(data, total, n=20, title="top"):
        print(f"### {title}")
        for k, v in sorted(data.items(), key=lambda x: -x[1])[:n]:
            # print(k, v, total)
            print(f"{k}\t{v}\t{100.0 * v / total:.3f} %")

    most_frequent(
        res.sc_seq["5'"],
        np.array(list(res.sc_seq["5'"].values())).sum(),
        title="most frequent 5' soft-clipped sequences",
    )
    most_frequent(
        res.sc_seq["3'"],
        np.array(list(res.sc_seq["3'"].values())).sum(),
        title="most frequent 3' soft-clipped sequences",
    )


if __name__ == "__main__":
    cmdline()
