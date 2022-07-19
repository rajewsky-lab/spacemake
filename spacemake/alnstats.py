__version__ = "0.9.6"
__author__ = ["Marvin Jens"]
__license__ = "GPL"
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
        "--skim",
        type=int,
        default=0,
        help="only examin every n-th alignment (default=0 -> all)",
    )
    parser.add_argument(
        "--intact-signature",
        type=str,
        default="",
        help="use for longread cDNA: analyze separately alignments with and without intact-signature in the read-name",
    )
    parser.add_argument(
        "--parse-oligos",
        default=False,
        action="store_true",
        help="use for longread cDNA: calculate association of oligo-presence/absence with map status and gene annotation",
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


def coarsegrain_CIGAR(cigar):
    """
    recursively remove insertions, deletions, and gaps (introns) (I,D,N)
    from CIGAR string. Keep only matching and soft-clipping info (M,S)
    """
    coarse = cigar.replace("MNM", "M").replace("MDM", "M").replace("MIM", "M")
    if coarse == cigar:
        return cigar
    else:
        return coarsegrain_CIGAR(coarse)


def scan_bam(fname, skim=0, intact_signature="", parse_oligos=False):
    res = Results()
    res.fname = fname
    res.skim = skim
    res.aln_types = defaultdict(int)
    res.cigar_types = defaultdict(int)
    res.tag_types = defaultdict(int)
    res.uniq_tag_types = defaultdict(int)
    res.match_len_by_cigtype = defaultdict(lambda: defaultdict(int))
    res.match_len_by_tag = defaultdict(lambda: defaultdict(int))
    res.match_len_by_utag = defaultdict(lambda: defaultdict(int))
    res.sc_lengths = defaultdict(lambda: defaultdict(int))
    res.sc_seq = defaultdict(lambda: defaultdict(int))
    res.read_lengths_by_mapstat = defaultdict(lambda: defaultdict(int))
    res.oligo_by_mapstat = defaultdict(lambda: defaultdict(int))
    res.oligo_by_tag = defaultdict(lambda: defaultdict(int))
    res.oligo_count = defaultdict(int)
    res.all_oligos = set()

    op_dict = list("MIDNSHP=XB")
    sam = pysam.Samfile(fname, "rb")
    for i, r in enumerate(sam.fetch(until_eof=True)):
        if skim:
            if (i % skim) != 0:
                continue

        if intact_signature:
            if not intact_signature in r.qname:
                continue

        parts = []
        if parse_oligos:
            try:
                parts = set(r.qname.split("__sig:")[1].split("__")[0].split(","))
            except (ValueError, KeyError):
                pass

        res.aln_types["N_reads"] += 1
        for p in parts:
            res.oligo_count[p] += 1
            res.all_oligos.add(p)

        l_read = len(r.query_sequence)
        res.read_lengths_by_mapstat["all"][l_read] += 1

        if r.is_unmapped:
            res.aln_types["unmapped"] += 1
            res.read_lengths_by_mapstat["unmapped"][l_read] += 1
            for p in parts:
                res.oligo_by_mapstat["unmapped"][p] += 1
                res.oligo_by_tag["NA"][p] += 1
            continue

        minus_strand = r.is_reverse
        if r.has_tag("gf"):
            tag = gf_prio(r.get_tag("gf"))
        else:
            tag = "INTERGENIC"

        res.tag_types[tag] += 1

        is_uniq = r.mapping_quality == 255
        if not is_uniq:
            mapstat = "multimapper"
        else:
            mapstat = "unique"
            res.uniq_tag_types[tag] += 1

        for p in parts:
            res.oligo_by_mapstat[mapstat][p] += 1
            res.oligo_by_tag[tag][p] += 1

        res.aln_types[mapstat] += 1
        res.read_lengths_by_mapstat[mapstat][l_read] += 1
        if mapstat in ["multimapper", "unique"]:
            res.aln_types["mapped"] += 1
            res.read_lengths_by_mapstat["mapped"][l_read] += 1
            for p in parts:
                res.oligo_by_mapstat["mapped"][p] += 1

        cigar = list(r.cigartuples)
        cigtype = "".join([op_dict[c[0]] for c in cigar])
        cigtype = coarsegrain_CIGAR(cigtype)  # ignore splicing
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

        if is_uniq:
            res.match_len_by_utag[tag][n_match] += 1
            res.match_len_by_utag["all"][n_match] += 1

    res.N_total = res.aln_types["N_reads"]
    return res


def make_plots(res, required_aln_types=["unmapped", "mapped", "multimapper", "unique"]):
    import matplotlib.pyplot as plt

    # preparing count dictionaries for the donut plots (pie charts)
    aln_types_d, _ = rep.count_dict_collapse_misc(
        res.aln_types, total=res.N_total, add_up="N_reads"
    )
    # ensure that these keys are always present!
    for t in required_aln_types:
        aln_types_d[t] = res.aln_types.get(t, 0)

    # cig_types_d, _ = rep.count_dict_collapse_misc(res.cigar_types, total=res.N_total)
    tag_types_d, _ = rep.count_dict_collapse_misc(res.tag_types, total=res.N_total)
    utag_types_d, _ = rep.count_dict_collapse_misc(
        res.uniq_tag_types, total=res.aln_types["unique"]
    )

    # preparing the plot
    fig, axes = plt.subplots(
        3, 3, figsize=(14, 12)
    )  # gridspec_kw={'height_ratios': [1, 1, 1]})
    axes = np.array(axes).flatten()

    # # First row: donut plots
    # ciglabels, cigcolors = rep.donut_plot(
    #     axes[0], cig_types_d, title="CIGAR types", cmap="Accent"
    # )
    alnlabels, alncolors = rep.donut_plot(
        axes[0],
        aln_types_d,
        title="Alignment status",
        labels=["unique", "multimapper", "unmapped"],
        cmap="Set2",
    )

    taglabels, tagcolors = rep.donut_plot(
        axes[1], tag_types_d, title="annotation", cmap="Dark2"
    )

    utaglabels, utagcolors = rep.donut_plot(
        axes[2], utag_types_d, title="annotation (uniq only)", cmap="Dark2"
    )
    mean_lens_d = {}
    for mapstat, lens in res.read_lengths_by_mapstat.items():
        # print(mapstat, lens)
        l, f = np.array(list(lens.items())).T
        mean_lens_d[mapstat] = (l * f).sum() / f.sum()

    print(mean_lens_d)

    # Second row: histograms
    rep.len_plot(
        axes[3],
        res.read_lengths_by_mapstat,
        title="alignment status",
        labels=alnlabels + ["all"],
        colors=alncolors + ["k"],
        xlabel="read length [nt]",
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
        res.match_len_by_utag,
        labels=utaglabels + ["all"],
        colors=utagcolors + ["k"],
        title="gf (annotation) types (uniq only)",
    )

    # Third row: cumulative histograms of no. matching bases
    rep.len_plot(
        axes[6],
        res.read_lengths_by_mapstat,
        title="alignment status",
        labels=alnlabels + ["all"],
        colors=alncolors + ["k"],
        xlabel="read length [nt]",
        cumulative=True,
        ylabel="cumulative fraction",
        legend=False,
    )

    # rep.len_plot(
    #     axes[8],
    #     res.match_len_by_cigtype,
    #     labels=ciglabels + ["all"],
    #     colors=cigcolors + ["k"],
    #     title="CIGAR types",
    #     cumulative=True,
    #     ylabel="cumulative fraction",
    #     legend=False,
    # )

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
        res.match_len_by_utag,
        labels=utaglabels + ["all"],
        colors=utagcolors + ["k"],
        title="gf (annotation) types (uniq only)",
        cumulative=True,
        ylabel="cumulative fraction",
        legend=False,
    )

    # rep.len_plot(
    #     axes[10],
    #     res.sc_lengths,
    #     title="read position",
    #     colors=["r", "b"],
    #     xlabel="soft-clipped bases",
    #     cumulative=True,
    #     ylabel="cumulative fraction",
    #     legend=False,
    # )

    # rep.len_plot(
    #     axes[4],
    #     res.match_len_by_cigtype,
    #     labels=ciglabels + ["all"],
    #     colors=cigcolors + ["k"],
    #     title="CIGAR types",
    # )

    # rep.len_plot(
    #     axes[6],
    #     res.sc_lengths,
    #     title="read position",
    #     xlabel="soft-clipped bases",
    #     colors=["r", "b"],
    # )

    # axes[3].axis("off")

    return fig, (aln_types_d, tag_types_d, utag_types_d, mean_lens_d)  # cig_types_d,


def oligo_analysis(res):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import fisher_exact

    all_oligos = sorted(res.all_oligos)
    N = res.N_total
    tests = []
    interactions = []
    int_pvals = []
    int_odds = []
    table_data = []
    for oli in all_oligos:
        nuniq_w_oli = res.oligo_by_mapstat["unique"][oli]
        nuniq_wo_oli = res.aln_types["unique"] - nuniq_w_oli

        nmulti_w_oli = res.oligo_by_mapstat["multimapper"][oli]
        nmulti_wo_oli = res.aln_types["multimapper"] - nmulti_w_oli

        nmap_w_oli = res.oligo_by_mapstat["mapped"][oli]
        nmap_wo_oli = res.aln_types["mapped"] - nmap_w_oli
        nunmap_w_oli = res.oligo_by_mapstat["unmapped"][oli]
        nunmap_wo_oli = res.aln_types["unmapped"] - nunmap_w_oli

        table_map_vs_unmapped = (
            np.array([[nmap_w_oli, nmap_wo_oli], [nunmap_w_oli, nunmap_wo_oli]]) + 1
        )  # +1 = pseudo count

        table_uniq_vs_multi = (
            np.array([[nuniq_w_oli, nuniq_wo_oli], [nmulti_w_oli, nmulti_wo_oli]]) + 1
        )  # +1 = pseudo count

        for test, table in zip(
            ["mapped_vs_unmapped", "unique_vs_multimapper"],
            [table_map_vs_unmapped, table_uniq_vs_multi],
        ):
            # print(test)
            # print(table)
            # plus 1 is pseudo count
            odds, pval = fisher_exact(table + 1)
            # pval = max(pval, 1e-100)
            if pval < 0.05:
                # print(f"{test} {oli} -> odds={odds:.3e} pval={pval:.3e}")
                tests.append(test)
                interactions.append(oli)
                int_pvals.append(pval)
                int_odds.append(odds)
                table_data.append(table.ravel())

    # fig, ax = plt.subplots(figsize=(10, 7))
    # ax.spines["right"].set_visible(False)
    # ax.spines["top"].set_visible(False)
    # x = np.array(int_pvals)
    # y = np.log2(int_odds)
    # print(y)
    # ax.semilogx(x, y, ".")
    # mx, Mx = x.min(), x.max()
    # my, My = y.min(), y.max()
    # My = max(abs(my), abs(My)) * 1.1
    # wx = Mx / mx
    # wy = 2 * My

    # for x, y, label in zip(int_pvals, np.log2(int_odds), interactions):
    #     dy = (np.random.random() - 0.5) * 0.1 * wy + 0.01 * wy
    #     ax.text(x, y + dy, str(label), horizontalalignment="center")

    # ax.axhline(0, color="k", linewidth=0.5, ls="dashed")
    # ax.set_xlabel("P-value (Fisher exact)")
    # ax.set_ylabel("log2 odds ratio (oligo present vs absent)")
    # ax.set_xlim(mx * 10 ** (-0.2 * np.log10(wx)), 1)
    # ax.set_ylim(-My, My)
    # fig.tight_layout(rect=[0.2, 0.2, 0.8, 0.8])

    if len(table_data):
        a, b, c, d = np.array(table_data).T
    else:
        a, b, c, d = [], [], [], []
    df = pd.DataFrame(
        dict(
            test=tests,
            name=interactions,
            p=int_pvals,
            odds=int_odds,
            a=a,
            b=b,
            c=c,
            d=d,
        )
    )
    return df


def cmdline():
    args = parse_args()
    path = args.fname.split("/")

    if len(path) >= 4 and path[-4].startswith("sts"):
        sample_name = path[-4]
    else:
        sample_name = os.path.basename(args.fname).rsplit(".", maxsplit=1)[0]

    logging.info(f"Starting analysis of sample_name='{sample_name}'")
    res = scan_bam(
        args.fname,
        skim=args.skim,
        intact_signature=args.intact_signature,
        parse_oligos=args.parse_oligos,
    )
    # print(res.aln_types)
    if args.parse_oligos:
        df = oligo_analysis(res)
        df["sample"] = sample_name
        # fig.suptitle(sample_name)
        # fig.savefig(f"{args.out_pdf}/{sample_name}_oligo_analysis.pdf")
        # fig.savefig(f"{args.out_png}/{sample_name}_oligo_analysis.png", dpi=300)
        df.to_csv(
            f"{args.out_csv}/{sample_name}.oligo_analysis.csv", sep="\t", index=False
        )
    import matplotlib

    matplotlib.use("Agg")

    fig, (aln_counts, tag_counts, utag_counts, mean_lens_d) = make_plots(res)
    if args.intact_signature:
        fig.suptitle(f"{sample_name} ({args.intact_signature})")
    else:
        fig.suptitle(sample_name)

    if args.intact_signature:
        suffix = f"_intact_sig"
    else:
        suffix = ""

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(f"{args.out_pdf}/{sample_name}{suffix}.pdf")
    fig.savefig(f"{args.out_png}/{sample_name}{suffix}.png", dpi=300)

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

    from spacemake.longread.report import count_dict_to_df

    df = pd.concat(
        [
            # count_dict_to_df(cig_counts, kind="CIGAR_types", n_total=res.N_total),
            count_dict_to_df(aln_counts, kind="align_types", n_total=res.N_total),
            count_dict_to_df(tag_counts, kind="tag_types", n_total=res.N_total),
            count_dict_to_df(
                utag_counts, kind="uniq_tag_types", n_total=aln_counts["unique"]
            ),
        ],
        ignore_index=True,
    )
    for mapstat, meanlen in mean_lens_d.items():
        df = df.append(
            dict(
                name=mapstat,
                count=meanlen,
                fraction=1.0,
                kind="mean_readlen",
            ),
            ignore_index=True,
        )

    df.to_csv(f"{args.out_csv}/{sample_name}.csv", sep="\t")


if __name__ == "__main__":
    cmdline()
