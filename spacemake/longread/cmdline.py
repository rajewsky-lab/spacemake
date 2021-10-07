__author__ = ["Marvin Jens"]
__license__ = "MIT"
__email__ = ["marvin.jens@mdc-berlin.de"]

import os
import argparse
import logging
import pandas as pd
import numpy as np
import spacemake.longread.util as util
import spacemake.longread.report as report
import spacemake.longread.cache as cache
import spacemake.longread.annotation as ann
from collections import defaultdict


def detect_sample(args):
    if args.sample is None:
        sample_name = os.path.splitext(os.path.basename(args.fname))[0]
        logging.info(f"auto-detected sample_name={sample_name}")
    else:
        sample_name = args.sample
    return sample_name.replace(".stats", "")


def aln_main(args):
    sample_name = detect_sample(args)
    blocks = util.load_oligos(args.blocks)

    cache.fill_caches(
        args.fname, sample_name, blocks, path=args.cache, n_proc=args.parallel
    )
    df = cache.annotate(args.fname, sample_name, blocks, path=args.cache)
    df.to_csv(
        os.path.join(
            util.ensure_path(args.annotation_out), f"{sample_name}.annotation.tsv"
        ),
        sep="\t",
        index=False,
    )


def ann_main(args):
    logger = logging.getLogger("longread.cmd.annotation")
    sample_name = detect_sample(args)
    blocks = util.load_oligos(args.blocks)

    sig_intact = args.intact_bead.split(",")
    sig_core, sig_core_order = util.process_intact_signature(args.intact_bead)

    logger.info(f"analyzing {sample_name} for intact signature {sig_intact}")
    annotation = ann.AnnotatedSequences(
        args.fname,
        os.path.join(args.annotation_out, f"{sample_name}.annotation.tsv"),
        sample_name,
        blocks,
        min_score=args.min_score,
        orient_by=sig_core[0],
    )
    n_total = len(annotation.raw_sequences)
    logging.info(f"total number of reads in {args.fname} ({sample_name}) is {n_total}")

    sig_counts = annotation.count_signatures()
    util.count_dict_out(sig_counts, "common signatures", total=n_total)

    df_sig = util.count_dict_to_df(sig_counts, "signatures", n_total=n_total)

    partial_counts, prefixes, suffixes, pT_counts = annotation.completeness(
        sig_core, polyT=args.polyT
    )

    partial_counts_simple, _ = util.count_dict_collapse_misc(
        partial_counts, sig_intact=sig_intact, total=n_total, misc_thresh=0.00001
    )
    util.count_dict_out(partial_counts_simple, "completeness", total=n_total)

    df_comp = util.count_dict_to_df(
        partial_counts_simple, kind="bead_complete", n_total=n_total
    ).sort_values("name")

    df_pT = util.count_dict_to_df(pT_counts, kind="polyT_after", n_total=n_total)

    fname = os.path.join(util.ensure_path(args.stats_out), f"{sample_name}.stats.tsv")
    logging.info(f"storing annotation signature counts as DataFrame '{fname}'")
    df = pd.concat([df_sig, df_comp, df_pT])
    df.to_csv(fname, sep="\t", index=False)

    # TODO: prefix/suffix counts add up to > 100%. Needs fix
    util.count_dict_out(pT_counts, "polyT after", total=n_total)
    util.count_dict_out(prefixes, "prefixes", total=n_total)
    util.count_dict_out(suffixes, "suffixes", total=n_total)

    # Gather statistics about the parts that make up intact oligos
    qintact, qL, qstarts, qends, qscores = annotation.query_dimensions(
        sig_core, substring=True
    )
    # print(qstarts.shape, qintact.shape, qL.shape)
    from collections import Counter

    data = []
    for part, starts, ends, scores in zip(sig_core, qstarts.T, qends.T, qscores.T):
        starts_hist = sorted(Counter(starts).items())
        ends_hist = sorted(Counter(ends).items())
        lens_hist = sorted(Counter(ends - starts).items())
        scores_hist = sorted(Counter(scores).items())

        for x, f in starts_hist:
            data.append((args.intact_bead, part, "start", x, f))

        for x, f in ends_hist:
            data.append((args.intact_bead, part, "end", x, f))

        for x, f in lens_hist:
            data.append((args.intact_bead, part, "len", x, f))

        for x, f in scores_hist:
            data.append((args.intact_bead, part, "score", x, f))

    # For reference, also gather statistics for each part of the intact signature,
    # regardless of whether it occurs in the context of an intact signature or not,
    # i.e. what do *all matches* look like
    for part in sig_intact:
        qnames, starts, ends, scores, qL = annotation.query_oligo_occurrences(part)
        starts_hist = sorted(Counter(starts).items())
        ends_hist = sorted(Counter(ends).items())
        lens_hist = sorted(Counter(ends - starts).items())
        scores_hist = sorted(Counter(scores).items())

        for x, f in starts_hist:
            data.append(("anywhere", part, "start", x, f))

        for x, f in ends_hist:
            data.append(("anywhere", part, "end", x, f))

        for x, f in lens_hist:
            data.append(("anywhere", part, "len", x, f))

        for x, f in scores_hist:
            data.append(("anywhere", part, "score", x, f))

    df_parts = pd.DataFrame(
        data, columns=["signature", "oligo", "attr", "value", "freq"]
    )
    fname = os.path.join(
        util.ensure_path(args.stats_out), f"{sample_name}.intact_parts.tsv"
    )
    logging.info(f"storing intact signature part-statistics as DataFrame '{fname}'")
    df_parts.to_csv(fname, sep="\t", index=False)

    # plot score distributions
    qnames, starts, ends, scores, qlens = annotation.query_oligo_occurrences(
        "bead_start"
    )

    # output representative examples
    eo_fname = util.ensure_path(os.path.join(args.examples_out, f"{sample_name}.txt"))
    with open(eo_fname, "wt") as eo:
        for signame, sigcount in sorted(sig_counts.items(), key=lambda x: -x[1]):
            try:
                qname, _, _ = next(
                    annotation.filter_signatures(tuple(signame.split(",")))
                )
                eo.write(f"# {signame} n={sigcount}\n{annotation.fmt(qname)}\n")
            except StopIteration:
                logger.warning(f"unable to find any reads with signature {signame}")


def rep_main(args):
    logger = logging.getLogger("longread.report")
    logger.info(f"generating reports from '{args.fname}'")
    sample_name = detect_sample(args)
    fname = os.path.join(util.ensure_path(args.stats_out), f"{sample_name}.stats.tsv")
    df = pd.read_csv(fname, sep="\t", comment="#").fillna("other")
    print(df)
    sig_counts = util.count_dict_from_df(df, "signatures")
    # print(sig_counts)
    util.count_dict_out(
        sig_counts, "signatures", misc_thresh=0.05, total=sig_counts["n_total"]
    )
    print(
        "summed up counts in sig_counts dict", np.array(list(sig_counts.values())).sum()
    )
    bead_related = args.bead_related
    if bead_related is None:
        sig_intact = tuple(args.intact_bead.split(","))
        bead_related = sig_intact[0]

    n_total = sig_counts["n_total"]
    print(f"oligo name flagging 'bead-related' signature: {bead_related}")

    ov_counts, bead_counts, found_part_counts, core_signature = util.digest_signatures(
        sig_counts,
        bead_related,
        args.intact_bead,
        # prefixes=args.prefixes,
        # suffixes=args.suffixes,
    )
    ov_simple, _ = util.count_dict_out(
        ov_counts, "overview counts", misc_thresh=0.01, total=sig_counts["n_total"]
    )
    bead_simple, _ = util.count_dict_out(
        bead_counts, "bead counts", misc_thresh=0.0, total=ov_counts["bead-related"]
    )
    print("found_part_counts")
    for k, v in found_part_counts.items():
        print(k, v)

    print("=" * 30)
    del ov_simple["n_total"]

    print(
        "summed up counts after splitting and dropping n_total",
        np.array(list(ov_simple.values())).sum(),
    )
    print(f"n_total={n_total}")

    ov_items = sorted(ov_simple.items(), key=lambda x: -x[1])
    ov_labels = [x[0] for x in ov_items]
    ov_counts = [x[1] for x in ov_items]

    bead_items = sorted(bead_simple.items(), key=lambda x: -x[1])

    bead_labels = [x[0] for x in bead_items]
    bead_counts = [x[1] for x in bead_items]

    fname = os.path.join(util.ensure_path(args.report_out), f"{sample_name}.donuts.pdf")
    syn_rates = report.plot_results(
        sig_counts,
        ov_labels,
        ov_counts,
        bead_labels,
        bead_counts,
        found_part_counts,
        all_parts=core_signature,
        fname=fname,
        suptitle=sample_name,
    )
    # print(core_signature)
    # print(syn_rates)
    fname = os.path.join(util.ensure_path(args.stats_out), f"{sample_name}.synth.tsv")
    df = pd.DataFrame(dict(segment=core_signature, rate=syn_rates))
    df["sample"] = sample_name
    print(df)
    df.to_csv(fname, sep="\t")

    fname = os.path.join(
        util.ensure_path(args.stats_out), f"{sample_name}.intact_parts.tsv"
    )
    df = pd.read_csv(fname, sep="\t", comment="#").fillna("other")
    fname = os.path.join(util.ensure_path(args.report_out), f"{sample_name}.hists.pdf")
    report.plot_histograms(
        df, fname, n_total=n_total, parts=args.intact_bead.split(",")
    )


def main_edits(args):
    sample_name = detect_sample(args)
    blocks = util.load_oligos(args.blocks)
    sig_intact = tuple(args.intact_bead.split(","))

    annotation = ann.AnnotatedSequences(
        args.fname,
        os.path.join(args.annotation_out, f"{sample_name}.annotation.tsv"),
        sample_name,
        blocks,
        min_score=args.min_score,
    )
    n_total = len(annotation.raw_sequences)
    logging.info(f"total number of reads in {args.fname} ({sample_name}) is {n_total}")

    data = []
    for part in sig_intact:
        qmatches = annotation.query_oligo_occurrences(part)
        if len(qmatches[0]) > args.n_samples:
            qmatches = ann.subsample(qmatches, n=args.n_samples)

        nmatch = len(qmatches[0])
        m, ed = ann.align_stats(annotation, blocks[part], qmatches)
        for x in np.arange(len(m)):
            # print(part, x, m)
            data.append((part, blocks[part], nmatch, x, m[x], ed[x]))

    df = pd.DataFrame(
        data, columns=["oligo", "seq", "nmatch", "pos", "fmatch", "ed_dict"]
    )
    # print(df)
    df.to_csv(
        os.path.join(args.stats_out, f"{sample_name}.oligo_edits.tsv"),
        sep="\t",
        index=False,
    )

    report.plot_edits(
        df,
        os.path.join(args.report_out, f"{sample_name}.oligo_edits.pdf"),
        parts=sig_intact,
    )


def main_extract(args):
    sample_name = detect_sample(args)
    blocks = util.load_oligos(args.blocks)

    annotation = ann.AnnotatedSequences(
        args.fname,
        os.path.join(args.annotation_out, f"{sample_name}.annotation.tsv"),
        sample_name,
        blocks,
        min_score=args.min_score,
    )
    n_total = len(annotation.raw_sequences)
    logging.info(f"total number of reads in {args.fname} ({sample_name}) is {n_total}")

    anchor_scores = defaultdict(float)
    barcodes = defaultdict(lambda: "NA")
    umis = defaultdict(lambda: "NA")

    cb_start, cb_end = args.CB.split(",")
    cb_start, cb_end = int(cb_start), int(cb_end)

    umi_start, umi_end = args.UMI.split(",")
    umi_start, umi_end = int(umi_start), int(umi_end)

    if args.barcode_after:
        hits = annotation.query_oligo_occurrences(args.barcode_after)
        for qname, start, end, score, L in zip(*hits):
            a_score = anchor_scores[qname]
            if score > a_score:
                anchor_scores[qname] = score
                seq = annotation.raw_sequences[qname]
                barcodes[qname] = seq[end + cb_start : end + cb_end]
                umis[qname] = seq[end + umi_start : end + umi_end]

    if args.top_barcodes:
        known = set([bc.strip()[::-1] for bc in open(args.top_barcodes).readlines()])
        rev = set([bc[::-1] for bc in known])
        detect = set(barcodes.values())
        nd = len(detect)
        ov = len(detect & known)
        ctrl = len(detect & rev)
        logging.info(
            f"loaded {len(known)} barcodes from '{args.top_barcodes}'. "
            f"{ov} / {nd} detected barcodes overlap ({ov/nd * 100:.3f}%). "
            f"Reverse BC control {ctrl/nd * 100:.3f}%"
        )
    else:
        known = set()

    n = 0
    for qname, sig in annotation.signatures.items():
        sig_str = ",".join(sig)
        if args.sig_include and not (args.sig_include in sig_str):
            continue

        if args.sig_exclude and (args.sig_exclude in sig_str):
            continue

        if args.cDNA_after in sig:
            (
                (cDNA_start, cDNA_end),
                (start_oli, end_oli),
                cDNA,
            ) = annotation.extract_cDNA(
                qname, after_oligo=args.cDNA_after, distal=args.distal
            )
            bc = barcodes[qname]
            if bc not in known:
                bc = bc.lower()
            n += 1
            seq = annotation.raw_sequences[qname]
            header = (
                f">{n}__CB:{bc}__UMI:{umis[qname]}__"
                + f"sig:{','.join(sig)}__cDNA:{cDNA_start}-{cDNA_end}__oli:{start_oli}-{end_oli}__"
                + f"L_read={len(seq)}__L_cDNA={cDNA_end - cDNA_start}"
            )

            # print(header[: (254 - 14)] + " 1:N:0:TCCTGAGC" + f"\n{cDNA}")
            print(header[: (254 - 14)] + " 1:N:0:TCCTGAGC" + f"\n{cDNA}")


def load_tables_overview(sample_fnames, pattern):
    dfs = []
    for fname in sample_fnames:
        sample_name = os.path.splitext(os.path.basename(fname))[0]
        df = pd.read_csv(
            pattern.format(sample_name=sample_name),
            sep="\t",
        )
        df["sample"] = sample_name
        dfs.append(df)

    df = pd.concat(dfs)
    samples = sorted(df["sample"].drop_duplicates())

    return df, samples


class SampleDB:
    def __init__(self, intact={}, color={}, label={}, prio={}):
        self.intact = intact
        self.color = color
        self.label = label
        self.prio = prio

    @classmethod
    def from_YAML(cls, fname="samples.yaml"):
        import yaml

        groups = yaml.load(open(fname), Loader=yaml.SafeLoader)
        print(groups)
        grp = groups["groups"]
        default = grp[groups["default"]]
        print("default THING", default)

        intact_lkup = defaultdict(lambda: default["intact_bead"])
        color_lkup = defaultdict(lambda: default["color"])
        label_lkup = defaultdict(lambda: default["label"])
        prio_lkup = defaultdict(lambda: default["prio"])
        for name, d in groups["groups"].items():
            print(f"name={name} d={d}")
            for sample in d["samples"]:
                intact_lkup[sample] = d["intact_bead"]
                color_lkup[sample] = d["color"]
                label_lkup[sample] = d["label"]
                prio_lkup[sample] = d["prio"]

        print(color_lkup)
        print("DEFAULT COLOR", color_lkup["TEST"])
        return cls(intact_lkup, color_lkup, label_lkup, prio_lkup)

    def sort_samples(self, samples):
        return sorted(samples, key=lambda x: (self.prio.get(x, np.inf), x))


def synthesis_rates_overview(args):
    df, samples = load_tables_overview(
        args.overview_samples, os.path.join(args.stats_out, "{sample_name}.synth.tsv")
    )
    df["segment"] = df["segment"].apply(lambda x: x.replace("_2s", ""))

    db = SampleDB.from_YAML("samples.yaml")
    samples = db.sort_samples(samples)

    row_labels = ["bead_start", "OP1", "OP2", "OP3", "polyT"]
    dfs = [df.query(f"segment == '{part}'") for part in row_labels]

    fig, axes = report.multi_row_barplots(
        dfs, row_labels, samples, "rate", color_dict=db.color
    )
    fig.suptitle("apparent synthesis rates")

    fig.tight_layout()
    fig.savefig("overview_syn_rates.pdf")
    return df


def main_overview(args):
    return synthesis_rates_overview(args)


def prepare_parser():
    # import pb
    # import overview as ov

    parser = argparse.ArgumentParser(prog="longread")
    parser.add_argument(
        "fname", default=None, help="file with pacbio reads (FASTQ or BAM format)"
    )
    parser.add_argument(
        "--sample",
        help="sample name (default=autodetect from fname)",
        default=None,
    )
    # parser.add_argument(
    #     "--aln-cache", help="alignment cache folder", default="./cache/"
    # )
    parser.add_argument(
        "--blocks",
        default="",
        help="FASTA file with known oligo sequences (default=built-in)",
    )
    parser.add_argument(
        "--cache",
        default="./cache/",
        help="path to alignment caches (default=./cache/)",
    )
    parser.add_argument(
        "--annotation-out",
        default="./annotation/",
        help="path to store annotation data in (default=./annotation/)",
    )
    parser.add_argument(
        "--parallel",
        default=16,
        type=int,
        help="number of parallel processes (default=16)",
    )
    parser.add_argument(
        "--min-score",
        default=0.6,
        type=float,
        help="minimal match alignment score to consider a match for annotation, relative to its size (default=0.6)",
    )
    parser.add_argument(
        "--debug", default=False, action="store_true", help="activate debug output"
    )
    parser.add_argument(
        "--intact-bead",
        default="bead_start,OP1,pT",
        help="sequence of oligos that correspond to intact bead sequence in correct order",
    )
    parser.add_argument(
        "--bead-related",
        default=None,
        help="name of oligo-block that is specific to the capture oligo. default=first element of --intact-bead",
    )

    parser.add_argument(
        "--stats-out",
        help="path to store statistics (pandas dataframes)",
        default="./stats/",
    )
    parser.add_argument(
        "--report-out", help="path to render graphical reports in", default="./reports/"
    )

    ## sub-parser setup ##
    subparsers = parser.add_subparsers(help="sub-command help")
    aln_parser = subparsers.add_parser(
        "align", help="align PacBio reads against oligos"
    )

    ann_parser = subparsers.add_parser(
        "annotate", help="create annotation from detected oligo matches"
    )
    ann_parser.add_argument(
        "--examples-out",
        help="path to store annotated example read sequences",
        default="./examples/",
    )
    ann_parser.add_argument(
        "--polyT",
        help="name of oligo-block that matches polyT (default='polyT')",
        default="polyT",
    )

    rep_parser = subparsers.add_parser("report", help="create PDF/PNG reports")

    ed_parser = subparsers.add_parser(
        "edits", help="gather mismatch and indel stats for oligo matches"
    )
    ed_parser.add_argument(
        "--n-samples",
        type=int,
        default=1000,
        help="number of sample alignments to gather edit statistics from (default=1000)",
    )

    xt_parser = subparsers.add_parser(
        "extract", help="extract sequences from the long read"
    )
    xt_parser.add_argument(
        "--barcode-after",
        type=str,
        default="bead_start",
        help="name of anchor match after which the barcodes follow (default='bead_start')",
    )
    xt_parser.add_argument(
        "--CB",
        type=str,
        default="8,20",
        help="bases downstream of anchor match at which the cell barcode starts and ends. Default='8,20'",
    )
    xt_parser.add_argument(
        "--UMI",
        type=str,
        default="0,8",
        help="bases downstream of anchor match at which the UMI starts and ends. Default='0,8'",
    )
    xt_parser.add_argument(
        "--cDNA-after",
        type=str,
        default="bead_start",
        help="excise sequence between this oligo and the last oligo match (if it exists)",
    )
    xt_parser.add_argument(
        "--sig-include",
        type=str,
        default="",
        help="extract cDNA only from long reads INCLUDING this substring in the signature",
    )
    xt_parser.add_argument(
        "--sig-exclude",
        type=str,
        default="",
        help="extract cDNA only from long reads EXCLUDING this substring in the signature",
    )
    xt_parser.add_argument(
        "--distal",
        default=150,
        type=int,
        help="number of nt from end of sequence that are considered for oligo matches when extracting cDNA",
    )

    xt_parser.add_argument(
        "--top-barcodes",
        type=str,
        default="",
        help="path to text file with known barcodes (e.g. from Illumina)",
    )

    ov_parser = subparsers.add_parser(
        "overview", help="make overview plots across samples"
    )
    ov_parser.add_argument(
        "overview_samples", nargs="+", help="sample names to aggregate"
    )

    parser.set_defaults(func=lambda args: parser.print_help())
    aln_parser.set_defaults(func=aln_main)
    ann_parser.set_defaults(func=ann_main)
    rep_parser.set_defaults(func=rep_main)
    ed_parser.set_defaults(func=main_edits)
    xt_parser.set_defaults(func=main_extract)
    ov_parser.set_defaults(func=main_overview)

    return parser


def cmdline():
    import logging

    logging.basicConfig(level=logging.INFO)

    parser = prepare_parser()
    args = parser.parse_args()
    return args.func(args)
