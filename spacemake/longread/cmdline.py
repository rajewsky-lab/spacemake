__author__ = ["Marvin Jens"]
__email__ = ["marvin.jens@mdc-berlin.de"]

import os
import argparse
import logging
import pandas as pd
import numpy as np
from spacemake.util import ensure_path
import spacemake.longread.report as report
import spacemake.longread.cache as cache
import spacemake.longread.annotation as ann
from spacemake.longread.signature import (
    get_signature_db,
    process_intact_signature,
    digest_signatures,
)
from collections import defaultdict


def detect_sample(args):
    if args.sample is None:
        if not hasattr(args, "fname"):
            raise ValueError("please provide a sample-identifier via --sample")

        sample_name, _ = os.path.splitext(os.path.basename(args.fname))
        logging.info(f"auto-detected sample_name={sample_name}")
    else:
        sample_name = args.sample
    return sample_name.replace(".stats", "")


def store_results(df, path, fname, logger, **kw):
    fpath = os.path.join(ensure_path(path), fname)
    df.to_csv(fpath, sep="\t", **kw)
    logger.info(f"storing {len(df)} rows to '{fpath}'")
    return fpath


def load_results(path, fname, logger, **kw):
    fpath = os.path.join(path, fname)
    df = pd.read_csv(fpath, sep="\t", comment="#", **kw).fillna("other")
    logger.info(f"loaded {len(df)} rows from '{fpath}'")
    return df


def setup_namespace(args, need_sample_name=True):
    """
    perform initialization code common to more than one sub-command.
    Augments the args object returned from argparse module with additional
    variables that are used by multiple sub-commands.

    :param args: cmdline parser args
    :type args: Namespace
    :return: augmented namespace
    :rtype: SimpleNamespace
    """
    d = vars(args)
    if need_sample_name:
        sample = detect_sample(args)
        d["sample_name"] = sample

    db = get_signature_db(args.config)

    d["signature_db"] = db
    d["CB"] = db.CB[args.signature]
    d["UMI"] = db.UMI[args.signature]
    d["read1_primer"] = db.read1_primer[args.signature]
    d["read2_primer"] = db.read2_primer[args.signature]
    d["cDNA_after"] = db.cDNA_after[args.signature]
    intact_bead = db.intact[args.signature]
    d["intact_bead"] = intact_bead
    d["relevant"] = (
        intact_bead.split(",")
        + db.other[args.signature].split(",")
        + db.prefixes[args.signature].split(",")
        + db.suffixes[args.signature].split(",")
    )
    d["sig_intact"] = intact_bead.split(",")
    d["sig_core"], d["sig_core_order"] = process_intact_signature(intact_bead)
    d["bead_related"] = d["sig_intact"][0]
    d["logger"] = logging.getLogger(f"spacemake.longread.cmd.{args.subcmd}")
    d["blocks"] = db.blocks

    from types import SimpleNamespace

    return SimpleNamespace(**d)


def initialize(args):
    args = setup_namespace(args)

    annotation = ann.AnnotatedSequences(
        args.fname,
        os.path.join(args.annotation_out, f"{args.sample_name}.annotation.tsv"),
        args.sample_name,
        args.blocks,
        min_score=args.min_score,
        relevant=args.relevant,
        orient_by=args.bead_related,
    )
    n_total = len(annotation.raw_sequences)

    args.logger.info(
        f"total number of reads in {args.fname} ({args.sample_name}) is {n_total}"
    )
    return args, annotation


## Subcommand 'align'
def cmd_align(args):
    args = setup_namespace(args)

    cache.fill_caches(
        args.fname,
        args.sample_name,
        args.blocks,
        relevant=args.relevant,
        path=ensure_path(args.cache),
        n_proc=args.parallel,
    )

    df = cache.annotate(
        args.fname,
        args.sample_name,
        args.blocks,
        path=ensure_path(args.cache),
        relevant=args.relevant,
    )
    df.to_csv(
        os.path.join(
            ensure_path(args.annotation_out), f"{args.sample_name}.annotation.tsv"
        ),
        sep="\t",
        index=False,
    )


def prepare_align_parser(subparsers):
    parser = subparsers.add_parser("align", help="align PacBio reads against oligos")
    parser.add_argument(
        "fname",
        default=None,
        help="file with pacbio reads (FASTQ or BAM format)",
    )
    parser.set_defaults(func=cmd_align)
    return parser


## subcommand 'annotate'
def cmd_annotate(args):
    args, annotation = initialize(args)
    n_total = len(annotation.raw_sequences)

    sig_counts, n_concat, n_reprimed = annotation.count_signatures()
    report.count_dict_out(sig_counts, "common signatures", total=n_total)
    print(
        f"n_concat={n_concat} ({100.0 * n_concat/n_total:.2f}%) "
        f"n_reprimed={n_reprimed} ({100.0 * n_reprimed/n_total:.2f}%)"
    )
    df_sig = report.count_dict_to_df(sig_counts, "signatures", n_total=n_total)

    concat_counts, n_occurrences = annotation.count_concatenations()
    report.count_dict_out(concat_counts, "concatenations", total=n_occurrences)
    df_concat = report.count_dict_to_df(
        concat_counts, "concatenations", n_total=n_occurrences
    )

    reprimed_counts = annotation.count_repriming()
    report.count_dict_out(reprimed_counts, "repriming", total=n_total)
    df_reprimed = report.count_dict_to_df(reprimed_counts, "repriming", n_total=n_total)

    partial_counts, prefixes, suffixes, pT_counts = annotation.completeness(
        args.sig_core, polyT=args.polyT
    )

    partial_counts_simple, _ = report.count_dict_collapse_misc(
        partial_counts, sig_intact=args.sig_intact, total=n_total, misc_thresh=0.00001
    )
    report.count_dict_out(partial_counts_simple, "completeness", total=n_total)

    df_comp = report.count_dict_to_df(
        partial_counts_simple, kind="bead_complete", n_total=n_total
    ).sort_values("name")

    df_pT = report.count_dict_to_df(pT_counts, kind="polyT_after", n_total=n_total)

    df = pd.concat([df_sig, df_concat, df_reprimed, df_comp, df_pT])
    store_results(df, args.stats_out, f"{args.sample_name}.stats.tsv", args.logger)

    # TODO: prefix/suffix counts add up to > 100%. Needs fix
    report.count_dict_out(pT_counts, "polyT after", total=n_total)
    report.count_dict_out(prefixes, "prefixes", total=n_total)
    report.count_dict_out(suffixes, "suffixes", total=n_total)

    # Gather statistics about the parts that make up intact oligos
    qintact, qL, qstarts, qends, qscores = annotation.query_dimensions(
        args.sig_core, substring=True
    )
    # print(qstarts.shape, qintact.shape, qL.shape)
    from collections import Counter

    data = []
    for part, starts, ends, scores in zip(args.sig_core, qstarts.T, qends.T, qscores.T):
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
    for part in args.sig_intact:
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
    store_results(
        df_parts, args.stats_out, f"{args.sample_name}.intact_parts.tsv", args.logger
    )

    # output representative examples
    eo_fname = os.path.join(ensure_path(args.examples_out), f"{args.sample_name}.txt")

    with open(eo_fname, "wt") as eo:
        for signame, sigcount in sorted(sig_counts.items(), key=lambda x: -x[1]):
            i = None
            for i, (qname, _, _) in enumerate(
                annotation.filter_signatures(tuple(signame.split(",")))
            ):
                eo.write(
                    f"# example={i} {signame} n={sigcount}\n{annotation.fmt(qname)}\n"
                )
                if i >= 10:
                    break

            if i is None:
                args.logger.warning(
                    f"unable to find any reads with signature {signame}"
                )


def prepare_annotate_parser(subparsers):
    parser = subparsers.add_parser(
        "annotate", help="create annotation from detected oligo matches"
    )
    parser.add_argument(
        "fname",
        default=None,
        help="file with pacbio reads (FASTQ or BAM format)",
    )
    parser.add_argument(
        "--polyT",
        help="name of oligo-block that matches polyT (default='polyT')",
        default="polyT",
    )
    parser.set_defaults(func=cmd_annotate)
    return parser


## subcommand 'report'
def get_synth_rates(found_part_counts, all_parts, n_total):
    rates = []

    n0 = n_total
    x = range(1, len(all_parts) + 1)
    for i in x:
        key = tuple(all_parts[:i])
        # print(found_part_counts)
        # print(key, found_part_counts[key], n0, rates)
        rates.append(found_part_counts[key] / float(max(1, n0)))
        # print(i, rates)
        n0 = found_part_counts[key]

    return np.array(rates)


def cmd_report(args):
    args = setup_namespace(args)
    args.logger.info(f"generating report plots for '{args.sample_name}'")
    args.logger.debug(f"'bead-related' if we detect: '{args.bead_related}'")

    df = load_results(args.stats_out, f"{args.sample_name}.stats.tsv", args.logger)
    sig_counts = report.count_dict_from_df(df, "signatures")
    report.count_dict_out(
        sig_counts, "signatures", misc_thresh=0.05, total=sig_counts["n_total"]
    )
    n_total = sig_counts["n_total"]
    args.logger.info(f"n_total={n_total}")

    # disentangle bead-related from other signatures
    ov_counts, bead_counts, found_part_counts, core_signature = digest_signatures(
        sig_counts,
        args.bead_related,
        args.intact_bead,
        # prefixes=args.prefixes,
        # suffixes=args.suffixes,
    )
    # group low abundance signatures into 'misc' for overview donut plot
    ov_simple, _ = report.count_dict_out(
        ov_counts, "overview counts", misc_thresh=0.01, total=sig_counts["n_total"]
    )
    # group low abundance signatures into 'misc' for bead completeness donut plot
    bead_simple, _ = report.count_dict_out(
        bead_counts, "bead counts", misc_thresh=0.0, total=ov_counts["bead-related"]
    )
    # store the donut plot values in ...report.tsv for cross-sample overviews and such
    df_ov = report.count_dict_to_df(ov_simple, kind="overview")
    df_bead = report.count_dict_to_df(bead_simple, kind="bead_related")
    df_rep = pd.concat([df_ov, df_bead])
    df_rep["sample"] = args.sample_name
    df_rep["signature"] = args.signature
    print(df_rep)
    store_results(df_rep, args.stats_out, f"{args.sample_name}.report.tsv", args.logger)

    # compute and store apparent synthesis rates
    syn_rates = get_synth_rates(found_part_counts, core_signature, n_total)
    df = pd.DataFrame(dict(segment=core_signature, rate=syn_rates))
    df["sample"] = args.sample_name
    store_results(df, args.stats_out, f"{args.sample_name}.synth.tsv", args.logger)

    # prepare for plotting
    del ov_simple["n_total"]
    # assert np.array(list(ov_simple.values())).sum() == n_total
    print(
        f"sum of ov_simple {np.array(list(ov_simple.values())).sum()} n_total={n_total}"
    )

    ov_items = sorted(ov_simple.items(), key=lambda x: -x[1])
    ov_labels = [x[0] for x in ov_items]
    ov_counts = [x[1] for x in ov_items]

    bead_items = sorted(bead_simple.items(), key=lambda x: -x[1])
    bead_labels = [x[0] for x in bead_items]
    bead_counts = [x[1] for x in bead_items]

    # render the donut plots #
    report.plot_results(
        sig_counts,
        ov_labels,
        ov_counts,
        bead_labels,
        bead_counts,
        syn_rates,
        all_parts=core_signature,
        fname=os.path.join(
            ensure_path(args.report_out), f"{args.sample_name}.donuts.pdf"
        ),
        suptitle=args.sample_name,
    )

    # render the pre-computed position/score/length histograms as a plot
    report.plot_histograms(
        df=load_results(
            args.stats_out, f"{args.sample_name}.intact_parts.tsv", args.logger
        ),
        fname=os.path.join(
            ensure_path(args.report_out), f"{args.sample_name}.hists.pdf"
        ),
        n_total=n_total,
        parts=args.intact_bead.split(","),
    )


def prepare_report_parser(subparsers):
    parser = subparsers.add_parser("report", help="create PDF/PNG reports")
    parser.set_defaults(func=cmd_report)
    return parser


## subcommand 'edits'
def cmd_edits(args):
    args, annotation = initialize(args)

    data = []
    for part in args.sig_intact:
        qmatches = annotation.query_oligo_occurrences(part)
        if len(qmatches[0]) > args.n_samples:
            qmatches = ann.subsample(qmatches, n=args.n_samples)

        nmatch = len(qmatches[0])
        m, ed = ann.align_stats(
            annotation, args.blocks[part], qmatches, min_score=args.min_score
        )
        for x in np.arange(len(m)):
            # print(part, x, m)
            data.append((part, args.blocks[part], nmatch, x, m[x], ed[x]))

    df = pd.DataFrame(
        data, columns=["oligo", "seq", "nmatch", "pos", "fmatch", "ed_dict"]
    )
    store_results(
        df, args.stats_out, f"{args.sample_name}.oligo_edits.tsv", args.logger
    )

    report.plot_edits(
        df,
        os.path.join(args.report_out, f"{args.sample_name}.oligo_edits.pdf"),
        parts=args.sig_intact,
    )


def prepare_edits_parser(subparsers):
    parser = subparsers.add_parser(
        "edits", help="gather mismatch and indel stats for oligo matches"
    )
    parser.add_argument(
        "--n-samples",
        type=int,
        default=1000,
        help="number of sample alignments to gather edit statistics from (default=1000)",
    )
    parser.add_argument(
        "fname",
        default=None,
        help="file with pacbio reads (FASTQ or BAM format)",
    )
    parser.add_argument(
        "--min-score",
        default=0.6,
        type=float,
        help="minimal match alignment score to consider a match for edit extraction, relative to its size (default=0.6)",
    )
    parser.set_defaults(func=cmd_edits)
    return parser


## subcommand 'extract'
def cmd_extract(args):
    args, annotation = initialize(args)

    # anchor_scores = defaultdict(float)
    CB_detect = set()
    # umis = defaultdict(lambda: "NA")

    args.logger.info(
        "expecting read1 primed by {args.read1_primer} and CB={args.CB} UMI={args.UMI}"
    )
    args.logger.info(
        f"expecting cDNA after {args.cDNA_after} and primed by {args.read2_primer}"
    )

    if args.top_barcodes:
        known = set([bc.strip()[::-1] for bc in open(args.top_barcodes).readlines()])
        rev = set([bc[::-1] for bc in known])
        logger.info(f"loaded {len(known)} barcodes from '{args.top_barcodes}'. ")
    else:
        known = set()

    n = 0
    for qname, sig in annotation.signatures.items():
        sig_str = ",".join(sig)
        if args.sig_include and not (args.sig_include in sig_str):
            continue

        if args.sig_exclude and (args.sig_exclude in sig_str):
            continue

        r1 = ""
        if args.read1_primer in sig:
            ((r1_start, r1_end), (start_oli, end_oli), r1,) = annotation.extract_cDNA(
                qname, after_oligo=args.read1_primer, distal=args.distal
            )

        r2 = ""
        if args.cDNA_after in sig:
            (
                (cDNA_start, cDNA_end),
                (start_oli, end_oli),
                r2,
            ) = annotation.extract_cDNA(
                qname, after_oligo=args.cDNA_after, distal=args.distal
            )
            if not r2.strip() or (cDNA_end - cDNA_start <= 0):
                r2 = ""

        CB = eval(args.CB, dict(r1=r1)) if r1 else "NA"
        UMI = eval(args.UMI, dict(r1=r1, r2=r2)) if r1 and r2 else "NA"

        CB_detect.add(CB)
        if CB not in known:
            CB = CB.lower()

        n += 1
        seq = annotation.raw_sequences[qname]
        if r2:
            header = (
                f">{n}__CB:{CB}__UMI:{UMI}__"
                + f"sig:{sig_str}__cDNA:{cDNA_start}-{cDNA_end}__oli:{start_oli}-{end_oli}__"
                + f"L_read={len(seq)}__L_cDNA={cDNA_end - cDNA_start}"
            )
            # print(header[: (254 - 14)] + " 1:N:0:TCCTGAGC" + f"\n{cDNA}")
            print(header[: (254 - 14)] + " 1:N:0:TCCTGAGC" + f"\n{r2}")

    if len(known):
        nd = len(CB_detect)
        ov = len(CB_detect & known)
        ctrl = len(CB_detect & rev)
        args.logger.info(
            f"{ov} / {nd} detected barcodes overlap ({ov/nd * 100:.3f}%). "
            f"Reverse BC control {ctrl/nd * 100:.3f}%"
        )


def prepare_extract_parser(subparsers):
    parser = subparsers.add_parser(
        "extract", help="extract sequences from the long read"
    )
    parser.add_argument(
        "fname",
        default=None,
        help="file with pacbio reads (FASTQ or BAM format)",
    )
    parser.add_argument(
        "--barcode-after",
        type=str,
        default="bead_start",
        help="name of anchor match after which the barcodes follow (default='bead_start')",
    )
    parser.add_argument(
        "--cDNA-after",
        type=str,
        default="bead_start",
        help="excise sequence between this oligo and the last oligo match (if it exists)",
    )
    parser.add_argument(
        "--sig-include",
        type=str,
        default="",
        help="extract cDNA only from long reads INCLUDING this substring in the signature",
    )
    parser.add_argument(
        "--sig-exclude",
        type=str,
        default="",
        help="extract cDNA only from long reads EXCLUDING this substring in the signature",
    )
    parser.add_argument(
        "--distal",
        default=150,
        type=int,
        help="number of nt from end of sequence that are considered for oligo matches when extracting cDNA",
    )
    parser.add_argument(
        "--top-barcodes",
        type=str,
        default="",
        help="path to text file with known barcodes (e.g. from Illumina)",
    )
    parser.set_defaults(func=cmd_extract)
    return parser


## subcommand 'overview'
def cmd_overview(args):
    args = setup_namespace(args, need_sample_name=False)

    dfs = []
    for fname in list(args.reports):
        print(f"loading {fname}")
        df = pd.read_csv(fname, sep="\t", index_col=0)
        dfs.append(df)

    df = pd.concat(dfs, ignore_index=True)

    dfsig = df[["sample", "signature"]].drop_duplicates().set_index("sample")
    print(dfsig)
    df = df[["name", "count", "sample"]].pivot(
        index="sample", columns="name", values="count"
    )  # .query("kind == 'overview'")
    df = df.fillna(0).div(df["n_total"], axis=0)
    df *= 100  # we'd like percentages

    repriming = [
        "TSO,TSO_RC",
        "dN-SMRT,dN-SMRT_RC",
    ]
    concatenation = [c for c in df.columns if "+" in c]
    bead = ["complete", "only_bead_start", "missing_polyT", "missing_OP1"][::-1]

    # # avoid crash if columns are missing
    for r in repriming + concatenation + bead:
        if r not in df.columns:
            df[r] = 0

    # print(df)
    # print(f"concat columns {concatenation}")
    # print(f"bead columns {bead}")
    df["reprimed"] = df[repriming].sum(axis=1)
    # df["complete"] = np.nan_to_num(df["complete"], nan=0.0)
    df["concat"] = df[concatenation].sum(axis=1)
    df["bead_related"] = np.nan_to_num(df[bead].sum(axis=1), nan=0.0)

    def ds(row):
        sig = dfsig.loc[row.name]
        res = (sig == "dropseq") * row.complete
        return res

    df["dropseq"] = df[["complete"]].apply(ds, axis=1)
    # df["complete"] = df[""]
    df["incomplete"] = df["bead-related"] - df["complete"] - df["dropseq"]
    df["non_bead"] = 100 - df["bead-related"]
    df["fidelity"] = 100 * df["complete"] / df["bead-related"]
    df = df.fillna(0)
    # print(df)

    df = df.reset_index().sort_values("sample")
    store_results(
        df,
        args.output,
        "overview.csv",
        float_format="%.2f",
        logger=args.logger,
    )
    report.overview_plots(
        df,
        path=args.output,
        bead=bead,
        concatenation=concatenation,
        repriming=repriming,
    )

    return df


def prepare_overview_parser(subparsers):
    parser = subparsers.add_parser(
        "overview", help="make overview plots across samples"
    )
    parser.add_argument("reports", nargs="+", help="sample reports to aggregate")
    parser.add_argument(
        "--output",
        default=".",
        help="path for detailed report PDFs and CSV overview table",
    )
    parser.set_defaults(func=cmd_overview)
    return parser


## main commandline interface
def prepare_toplevel_parser():
    parser = argparse.ArgumentParser(prog="longread")
    parser.add_argument(
        "--debug", default=False, action="store_true", help="activate debug output"
    )
    parser.add_argument(
        "--sample",
        help="sample name (default=autodetect from fname)",
        default=None,
    )
    parser.add_argument(
        "--config",
        help="YAML describing expected longread signature (default=./longread.yaml if detected or built-in)",
        default="longread.yaml",
    )
    parser.add_argument(
        "--signature",
        help="expected long read signature (e.g. dropseq, chromium, noUMI, withUMI,...)",
        default=None,
    )
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
        "--examples-out",
        help="path to store annotated example read sequences",
        default="./examples/",
    )
    parser.add_argument(
        "--stats-out",
        help="path to store statistics (pandas dataframes)",
        default="./stats/",
    )
    parser.add_argument(
        "--report-out", help="path to render graphical reports in", default="./reports/"
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

    ## sub-parser setup ##
    sp = parser.add_subparsers(help="sub-command help", dest="subcmd")
    prepare_align_parser(sp)
    prepare_annotate_parser(sp)
    prepare_edits_parser(sp)
    prepare_extract_parser(sp)
    prepare_report_parser(sp)
    prepare_overview_parser(sp)

    parser.set_defaults(func=lambda args: parser.print_help())

    return parser


def cmdline():
    """
    Commandline interface to handle long-read (PacBio, nanopore, ...)
    annotation and analysis.

    :return: depends on sub-command used
    :rtype: depends on sub-command used
    """
    import logging

    parser = prepare_toplevel_parser()
    args = parser.parse_args()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    return args.func(args)


if __name__ == "__main__":
    result = cmdline()
