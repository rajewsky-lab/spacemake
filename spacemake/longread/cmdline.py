__author__ = ["Marvin Jens"]
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
        if not hasattr(args.fname):
            raise ValueError("please provide a sample-identifier via --sample")

        sample_name, _ = os.path.splitext(os.path.basename(args.fname))
        logging.info(f"auto-detected sample_name={sample_name}")
    else:
        sample_name = args.sample
    return sample_name.replace(".stats", "")


def store_results(df, path, fname, logger, **kw):
    fpath = os.path.join(util.ensure_path(path), fname)
    df.to_csv(fpath, sep="\t", **kw)
    logger.info(f"storing {len(df)} rows to '{fpath}'")
    return fpath


def load_results(path, fname, logger, **kw):
    fpath = os.path.join(path, fname)
    df = pd.read_csv(fpath, sep="\t", comment="#", **kw).fillna("other")
    logger.info(f"loaded {len(df)} rows from '{fpath}'")
    return df


class SignatureDB:
    def __init__(self, intact={}, related={}, color={}, label={}, prio={}):
        self.intact = intact
        self.related = related
        self.color = color
        self.label = label
        self.prio = prio

    @classmethod
    def from_YAML(cls, fname="samples.yaml"):
        import yaml

        logger = logging.getLogger("spacemake.longread.SampleDB.from_YAML")
        logger.info(f"reading longread signature definitions from '{fname}'")

        groups = yaml.load(open(fname), Loader=yaml.SafeLoader)
        grp = groups["signatures"]
        default = grp[groups["default"]]
        # print("default THING", default)

        intact_lkup = defaultdict(lambda: default["intact_bead"])
        brelated_lkup = defaultdict(lambda: default["bead_related"])
        color_lkup = defaultdict(lambda: default["color"])
        label_lkup = defaultdict(lambda: default["label"])
        prio_lkup = defaultdict(lambda: default["prio"])
        for name, d in groups["signatures"].items():
            # print(f"name={name} d={d}")
            intact_lkup[name] = d["intact_bead"]
            brelated_lkup[name] = d["bead_related"]
            color_lkup[name] = d["color"]
            label_lkup[name] = d["label"]
            prio_lkup[name] = d["prio"]

        # print(color_lkup)
        # print("DEFAULT COLOR", color_lkup["TEST"])
        logger.info(
            f"found {len(intact_lkup)} signature definitions: {sorted(intact_lkup.keys())}."
        )
        return cls(intact_lkup, brelated_lkup, color_lkup, label_lkup, prio_lkup)

    def sort_samples(self, samples, signatures):
        return sorted(
            zip(samples, signatures), key=lambda x: (self.prio.get(x[1], np.inf), x[0])
        )


def get_signature_db(args):
    if os.access(args.config, os.R_OK):
        cfg = args.config
    else:
        cfg = os.path.join(os.path.dirname(__file__), "../data/config/longread.yaml")

    return SignatureDB.from_YAML(cfg)


def setup_namespace(args):
    "perform initialization code common to more than one sub-command"
    d = vars(args)
    sample = detect_sample(args)
    d["sample_name"] = sample

    db = get_signature_db(args)
    d["signature_db"] = db

    intact_bead = db.intact[args.signature]
    d["intact_bead"] = intact_bead
    d["sig_intact"] = intact_bead.split(",")
    d["sig_core"], d["sig_core_order"] = util.process_intact_signature(intact_bead)

    d["bead_related"] = db.related[args.signature]
    # bead_related = args.bead_related
    # if bead_related is None:
    #     sig_intact = tuple(args.intact_bead.split(","))
    #     bead_related = sig_intact[0]

    d["logger"] = logging.getLogger(f"spacemake.longread.cmd.{args.subcmd}")
    d["blocks"] = util.load_oligos(args.blocks)

    from types import SimpleNamespace

    return SimpleNamespace(**d)


def aln_main(args):
    args = setup_namespace(args)

    cache.fill_caches(
        args.fname,
        args.sample_name,
        args.blocks,
        path=util.ensure_path(args.cache),
        n_proc=args.parallel,
    )
    df = cache.annotate(
        args.fname, args.sample_name, args.blocks, path=util.ensure_path(args.cache)
    )
    df.to_csv(
        os.path.join(
            util.ensure_path(args.annotation_out), f"{args.sample_name}.annotation.tsv"
        ),
        sep="\t",
        index=False,
    )


def ann_main(args):
    args = setup_namespace(args)
    args.logger.info(
        f"analyzing {args.sample_name} for intact signature {args.sig_intact}"
    )
    annotation = ann.AnnotatedSequences(
        args.fname,
        os.path.join(args.annotation_out, f"{args.sample_name}.annotation.tsv"),
        args.sample_name,
        args.blocks,
        min_score=args.min_score,
        orient_by=args.bead_related,
    )
    n_total = len(annotation.raw_sequences)
    logging.info(
        f"total number of reads in {args.fname} ({args.sample_name}) is {n_total}"
    )

    sig_counts, n_concat, n_reprimed = annotation.count_signatures()
    util.count_dict_out(sig_counts, "common signatures", total=n_total)
    print(
        f"n_concat={n_concat} ({100.0 * n_concat/n_total:.2f}%) "
        f"n_reprimed={n_reprimed} ({100.0 * n_reprimed/n_total:.2f}%)"
    )
    df_sig = util.count_dict_to_df(sig_counts, "signatures", n_total=n_total)

    concat_counts, n_occurrences = annotation.count_concatenations()
    util.count_dict_out(concat_counts, "concatenations", total=n_occurrences)
    df_concat = util.count_dict_to_df(
        concat_counts, "concatenations", n_total=n_occurrences
    )

    reprimed_counts = annotation.count_repriming()
    util.count_dict_out(reprimed_counts, "repriming", total=n_total)
    df_reprimed = util.count_dict_to_df(reprimed_counts, "repriming", n_total=n_total)

    partial_counts, prefixes, suffixes, pT_counts = annotation.completeness(
        args.sig_core, polyT=args.polyT
    )

    partial_counts_simple, _ = util.count_dict_collapse_misc(
        partial_counts, sig_intact=args.sig_intact, total=n_total, misc_thresh=0.00001
    )
    util.count_dict_out(partial_counts_simple, "completeness", total=n_total)

    df_comp = util.count_dict_to_df(
        partial_counts_simple, kind="bead_complete", n_total=n_total
    ).sort_values("name")

    df_pT = util.count_dict_to_df(pT_counts, kind="polyT_after", n_total=n_total)

    df = pd.concat([df_sig, df_concat, df_reprimed, df_comp, df_pT])
    store_results(df, args.stats_out, f"{args.sample_name}.stats.tsv", args.logger)

    # TODO: prefix/suffix counts add up to > 100%. Needs fix
    util.count_dict_out(pT_counts, "polyT after", total=n_total)
    util.count_dict_out(prefixes, "prefixes", total=n_total)
    util.count_dict_out(suffixes, "suffixes", total=n_total)

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
    eo_fname = os.path.join(
        util.ensure_path(args.examples_out), f"{args.sample_name}.txt"
    )

    with open(eo_fname, "wt") as eo:
        for signame, sigcount in sorted(sig_counts.items(), key=lambda x: -x[1]):
            try:
                qname, _, _ = next(
                    annotation.filter_signatures(tuple(signame.split(",")))
                )
                eo.write(f"# {signame} n={sigcount}\n{annotation.fmt(qname)}\n")
            except StopIteration:
                args.logger.warning(
                    f"unable to find any reads with signature {signame}"
                )


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


def rep_main(args):
    args = setup_namespace(args)
    args.logger.info(f"generating report plots for '{args.sample_name}'")
    args.logger.debug(f"'bead-related' if we detect: '{args.bead_related}'")

    df = load_results(args.stats_out, f"{args.sample_name}.stats.tsv", args.logger)
    sig_counts = util.count_dict_from_df(df, "signatures")
    util.count_dict_out(
        sig_counts, "signatures", misc_thresh=0.05, total=sig_counts["n_total"]
    )
    n_total = sig_counts["n_total"]
    args.logger.info(f"n_total={n_total}")

    # disentangle bead-related from other signatures
    ov_counts, bead_counts, found_part_counts, core_signature = util.digest_signatures(
        sig_counts,
        args.bead_related,
        args.intact_bead,
        # prefixes=args.prefixes,
        # suffixes=args.suffixes,
    )
    # group low abundance signatures into 'misc' for overview donut plot
    ov_simple, _ = util.count_dict_out(
        ov_counts, "overview counts", misc_thresh=0.01, total=sig_counts["n_total"]
    )
    # group low abundance signatures into 'misc' for bead completeness donut plot
    bead_simple, _ = util.count_dict_out(
        bead_counts, "bead counts", misc_thresh=0.0, total=ov_counts["bead-related"]
    )
    # store the donut plot values in ...report.tsv for cross-sample overviews and such
    df_ov = util.count_dict_to_df(ov_simple, kind="overview")
    df_bead = util.count_dict_to_df(bead_simple, kind="bead_related")
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
            util.ensure_path(args.report_out), f"{args.sample_name}.donuts.pdf"
        ),
        suptitle=args.sample_name,
    )

    # render the pre-computed position/score/length histograms as a plot
    report.plot_histograms(
        df=load_results(
            args.stats_out, f"{args.sample_name}.intact_parts.tsv", args.logger
        ),
        fname=os.path.join(
            util.ensure_path(args.report_out), f"{args.sample_name}.hists.pdf"
        ),
        n_total=n_total,
        parts=args.intact_bead.split(","),
    )


def main_edits(args):
    args = setup_namespace(args)

    annotation = ann.AnnotatedSequences(
        args.fname,
        os.path.join(args.annotation_out, f"{args.sample_name}.annotation.tsv"),
        args.sample_name,
        args.blocks,
        min_score=args.min_score,
    )
    n_total = len(annotation.raw_sequences)
    logging.info(
        f"total number of reads in {args.fname} ({args.sample_name}) is {n_total}"
    )

    data = []
    for part in args.sig_intact:
        qmatches = annotation.query_oligo_occurrences(part)
        if len(qmatches[0]) > args.n_samples:
            qmatches = ann.subsample(qmatches, n=args.n_samples)

        nmatch = len(qmatches[0])
        m, ed = ann.align_stats(annotation, args.blocks[part], qmatches)
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


def main_overview(args):
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
    if args.csv_out:
        df.to_csv(args.csv_out, float_format="%.2f", sep="\t")

    def clean(txt):
        txt = os.path.basename(txt)
        t = (
            txt.replace("source/", "")
            .replace("sts_", "")
            .replace("pb_", "")
            .replace("ds_", "")
            .replace(".fq", "")
            .replace(".bam", "")
            .replace("lima.", "")
        )

        if t.count("_") > 1:
            t = "_".join(t.split("_")[:2])

        return t

    df = df.reset_index()
    # df = df.sort_values('bead_related')
    df = df.sort_values("sample")

    def guess_rRNA_file(path):
        # print("guessrRNA raw path", path)
        name = os.path.basename(path).replace(".summary", ".rRNA")

        if args.rRNA_same_place:
            place = os.path.dirname(path)
        else:
            place = args.rRNA

        return [
            os.path.join(place, name.replace(".fq", ".txt")),
            os.path.join(place, name.replace(".fq", ".txt")).replace(
                ".rRNA.tsv", ".txt"
            ),
            os.path.join(place, name.replace(".fq", ".txt")).replace(
                ".rRNA.tsv", ".rRNA.txt"
            ),
            os.path.join(place, name.replace(".bam", ".txt").replace("lima.", "")),
            os.path.join(
                place, name.replace(".bam", ".txt").replace("lima.", "")
            ).replace(".rRNA.tsv", ".txt"),
            os.path.join(
                place, name.replace(".bam", ".txt").replace("lima.", "")
            ).replace(".rRNA.tsv", ".rRNA.txt"),
        ]

    # rRNA_fracs = []
    # for row in df[["stats_file", "N_reads"]].itertuples():
    #     rcount = np.nan
    #     for fname in guess_rRNA_file(row.stats_file):
    #         print(fname)
    #         try:
    #             rcount = int(open(fname).read())
    #         except (FileNotFoundError, ValueError):
    #             pass
    #         else:
    #             break
    #     if rcount == np.nan:
    #         raise ValueError

    #     rRNA_fracs.append(100.0 * rcount / row.N_reads)

    # df["rRNA"] = rRNA_fracs
    # print(df[['qfa', 'rRNA']])

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 8})

    def make_bars(
        ax, df, kinds, labels, cmap=plt.get_cmap("tab10"), w=0.9, colors=None
    ):
        n = len(kinds)
        if colors is None:
            colors = cmap(np.linspace(0, 1, n))

        x = np.arange(len(df)) - w / 2.0
        y0 = np.zeros(len(x), dtype=float)
        for kind, label, color in zip(kinds, labels, colors):
            y = np.nan_to_num(df[kind], nan=0.0)
            # print(kind)
            # print(y)
            ax.bar(x, y, bottom=y0, label=label, width=w, color=color)
            y0 += y

        ax.set_ylabel("fraction of library")
        ax.set_xticks(x)
        labels = df["sample"]  # [clean(fq) for fq in df['qfa']]
        ax.set_xticklabels(labels, rotation=90)
        ax.set_ylim(0, 100)

    marie = [
        "non_bead",
        "incomplete",
        "dropseq",
        "complete",
    ]
    marie_colors = ["gray", "royalblue", "green", "gold"]

    w = max(8 / 25.0 * len(df), 5)
    if args.multi_page:
        pdf = PdfPages(args.breakdown)
        fig, ax1 = plt.subplots(1, figsize=(w, 8))
    else:
        fig, (ax1, ax2) = plt.subplots(2, figsize=(w, 8), sharex=True)

    make_bars(
        ax1,
        df,
        marie,
        labels=marie,
        colors=marie_colors,
    )
    ax1.legend(title="Marie-stats", ncol=len(marie))
    if args.multi_page:
        fig.tight_layout()
        pdf.savefig()
        plt.close()
        fig, ax2 = plt.subplots(1, figsize=(w, 4))

    make_bars(ax2, df, ["fidelity"], labels=["bead fidelity"])
    ax2.set_ylabel("bead fidelity")
    if args.multi_page:
        fig.tight_layout()
        pdf.savefig()
        pdf.close()
    else:
        fig.tight_layout()
        plt.savefig(args.breakdown)

    plt.close()

    if args.multi_page:
        pdf = PdfPages(args.output)
        fig, ax1 = plt.subplots(1, figsize=(w, 4))
    else:
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize=(w, 12), sharex=True)

    # print("bead related", bead)
    make_bars(ax1, df, bead, labels=bead)
    ax1.legend(title="bead-related", ncol=len(bead))
    if args.multi_page:
        fig.tight_layout()
        pdf.savefig()
        plt.close()
        fig, ax2 = plt.subplots(1, figsize=(w, 4))

    # print("repriming events", repriming)
    make_bars(
        ax2,
        df,
        repriming,
        labels=[r.split(",")[0] for r in repriming],
        cmap=plt.get_cmap("tab20c"),
    )
    ax2.legend(title="repriming", ncol=len(repriming))
    if args.multi_page:
        fig.tight_layout()
        pdf.savefig()
        plt.close()
        fig, ax3 = plt.subplots(1, figsize=(w, 4))

    # print("concat events", concatenation)
    make_bars(ax3, df, concatenation, labels=concatenation, cmap=plt.get_cmap("tab20b"))
    ax3.legend(title="concatamers", ncol=len(concatenation))
    if args.multi_page:
        fig.tight_layout()
        pdf.savefig()
        plt.close()
        fig, ax4 = plt.subplots(1, figsize=(w, 4))

    # make_bars(
    #     ax4,
    #     df,
    #     [
    #         "rRNA",
    #     ],
    #     labels=["rRNA"],
    #     cmap=plt.get_cmap("tab20c"),
    # )
    ax4.legend(title="human rRNA", ncol=1)
    if args.multi_page:
        fig.tight_layout()
        pdf.savefig()
        pdf.close()
    else:
        fig.tight_layout()
        plt.savefig(args.output)

    plt.close()
    return df


def prepare_parser():
    parser = argparse.ArgumentParser(prog="longread")
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
        help="expected long read signature (e.g. dropseq, noUMI, withUMI,...)",
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
    parser.add_argument(
        "--debug", default=False, action="store_true", help="activate debug output"
    )

    ## sub-parser setup ##
    subparsers = parser.add_subparsers(help="sub-command help", dest="subcmd")
    aln_parser = subparsers.add_parser(
        "align", help="align PacBio reads against oligos"
    )
    aln_parser.add_argument(
        "fname",
        default=None,
        help="file with pacbio reads (FASTQ or BAM format)",
    )
    ann_parser = subparsers.add_parser(
        "annotate", help="create annotation from detected oligo matches"
    )
    ann_parser.add_argument(
        "fname",
        default=None,
        help="file with pacbio reads (FASTQ or BAM format)",
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
    ed_parser.add_argument(
        "fname",
        default=None,
        help="file with pacbio reads (FASTQ or BAM format)",
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
    ov_parser.add_argument("reports", nargs="+", help="sample reports to aggregate")
    ov_parser.add_argument(
        "--output", default="pb_overview.pdf", help="path/name of detailed report PDF"
    )
    ov_parser.add_argument(
        "--csv-out", default="all_pb_stats.csv", help="path/name of detailed report PDF"
    )
    ov_parser.add_argument(
        "--breakdown",
        default="bead_overview.pdf",
        help="path/name of bead report (Marie style) PDF",
    )
    # ov_parser.add_argument(
    #     "--rRNA",
    #     default="rRNA/",
    #     help="path to search for rRNA counts corresponding to samples",
    # )
    # ov_parser.add_argument(
    #     "--rRNA-same-place",
    #     default=False,
    #     action="store_true",
    #     help="If set, look for rRNA txt file with same sample name in same directory",
    # )
    ov_parser.add_argument(
        "--multi-page",
        default=False,
        action="store_true",
        help="If set, generate multiple PDF pages instead of subplots",
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
