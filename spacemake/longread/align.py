import os
import sys
import re
import pandas as pd
import numpy as np

# import pysam
import multiprocessing as mp
from collections import defaultdict, OrderedDict, Counter
from more_itertools import grouper
from Bio import pairwise2
from Bio import SeqIO
from spacemake.util import rev_comp
from spacemake.util import fasta_chunks
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

blocks = OrderedDict()
for fa_id, seq in fasta_chunks(open("blocks_combv2.fa")):
    blocks[fa_id] = seq
    blocks[fa_id + "_RC"] = rev_comp(seq)

block_names = list(blocks.keys())


def read_fq(fname):
    for name, seq, _, qual in grouper(open(fname), 4):
        yield name.rstrip()[1:], seq.rstrip(), qual.rstrip()


# def orient(seq, fw=blocks['P5'][10:25], rc=blocks['N70X_RC'][31:46]):
#     start = seq[:50]
#     f = fw in start
#     r = rc in start
#     if f and not r:
#         return seq
#     elif r and not f:
#         return rev_comp(seq)
#     else:
#         print(start, f, r)
#         return None


# sample = "sts_067_1"
sample = sys.argv[1]
# base_path = (
#     f"/data/rajewsky/projects/combbeads/three_segments/pacbio/sts_083_{sample}.fq"
# )
base_path = (
    f"/data/rajewsky/projects/combbeads/three_segments/pacbio/sts_067_{sample}.fq"
)

# base_path = f"/data/rajewsky/projects/slide_seq/projects/sts_067/raw_data/pacbio/{sample}.fastq"
# oriented = []
# names = []
# for name, seq, qual in read_fq(base_path):
#     s = orient(seq)
#     if s is not None:
#         oriented.append(s)

# lengths = np.array([len(s) for s in oriented])
# L_max = lengths.max()

aln_caches = {}
import logging

logging.basicConfig(level=logging.DEBUG)
import cache


# # Run this block to compute all alignments for later re-use
# o = sys.argv[2]
# aln = cache.CachedAlignments(sample, o, blocks[o])

# # for name, seq, qual in read_fq(base_path):
# for name, seq, qual in read_fq(base_path):
#     aln.align_or_load(name, seq)

# sys.exit(0)


def annotate():
    multi = cache.MultiAlignments(sample, blocks)
    for o in ["SMART", "P5", "OP1", "OP2", "OP3", "pT", "N70X"]:
        aln_caches[o] = cache.CachedAlignments(sample, o, blocks[o])

    print("qname\tL\toligo\tstart\tend\tscore")
    for name, seq, qual in read_fq(base_path):
        for start, end, label, score in multi.annotate(name, seq):
            # pass
            print(f"{name}\t{len(seq)}\t{label}\t{start}\t{end}\t{score}")


# # Run this block next, to create the annotation flat files
# annotate()
# sys.exit(0)


class AnnotatedSequences:
    def __init__(
        self,
        sample_name,
        blocks,
        path=base_path,
    ):
        self.sample_name = sample_name
        self.path = path
        self.logger = logging.getLogger("AnnotatedSequences")
        self.raw_sequences = self.load_raw_sequences()
        self.ann_db = self.load_annotation()
        self.signatures = self.extract_signatures_and_orient()

    def load_raw_sequences(self):
        self.logger.info(f"loading raw sequences for {self.sample_name}")
        raw_sequences = {}
        for qname, seq, qual in read_fq(self.path.format(sample_name=self.sample_name)):
            raw_sequences[qname] = seq
        return raw_sequences

    def load_annotation(self):
        self.logger.info(f"loading oligo annotation for {self.sample_name}")
        df = pd.read_csv(f"{self.sample_name}.annotation.tsv", sep="\t")
        qdata = {}
        for qname, grp in df.groupby("qname"):
            qdata[qname] = (
                tuple(grp["oligo"]),
                grp["start"].values,
                grp["end"].values,
                grp["score"].values,
            )
        return qdata

    def extract_signatures_and_orient(self):
        self.logger.info(f"extracting signatures for {self.sample_name}")
        signatures = {}
        for qname, (names, starts, ends, scores) in self.ann_db.items():
            seq = self.raw_sequences[qname]
            L = len(seq)
            if names[0].endswith("_RC"):
                # reverse complement the whole thing
                names = tuple([(n + "_RC").replace("_RC_RC", "") for n in names[::-1]])
                scores = scores[::-1]
                starts, ends = L - ends[::-1], L - starts[::-1]
                seq = rev_comp(seq)
                # store the RC'ed data
                self.raw_sequences[qname] = seq
                self.ann_db[qname] = (names, starts, ends, scores)

            signatures[qname] = names

        return signatures

    def filter_signatures(self, qsig):
        for qname, sig in self.signatures.items():
            if sig == qsig:
                yield qname

    def query_dimensions(self, qsig):
        qnames = []
        L = []
        starts = []
        ends = []
        scores = []

        for qname in self.filter_signatures(qsig):
            qnames.append(qname)
            L.append(len(self.raw_sequences[qname]))
            names, s, e, sc = self.ann_db[qname]
            starts.append(s)
            ends.append(e)
            scores.append(sc)

        return (
            np.array(qnames),
            np.array(L),
            np.array(starts),
            np.array(ends),
            np.array(scores),
        )

    def query_oligo_occurrences(self, oligo):
        qnames = []
        L = []
        starts = []
        ends = []
        scores = []

        for qname, sig in self.signatures.items():
            if oligo in sig:
                i = list(sig).index(oligo)

                qnames.append(qname)
                L.append(len(self.raw_sequences[qname]))
                names, s, e, sc = self.ann_db[qname]
                starts.append(s[i])
                ends.append(e[i])
                scores.append(sc[i])

        return (
            np.array(qnames),
            np.array(starts),
            np.array(ends),
            np.array(scores),
            np.array(L),
        )

    def fmt(self, qname):
        def render_label(label, start, end):
            return list(f"|{label.center(end - start - 2)}|")

        #
        seq = self.raw_sequences[qname]
        L = len(seq)
        buf = [" "] * L
        for row in self.ann_db[qname].itertuples():
            buf[row.start : row.end + 1] = render_label(
                row.oligo, row.start, row.end + 1
            )
        #
        return f"{seq}\n{''.join(buf)}"


sig_intact = ("P5", "SMART", "OP1", "OP2", "OP3", "pT", "N70X")
S = AnnotatedSequences(
    sys.argv[1],
    blocks,
    path=base_path,
)
qintact, qL, qstarts, qends, qsores = S.query_dimensions(sig_intact)


# estimate extension efficiencies
# n_reads = 0
stages = np.zeros(5, dtype=float)
pT = np.zeros(4)
for qname, sig in S.signatures.items():
    stages[0] += 1
    # print(sig)
    l = list(sig)
    if "SMART" in sig:
        l = l[sig.index("SMART") + 1 :]
        stages[1] += 1
        try:
            block = l.pop(0)
            if block == "pT":
                pT[0] += 1
            elif block == "OP1":
                stages[2] += 1
                block = l.pop(0)
                if block == "pT":
                    pT[1] += 1
                elif block == "OP2":
                    stages[3] += 1
                    block = l.pop(0)
                    if block == "pT":
                        pT[2] += 1
                    elif block == "OP3":
                        stages[4] += 1
                        block = l.pop(0)
                        if block == "pT":
                            pT[3] += 1
        except IndexError:
            # this just happens when l is empty
            pass

# print(stages)
# stages /= n_reads
freqs = stages[1:] / stages[0]
efficiencies = stages[2:] / stages[1:-1]
pT_rates = pT / stages[0]
# print(stages)
print("## estimated extension efficiencies")
print(
    "sample\tf_SMART\tf_OP1\tf_OP2\tf_OP3\tr_seg1\tr_seg2\tr_seg3\tpT_0\tpT_1\tpT_2\tpT_3"
)
print(
    f"{sample}\t"
    + "\t".join([f"{e:.3f}" for e in np.concatenate((freqs, efficiencies, pT_rates))])
)
sys.exit(0)


def analyze_recurrence(AS, oligo):
    # collect start positions of additional oligo occurrences in "intact" reads
    extra_dstart = []
    extra_dend = []
    extra_scores = []
    extra_L = []
    for qname, sig in AS.signatures.items():
        if sig[:5] == sig_intact[:5]:
            if oligo in sig[5:]:
                i = 5 + list(sig[5:]).index(oligo)
                names, starts, ends, scores = AS.ann_db[qname]
                L = len(AS.raw_sequences[qname])
                extra_dstart.append(starts[i])
                extra_dend.append(L - ends[i])
                extra_scores.append(scores[i])
                extra_L.append(L)
    #
    return (
        np.array(extra_dstart),
        np.array(extra_dend),
        np.array(extra_scores),
        np.array(extra_L),
    )


def align_stats(ann, oname, qmatches, pad=1):
    oseq = blocks[oname]
    n = 0
    matches = np.zeros(len(oseq), dtype=float)
    edits = defaultdict(lambda: defaultdict(int))
    for qn, s, e, score, l in zip(*qmatches):
        s = max(0, s - pad)
        e = min(l, e + pad)
        n += 1
        rseq = ann.raw_sequences[qn][s:e]
        res = cache.align(rseq, oseq, min_score=0)
        if not res:
            continue
        aln = res[0]
        aln_str = pairwise2.format_alignment(*aln, full_sequences=True).split("\n")
        i = 0
        j = 0
        # print(aln_str)
        for r, m, q in zip(aln_str[0], aln_str[1], aln_str[2]):
            # print(r, m, q, i, len(oseq))
            if q == " ":
                # outside of oligo match range
                continue
            if m == "|":
                # print(matches[i])
                matches[i] += 1
            else:
                edits[i][q + r] += 1
            if q != "-":
                i += 1
    return matches / n, edits


def select(
    qmatches,
    from_start_min=None,
    from_start_max=None,
    from_end_min=None,
    from_end_max=None,
    min_score=None,
):
    qnames, starts, ends, scores, L = qmatches
    mask = np.ones(len(qnames), dtype=np.bool)
    if from_start_min is not None:
        mask &= starts >= from_start_min
    if from_start_max is not None:
        mask &= starts <= from_start_max
    fe = L - ends
    if from_end_min is not None:
        mask &= fe >= from_end_min
    if from_end_max is not None:
        mask &= fe <= from_end_max
    if min_score is not None:
        mask &= scores >= min_score
    return (qnames[mask], starts[mask], ends[mask], scores[mask], L[mask])


def subsample(qmatches, n=100):
    data = np.array(qmatches, dtype=np.object)
    n_cols, N = data.shape
    I = np.random.choice(N, size=n, replace=False)
    return data[:, I]


def plot_edits_heatmap(ax, oname, edits):
    oseq = blocks[oname]
    N = len(oseq)
    counts = np.zeros((N, 5))
    x = np.arange(N)
    for i in x:
        for j, nn in enumerate("ACGT-"):
            if nn == oseq[i]:
                counts[i, j] = np.nan
            counts[i, j] = edits[i][oseq[i] + nn]
    im = ax.imshow(counts.T, interpolation="none", cmap="viridis")
    plt.colorbar(im, ax=ax, shrink=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(list(oseq))
    ax.set_yticks(range(5))
    ax.set_yticklabels(list("ACGT-"))


fig, (ax1, ax2) = plt.subplots(2, figsize=(5, 8), sharex=True)
# plot occurrences of OP1,OP2 regardless of other matches
exp1 = 1 + len(blocks["P5"] + blocks["SMART"]) + 9.5 + 2
exp2 = 1 + len(blocks["P5"] + blocks["SMART"] + blocks["OP1"]) + 9.5 + 8
smart = op1 = S.query_oligo_occurrences("SMART")
op1 = S.query_oligo_occurrences("OP1")
op2 = S.query_oligo_occurrences("OP2")
op3 = S.query_oligo_occurrences("OP3")
pT = S.query_oligo_occurrences("pT")
op1rc = S.query_oligo_occurrences("OP1_RC")
op2rc = S.query_oligo_occurrences("OP2_RC")
op3rc = S.query_oligo_occurrences("OP3_RC")

r2_bound = 43 + len(blocks["N70X"])


def prepare_analyses(min_score=None):
    analyses = [
        # # Opseq1
        # (
        #     "OP1",
        #     "correct",
        #     select(
        #         op1,
        #         from_start_min=exp1 - 5,
        #         from_start_max=exp1 + 5,
        #         min_score=min_score,
        #     ),
        # ),
        # (
        #     "OP1",
        #     "prox",
        #     select(
        #         op1,
        #         from_start_min=exp2 + len(blocks["OP2"]),
        #         from_start_max=exp2 + len(blocks["OP2"]) + 100,
        #         from_end_min=r2_bound,
        #         min_score=min_score,
        #     ),
        # ),
        # (
        #     "OP1",
        #     "internal",
        #     select(
        #         op1,
        #         from_start_min=exp2 + len(blocks["OP2"]) + 100,
        #         from_end_min=r2_bound,
        #         min_score=min_score,
        #     ),
        # ),
        # (
        #     "OP1",
        #     "read2",
        #     select(op1, from_end_max=r2_bound, min_score=min_score),
        # ),
        # # Opseq2
        # (
        #     "OP2",
        #     "correct",
        #     select(
        #         op2,
        #         from_start_min=exp2 - 5,
        #         from_start_max=exp2 + 5,
        #         min_score=min_score,
        #     ),
        # ),
        # (
        #     "OP2",
        #     "prox",
        #     select(
        #         op2,
        #         from_start_min=exp2 + len(blocks["OP2"]) + 5,
        #         from_start_max=exp2 + len(blocks["OP2"]) + 105,
        #         from_end_min=r2_bound,
        #         min_score=min_score,
        #     ),
        # ),
        # (
        #     "OP2",
        #     "internal",
        #     select(
        #         op2,
        #         from_start_min=exp2 + len(blocks["OP2"]) + 105,
        #         from_end_min=r2_bound,
        #         min_score=min_score,
        #     ),
        # ),
        # (
        #     "OP2",
        #     "read2",
        #     select(op2, from_end_max=r2_bound, min_score=min_score),
        # ),
        # SMART
        ("SMART", "anywhere", select(smart, min_score=min_score)),
        # Opseq1
        ("OP1", "anywhere", select(op1, min_score=min_score)),
        # Opseq2
        ("OP2", "anywhere", select(op2, min_score=min_score)),
        # Opseq3
        ("OP3", "anywhere", select(op3, min_score=min_score)),
        # polyT
        ("pT", "anywhere", select(pT, min_score=min_score)),
        # Opseq1 RC
        ("OP1_RC", "anywhere", select(op1rc, min_score=min_score)),
        # Opseq2 RC
        ("OP2_RC", "anywhere", select(op2rc, min_score=min_score)),
        # Opseq3 RC
        ("OP3_RC", "anywhere", select(op3rc, min_score=min_score)),
    ]
    return analyses


def plot_scores(analyses, suffix=""):
    ### match score statistics
    fig, ax = plt.subplots()
    for oname, region, matches in analyses:
        qnames, starts, ends, scores, L = matches
        x = sorted(scores)
        y = np.linspace(0, 1, len(scores))
        ax.plot(x, y, label=f"{oname} {region} n={len(scores)}")
        ax.set_ylabel("cumulative rel. freq")
        ax.set_xlabel("alignment score")
    #
    ax.legend()
    fig.suptitle(sample)
    # fig.tight_layout()
    fig.savefig(f"scores_{sample}{suffix}.pdf")
    plt.close("all")


def plot_starts(analyses, suffix=""):
    ### start position statistics
    # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 12))
    fig, ax = plt.subplots()
    for oname, region, matches in analyses:
        qnames, starts, ends, scores, L = matches
        x = sorted(starts)
        y = np.linspace(0, 1, len(scores))
        ax.plot(x, y, label=f"{oname} {region} n={len(scores)}")
        ax.set_ylabel("cumulative rel. freq")
        ax.set_xlabel("start positions")
    #
    ax.legend()
    fig.suptitle(sample)
    # fig.tight_layout()
    fig.savefig(f"starts_{sample}{suffix}.pdf")
    plt.close("all")


def plot_readlengths(analyses, suffix=""):
    ### start position statistics
    # fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 12))
    fig, ax = plt.subplots()
    for oname, region, matches in analyses:
        qnames, starts, ends, scores, L = matches
        x = sorted(L)
        y = np.linspace(0, 1, len(x))
        ax.plot(x, y, label=f"{oname} {region} n={len(scores)}")
        ax.set_ylabel("cumulative rel. freq")
        ax.set_xlabel("original read length")
    #
    ax.legend()
    fig.suptitle(sample)
    # fig.tight_layout()
    fig.savefig(f"readlengths_{sample}{suffix}.pdf")
    fig.savefig(f"readlengths_{sample}{suffix}.svg")
    plt.close("all")


def plot_mismatches(analyses, suffix=""):
    fig, ax = plt.subplots()
    fig2, axes = plt.subplots(len(analyses), figsize=(8, 12), sharex=True, sharey=True)
    if len(analyses) == 1:
        axes = [
            axes,
        ]
    #
    for (oname, region, matches), ax2 in zip(analyses, axes):
        if len(matches[0]) > 500:
            matches = subsample(matches, n=500)
        #
        qnames, starts, ends, scores, L = matches
        m, ed = align_stats(S, oname, matches)
        x = np.arange(len(m))
        ax.step(x, m, where="mid", lw=2, label=f"{oname} {region} n={len(qnames)}")
        #
        plot_edits_heatmap(ax2, oname, ed)
        ax2.title.set_text(f"{region} n={len(qnames)}")
        ax2.set_xlabel(f"{oname} sequence")
    #
    ax.set_xticks(x)
    ax.set_ylim(0.5, 1.0)
    ax.legend()
    ax.set_xticklabels(list(blocks[oname]))
    ax.set_ylabel("match frequency")
    ax.set_xlabel(f"{oname} oligo")
    fig.tight_layout()
    fig.savefig(f"match_{oname}_{sample}{suffix}.pdf")
    fig.savefig(f"match_{oname}_{sample}{suffix}.svg")
    fig2.suptitle(f"{sample} {oname} mismatch, deletion statistics by location")
    fig2.savefig(f"edits_{oname}_{sample}{suffix}.pdf")
    fig2.savefig(f"edits_{oname}_{sample}{suffix}.svg")
    plt.close("all")


def plot_startpos_hist(ax, analyses):
    ax.set_yscale("log")
    for oname, region, matches in analyses:
        qnames, starts, ends, scores, L = matches
        x = np.arange(starts.max() + 1)
        pos = np.bincount(starts)
        ax.step(
            x,
            pos,
            where="mid",
            label=f"{oname} {region} n={len(qnames)}",
            lw=1.5,
            alpha=0.75,
        )
    #
    ax.set_ylabel("match count")
    ax.set_xlabel("start position")
    ax.legend()


def plot_fromend_hist(ax, analyses):
    ax.set_yscale("log")
    for oname, region, matches in analyses:
        qnames, starts, ends, scores, L = matches
        fe = L - ends
        x = np.arange(fe.max() + 1)
        pos = np.bincount(fe)
        ax.step(
            x,
            pos,
            where="mid",
            label=f"{oname} {region} n={len(qnames)}",
            lw=1.5,
            alpha=0.75,
        )
    #
    ax.set_ylabel("match count")
    ax.set_xlabel("distance from distal end")
    ax.legend()


analyses = prepare_analyses()
# plot_scores(analyses)
# plot_mismatches(analyses[-2:-1])
# plot_mismatches(analyses[-1:])

# fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
#     2, 2, figsize=(8, 8), sharex=True, sharey=True
# )
# plot_startpos_hist(ax1, analyses[-4:-2])
# plot_fromend_hist(ax2, analyses[-4:-2])
# plot_startpos_hist(ax3, analyses[-2:])
# plot_fromend_hist(ax4, analyses[-2:])
# fig.suptitle(sample)
# fig.savefig(f"aln_pos_{sample}.pdf")
# plt.close("all")

# # using 67_1 as negative control (should not have OP2) 25 is a good cutoff
# analyses = prepare_analyses(min_score=25.0)
# plot_scores(analyses, suffix="_ms")
# plot_starts(analyses, suffix="_ms")
# plot_readlengths(analyses, suffix="_ms")
# plot_mismatches(analyses[:4], suffix="_ms")
# plot_mismatches(analyses[4:8], suffix="_ms")
# plot_mismatches(analyses[-2:-1], suffix="_ms")
# plot_mismatches(analyses[-1:], suffix="_ms")

# fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
#     2, 2, figsize=(8, 8), sharex=True, sharey=True
# )
# plot_startpos_hist(ax1, analyses[-4:-2])
# plot_fromend_hist(ax2, analyses[-4:-2])
# plot_startpos_hist(ax3, analyses[-2:])
# plot_fromend_hist(ax4, analyses[-2:])
# fig.suptitle(sample)
# fig.savefig(f"aln_pos_{sample}_ms.pdf")
# plt.close("all")
plot_scores(analyses)
plot_mismatches(analyses)
maxX = 1000

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
    2, 2, figsize=(8, 8), sharex=True, sharey=True
)
plot_startpos_hist(ax1, analyses[:-3])
plot_fromend_hist(ax2, analyses[:-3])
plot_startpos_hist(ax3, analyses[-3:])
plot_fromend_hist(ax4, analyses[-3:])
fig.suptitle(sample)
ax1.set_xlim(0, maxX)
ax2.set_xlim(0, maxX)
ax3.set_xlim(0, maxX)
ax4.set_xlim(0, maxX)
fig.savefig(f"aln_pos_{sample}.pdf")
plt.close("all")

# using 67_1 as negative control (should not have OP2) 25 is a good cutoff
analyses = prepare_analyses(min_score=25.0)
plot_scores(analyses, suffix="_ms")
plot_starts(analyses, suffix="_ms")
plot_readlengths(analyses, suffix="_ms")
plot_mismatches(analyses, suffix="_ms")
# plot_mismatches(analyses[4:8], suffix="_ms")
# plot_mismatches(analyses[-2:-1], suffix="_ms")
# plot_mismatches(analyses[-1:], suffix="_ms")

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
    2, 2, figsize=(8, 12), sharex=True, sharey=True
)
plot_startpos_hist(ax1, analyses[:-3])
plot_fromend_hist(ax2, analyses[:-3])
plot_startpos_hist(ax3, analyses[-3:])
plot_fromend_hist(ax4, analyses[-3:])
ax1.set_xlim(0, maxX)
ax2.set_xlim(0, maxX)
ax3.set_xlim(0, maxX)
ax4.set_xlim(0, maxX)
fig.suptitle(sample)
fig.savefig(f"aln_pos_{sample}_ms.pdf")
plt.close("all")

### mismatch and deletion statistics for OP1
# op1_correct = select(
#     op1, from_start_min=exp1 - 5, from_start_max=exp1 + 5, min_score=ms_op1
# )
# m_corr, ed_corr = align_stats(S, "OP1", subsample(op1_correct, n=500))

# op1_read2 = select(op1, from_end_max=43 + len(blocks["N70X"]), min_score=ms_op1)
# m_r2, ed_r2 = align_stats(S, "OP1", np.array(op1_read2, dtype=np.object))

# op1_prox = select(
#     op1,
#     from_start_min=exp2 + len(blocks["OP2"]),
#     from_start_max=exp2 + len(blocks["OP2"]) + 100,
#     min_score=ms_op1,
# )
# m_prox, ed_prox = align_stats(S, "OP1", op1_prox)

# op1_internal = select(
#     op1,
#     from_start_min=exp2 + len(blocks["OP2"]) + 100,
#     from_end_min=43 + len(blocks["N70X"]),
#     min_score=ms_op1,
# )
# m_int, ed_int = align_stats(S, "OP1", op1_internal)


### occurrence statistics for OP1, OP2, OP1_RC and OP2_RC
fig, ((ax, ax3), (ax2, ax4)) = plt.subplots(2, 2)
ax.set_yscale("log")
ax2.set_yscale("log")
ax3.set_yscale("log")
ax4.set_yscale("log")

L_max = max(op1[4].max(), op2[4].max(), op3[4].max())
x = np.arange(L_max)
ax.axvline(exp1, label="expect OP1", lw=0.1, color="blue", ls="dashed")
ax.axvline(exp2, label="expect OP2", lw=0.1, color="orange", ls="dashed")
ax.step(
    x,
    np.bincount(op1[1], minlength=L_max),
    where="mid",
    label="OP1",
    alpha=0.5,
    color="blue",
)
ax.step(
    x,
    np.bincount(op2[1], minlength=L_max),
    where="mid",
    label="OP2",
    alpha=0.5,
    color="orange",
)
ax.step(
    x,
    np.bincount(op3[1], minlength=L_max),
    where="mid",
    label="OP3",
    alpha=0.5,
    color="red",
)
ax.set_ylabel("match count")
ax.set_xlabel("position from start of cDNA")
ax.legend()

L_max = max(op1rc[4].max(), op2rc[4].max())
x = np.arange(L_max)
ax2.axvline(exp1, label="expect OP1", lw=0.1, color="blue", ls="dashed")
ax2.axvline(exp2, label="expect OP2", lw=0.1, color="orange", ls="dashed")
ax2.step(
    x,
    np.bincount(op1rc[1], minlength=L_max),
    where="mid",
    label="OP1_RC",
    alpha=0.5,
    color="blue",
)
ax2.step(
    x,
    np.bincount(op2rc[1], minlength=L_max),
    where="mid",
    label="OP2_RC",
    alpha=0.5,
    color="orange",
)
ax2.step(
    x,
    np.bincount(op3rc[1], minlength=L_max),
    where="mid",
    label="OP3_RC",
    alpha=0.5,
    color="red",
)
ax2.set_ylabel("match count")
ax2.set_xlabel("position from start of cDNA")
ax2.legend()


L_max = max(op1[4].max(), op2[4].max(), op3[4].max())
x = np.arange(L_max)
ax3.step(
    x,
    np.bincount(op1[4] - op1[2], minlength=L_max),
    where="mid",
    label="OP1",
    alpha=0.5,
    color="blue",
)
ax3.step(
    x,
    np.bincount(op2[4] - op2[2], minlength=L_max),
    where="mid",
    label="OP2",
    alpha=0.5,
    color="orange",
)
ax3.step(
    x,
    np.bincount(op3[4] - op3[2], minlength=L_max),
    where="mid",
    label="OP3",
    alpha=0.5,
    color="red",
)

ax3.set_ylabel("match count")
ax3.set_xlabel("position from END of cDNA")
ax3.legend()

L_max = max(op1rc[4].max(), op2rc[4].max(), op3rc[4].max())
x = np.arange(L_max)
ax4.step(
    x,
    np.bincount(op1rc[4] - op1rc[2], minlength=L_max),
    where="mid",
    label="OP1_RC",
    alpha=0.5,
    color="blue",
)
ax4.step(
    x,
    np.bincount(op2rc[4] - op2rc[2], minlength=L_max),
    where="mid",
    label="OP2_RC",
    alpha=0.5,
    color="orange",
)
ax4.step(
    x,
    np.bincount(op3rc[4] - op3rc[2], minlength=L_max),
    where="mid",
    label="OP3_RC",
    alpha=0.5,
    color="red",
)

ax4.set_ylabel("match count")
ax4.set_xlabel("position from END of cDNA")
ax4.legend()

# ax2.hist(lengths, bins=bins, label="all", histtype='step', density=True)
# ax2.hist(lengths[m1 & m2], bins=bins, label="with OP1 & OP2", histtype='step', density=True)
# for func, name in analyses:
#     _, m1 = detect_matches(filter_func=func)
#     # ax2.hist(lengths[m1], bins=bins, label=f"with OP1 {name}", histtype='step', density=True)
#     ax2.step(sorted(lengths[m1]), np.linspace(0,1, m1.sum()), label=f"with OP1 {name} (n={m1.sum()})", where='mid')

# ax2.set_ylabel("cumulative rel. frequency")
# ax2.set_xlabel("PB read length")
# ax2.legend()
fig.tight_layout()
fig.savefig(f"aln_pos_match_{sample}.pdf")
plt.close("all")


rcOP1 = analyze_recurrence(S, "OP1")
rcOP1_RC = analyze_recurrence(S, "OP1_RC")
rcOP2_RC = analyze_recurrence(S, "OP2_RC")

fig, ((ax, ax2), (ax3, ax4)) = plt.subplots(2, 2)
ax.hist(rcOP1_RC[0], bins=100, histtype="step", label="OP1_RC", cumulative=True)
ax.hist(rcOP2_RC[0], bins=100, histtype="step", label="OP2_RC", cumulative=True)
ax.set_xlabel("distance from bead oligo start")
ax.set_ylabel("cumulative frequency")
ax.legend()

ax2.hist(rcOP1_RC[1], bins=100, histtype="step", label="OP1_RC", cumulative=True)
ax2.hist(rcOP2_RC[1], bins=100, histtype="step", label="OP2_RC", cumulative=True)
ax2.set_xlabel("distance from bead oligo end")
ax2.set_ylabel("cumulative frequency")
ax2.legend()

ax3.hist(rcOP1_RC[2], bins=100, histtype="step", label="OP1_RC", cumulative=True)
ax3.hist(rcOP2_RC[2], bins=100, histtype="step", label="OP2_RC", cumulative=True)
ax3.set_xlabel("oligo match score")
ax3.set_ylabel("cumulative frequency")
ax3.legend()

ax4.hist(rcOP1_RC[3], bins=100, histtype="step", label="OP1_RC", cumulative=True)
ax4.hist(rcOP2_RC[3], bins=100, histtype="step", label="OP2_RC", cumulative=True)
ax4.set_xlabel("length of PB-read with a match")
ax4.set_ylabel("cumulative frequency")
ax4.legend()

fig.suptitle(sample)
fig.tight_layout()
fig.savefig(f"recurrent_OPx_{sample}.pdf")
plt.close("all")


from collections import Counter

for sig, count in Counter(S.signatures.values()).most_common(10):
    print(f"{count}\t{sig}")

sys.exit(0)


def select_start(seq, x):
    return x < 110


def select_prox(seq, x):
    return (x > 110) and (x < 250)


def select_ds(seq, x):
    return x >= 250


def select_read2(seq, x):
    return x > len(seq) - (43 + len(blocks["N70X"]))


analyses = [
    (select_start, "start"),
    (select_prox, "prox"),
    (select_ds, "ds"),
    (select_read2, "read2"),
]


def detect_matches(pattern=blocks["OP1"][2:12], filter_func=lambda seq, x: True):
    pos = np.zeros(L_max)
    mask = np.zeros(len(oriented), dtype=np.bool)
    for i, seq in enumerate(oriented):
        for m in re.finditer(pattern, seq):
            x = m.start()
            if not filter_func(seq, x):
                continue
            pos[x] += 1
            mask[i] = True
    #
    return pos, mask


def detect_matches_reverse(
    pattern=blocks["OP1_RC"][-12:-2], filter_func=lambda seq, x: True
):
    pos = np.zeros(L_max)
    mask = np.zeros(len(oriented), dtype=np.bool)
    for i, seq in enumerate(oriented):
        seq = rev_comp(seq)
        for m in re.finditer(pattern, seq):
            x = m.start()
            if not filter_func(seq, x):
                continue
            pos[x] += 1
            mask[i] = True
    #
    return pos, mask


exp1 = 1 + len(blocks["P5"] + blocks["SMART"]) + 9.5 + 2
exp2 = 1 + len(blocks["P5"] + blocks["SMART"] + blocks["OP1"]) + 9.5 + 8

op1_pos, m1 = detect_matches()
op2_pos, m2 = detect_matches(pattern=blocks["OP2"][:10])

fig, (ax, ax2) = plt.subplots(2, sharex=True)
ax.set_yscale("log")
x = np.arange(L_max)
ax.axvline(exp1, label="expect OP1", lw=0.1, color="blue", ls="dashed")
ax.axvline(exp2, label="expect OP2", lw=0.1, color="orange", ls="dashed")
ax.step(x, op1_pos, where="mid", label="OP1", alpha=0.5, color="blue")
ax.step(x, op2_pos, where="mid", label="OP2", alpha=0.5, color="orange")
ax.set_ylabel("match count")
ax.set_xlabel("position in PB read")
ax.legend()

# subset_lengths1 = []
# subset_lengths2 = []
# for func, name in analyses:
#     _, m1 = detect_matches(filter_func=func)
#     _, m2 = detect_matches(filter_func=func, pattern=blocks['OP2'][:10])
#     subset_lengths1.append(lengths[m1])
#     subset_lengths2.append(lengths[m2])

bins = 100
# ax2.hist(lengths, bins=bins, label="all", histtype='step', density=True)
# ax2.hist(lengths[m1 & m2], bins=bins, label="with OP1 & OP2", histtype='step', density=True)
for func, name in analyses:
    _, m1 = detect_matches(filter_func=func)
    # ax2.hist(lengths[m1], bins=bins, label=f"with OP1 {name}", histtype='step', density=True)
    ax2.step(
        sorted(lengths[m1]),
        np.linspace(0, 1, m1.sum()),
        label=f"with OP1 {name} (n={m1.sum()})",
        where="mid",
    )

ax2.set_ylabel("cumulative rel. frequency")
ax2.set_xlabel("PB read length")
ax2.legend()
fig.tight_layout()
fig.savefig(f"pos_{sample}.pdf")
plt.close("all")


op1_distal, m1 = detect_matches_reverse()
op2_distal, m2 = detect_matches_reverse(pattern=blocks["OP2_RC"][-10:])

ill_start = len(blocks["N70X"])
ill_end = ill_start + 43

fig, ax = plt.subplots()
ax.set_yscale("log")
x = np.arange(L_max)
ax.axvline(ill_start, label="ill r2 start", lw=0.1, color="k", ls="dashed")
ax.axvline(ill_end, label="ill r2 end (43nt)", lw=0.1, color="k", ls="dashed")
ax.bar(x, op1_distal, label="OP1_RC", alpha=0.5, color="blue")
ax.bar(x, op2_distal, label="OP2_RC", alpha=0.5, color="orange")
ax.set_xlabel("position from end of PB read")
ax.set_ylabel("total count")
ax.legend()
fig.tight_layout()
fig.savefig(f"pos_distal_{sample}.pdf")
plt.close("all")


def nuc_count(mat, seq, pos, flank=100):
    idx = {
        "A": 0,
        "C": 1,
        "G": 2,
        "T": 3,
    }
    for x in range(max(0, pos - flank), min(pos + flank + 1, len(seq))):
        j = x - pos + flank
        nt = seq[x]
        if nt in idx:
            mat[idx[nt], j] += 1


def scan_sequence_matches(
    seqs, pos_filter=lambda seq, x: x > 250, pattern=blocks["OP1"][2:12]
):
    mat = np.ones((4, 201)) * 5  # 5 pseudo counts
    lengths = []
    for i, seq in enumerate(seqs):
        for m in re.finditer(pattern, seq):
            x = m.start()
            if pos_filter(seq, x):
                nuc_count(mat, seq, x)
                lengths.append(len(seq))
    #
    counts = mat.sum(axis=0)
    freq = mat / counts[np.newaxis, :]
    ent = 2 + np.where(freq > 0, freq * np.log2(freq), 0).sum(axis=0)
    return freq * ent[np.newaxis, :], np.array(lengths)


from RBPamp.affinitylogo import _draw_logo


def plot_seqlogo(ax, pfm, info=False, charwidth=1.0, **kwargs):
    if info:
        info_content = 2 - pfm.apply(lambda p: (-p * np.log2(p)).sum(), axis=1)
        matrix = pfm.mul(info_content, axis=0)
    else:
        matrix = pfm
    #
    _draw_logo(ax, matrix, charwidth, **kwargs)


for func, name in analyses:
    flank = 100
    mat1, L1 = scan_sequence_matches(
        np.array(oriented)[m1], pattern=blocks["OP1"][2:12], pos_filter=func
    )
    mat2, L2 = scan_sequence_matches(
        np.array(oriented)[m2], pattern=blocks["OP2"][:10], pos_filter=func
    )
    #
    fig, (ax1, ax2) = plt.subplots(2, figsize=(16, 4), sharey=True, sharex=True)
    fig.suptitle(f"{sample} {name} (n1={len(L1)} n2={len(L2)})")
    plot_seqlogo(ax1, pd.DataFrame(mat1.T, columns=list("ACGT")))
    plot_seqlogo(ax2, pd.DataFrame(mat2.T, columns=list("ACGT")))
    #
    x = np.arange(-flank, flank + 1, 20)
    ax2.set_xticks(x + flank)
    ax2.set_xticklabels(x)
    ax2.set_ylim(0, 2)
    ax2.set_xlim(0, 201)
    ax2.set_ylabel("bits")
    ax1.set_ylabel("bits")
    ax1.set_xlabel(f"pos rel to {name} OP1 ({blocks['OP1']})")
    ax2.set_xlabel(f"pos rel to {name} OP2 ({blocks['OP2']})")
    # ax2.set_ylim(0, 2)
    fig.tight_layout()
    fig.savefig(f"context_seqlogo_{sample}_{name}.pdf")
    plt.close()
