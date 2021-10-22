import pandas as pd
import numpy as np
import logging
from spacemake.util import rev_comp
from spacemake.util import read_fq
from collections import defaultdict


def sig2str(sig, max_repeat_only=True):
    """
    basically merge oligo block names with ',' but contract homopolymeric repeats by
    adding a '+' suffix instead of repeating the thing over and over
    """
    if not len(sig):
        return ""

    parts = [sig[0]]
    poly_runs = set()
    repeat_counts = defaultdict(int)
    for s in sig[1:]:
        if s != parts[-1]:
            parts.append(s)
        else:
            poly_runs.add(len(parts) - 1)
            repeat_counts[s] += 1

    for i in poly_runs:
        parts[i] += "+"

    if poly_runs and max_repeat_only:
        # This is a more drastic reduction of the original signature: we drop
        # all other building blocks and disregard orientation. If any block X
        # is tandem repeated we call the signature of this just 'X+'.
        max_rep, nrep = sorted(repeat_counts.items(), key=lambda x: -x[1])[0]
        return f"?{max_rep.replace('_RC', '')}+"
    else:
        return ",".join(parts)


class AnnotatedSequences:
    def __init__(
        self,
        fastq_path,
        ann_path,
        sample_name,
        blocks,
        min_score=0,
        orient_by="bead_start",
    ):
        self.sample_name = sample_name
        self.logger = logging.getLogger("AnnotatedSequences")
        self.orient_by = orient_by
        self.orient_by_RC = orient_by + "_RC"
        self.raw_sequences = self.load_raw_sequences(fastq_path)
        self.oligo_blocks = blocks
        self.min_oligo_scores = {}
        for name, seq in blocks.items():
            self.min_oligo_scores[name] = (2 * len(seq)) * min_score

        self.ann_db = self.load_annotation(ann_path)
        # self.ann_db = self.cleanup_overlaps(self.ann_db)

        self.signatures = self.extract_signatures_and_orient()

    def load_raw_sequences(self, fastq_path):
        self.logger.info(
            f"loading raw sequences for {self.sample_name} from {fastq_path}"
        )
        raw_sequences = {}
        for qname, seq, qual in read_fq(fastq_path):
            raw_sequences[qname] = seq

        return raw_sequences

    def load_annotation(self, ann_path):
        self.logger.info(
            f"loading oligo annotation for {self.sample_name} from {ann_path}"
        )
        df = pd.read_csv(ann_path, sep="\t")
        df["min_score"] = df["oligo"].apply(lambda x: self.min_oligo_scores.get(x, 22))
        df = df.query(f"score > min_score")
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
            lnames = list(names)
            fw_start = lnames.count(self.orient_by)
            rc_start = lnames.count(self.orient_by_RC)

            if (fw_start < rc_start) or (
                (fw_start == rc_start) and names[0].endswith("_RC")
            ):
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

    def filter_signatures(self, qsig, substring=False):
        # print(f"qsig={qsig}")
        nmatch = 0
        qstr = ",".join(qsig)
        lq = len(qsig)
        for qname, sig in self.signatures.items():
            sstr = sig2str(sig)
            if substring:
                # demand match somewhere in signature (substring)
                i = sstr.find(qstr)
                if i > -1:
                    ofs = sstr[:i].count(",")
                    # print(qstr, sstr)
                    nmatch += 1
                    yield qname, ofs, ofs + lq
            else:
                # demand *exact* identity for a match
                # print(f"sig='{sig}' qsig='{qsig}'")
                if sstr == qstr:
                    nmatch += 1
                    yield qname, 0, lq

        self.logger.info(
            f"filter_signatures({qstr}, substring={substring}) -> {nmatch} ({100*nmatch/len(self.signatures):.2f}%) hits"
        )

    def count_signatures(self):
        sig_counts = defaultdict(int)
        concat = 0
        reprimed = 0
        for sig in self.signatures.values():
            sstr = sig2str(sig)
            sig_counts[sstr] += 1
            if "+" in sstr:
                concat += 1
            if len(sig) > 1:
                first = sig[0]
                last = sig[-1]
                if (first + "_RC").replace("_RC_RC", "") == last:
                    reprimed += 1

        return sig_counts, concat, reprimed

    def count_concatenations(self):
        concat = defaultdict(int)
        n_occ = 0
        for sig in self.signatures.values():
            sstr = sig2str(sig)
            nc = sstr.count("+")
            if nc:
                n_occ += nc
                for part in sstr.split(","):
                    if part.endswith("+"):
                        concat[part] += 1

        return concat, n_occ

    def count_repriming(self):
        reprimed = defaultdict(int)
        for sig in self.signatures.values():
            first = sig[0]
            last = sig[-1]
            first_rc = (first + "_RC").replace("_RC_RC", "")
            if first_rc == last:
                reprimed[first.replace("_RC", "")] += 1

        return reprimed

    def extract_cDNA(self, qname, after_oligo="bead_start", distal=150):
        """
        Ensure the read has a match to <after_oligo>. Find the last
        oligo-match following <after_oligo> and still within <distal> nt
        from start. The end of that match is the start of cDNA.
        On the other end, select the innermost oligo within <distal> nt from
        the end. The start of that match is the end of cDNA. Lastly, rev_comp
        the cDNA as Illumina would sequence from the 3'end.
        """
        read_seq = self.raw_sequences[qname]
        names, starts, ends, _ = self.ann_db[qname]
        if after_oligo in names:
            i = names.index(after_oligo)
        else:
            i = 0

        cDNA_start = 0
        L = len(read_seq)
        cDNA_end = L
        cDNA_end_oli = ""
        cDNA_start_oli = ""

        for name, start, end in zip(names[i:], starts[i:], ends[i:]):
            if start <= distal:
                cDNA_start = end
                cDNA_start_oli = name
            if end >= L - distal:
                cDNA_end = min(cDNA_end, start)
                cDNA_end_oli = name

        return (
            (cDNA_start, cDNA_end),
            (cDNA_start_oli, cDNA_end_oli),
            rev_comp(read_seq[cDNA_start:cDNA_end]),
        )

    def extract_between(
        self,
        qname,
        after="bead_start",
        before="polyT",
        srange=(0, np.inf),
        min_L=0,
        max_L=np.inf,
    ):
        """
        Useful for extracting barcodes and potentially also cDNA (TODO: unify).
        Looks for matches to <after> and <before> within the <range> interval of the
        oriented long read. Then excises the intervening sequence,
        if min_L <= length <= max_L .
        """

        read_seq = self.raw_sequences[qname]
        names, starts, ends, _ = self.ann_db[qname]

        extracted = ""
        s = None
        e = None
        if len(names) > 1:
            # iterate over consecutive pairs of detected oligo
            # matches and their start & end coordinates
            for (n1, s1, e1), (n2, s2, e2) in zip(
                zip(names[:-1], starts[:-1], ends[:-1]),
                zip(names[1:], starts[1:], ends[1:]),
            ):

                if s1 > srange[1] or e2 < srange[0]:
                    # outside the area we are scanning
                    continue

                if (n1 == after) and (n2 == before):
                    L = s2 - e1
                    if min_L <= L <= max_L:
                        s = e1
                        e = s2
                        extracted = read_seq[e1:s2]

                        break

        return extracted, (s, e)

    def query_dimensions(self, qsig, substring=False):
        qnames = []
        L = []
        starts = []
        ends = []
        scores = []

        for qname, i, j in self.filter_signatures(qsig, substring=substring):
            qnames.append(qname)
            L.append(len(self.raw_sequences[qname]))
            names, s, e, sc = self.ann_db[qname]
            starts.append(s[i:j])
            ends.append(e[i:j])
            scores.append(sc[i:j])

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
        names, starts, ends, scores = self.ann_db[qname]
        # print(names)
        # print(starts)
        # print(ends)
        # print(scores)
        for oligo, start, end in zip(names, starts, ends):
            buf[start : end + 1] = render_label(oligo, start, end + 1)
        #
        return f"# qname={qname} oligo_matches={names} match_scores={scores}\n{seq}\n{''.join(buf)}"

    def completeness(self, sig_intact, polyT="polyT"):
        partial = []
        complete = ",".join(sig_intact)
        for i in range(len(sig_intact)):
            partial.append(",".join(sig_intact[: i + 1]))
            # print("partial append:", partial[-1])

        partial = partial[::-1]

        partial_counts = defaultdict(int)
        pT_counts = defaultdict(int)
        prefixes = defaultdict(int)
        suffixes = defaultdict(int)

        n = 0
        pT = f",{polyT}"
        # print(f"name of polyT oligo is '{polyT}' so scanning fot '{pT}'")
        for sig in self.signatures.values():
            n += 1
            sig_str = ",".join(sig)
            for i, p in enumerate(partial):
                # because we go in reverse order,
                # the first match is maximal. After that
                # we stop
                if p in sig_str:  #  and not (p + "_RC" in sig_str)
                    partial_counts[p] += 1
                    if pT in p:
                        # pT is part of the complete signature
                        pT_counts[p.rsplit(pT)[0]] += 1

                    assert partial_counts[p] <= n

                    if i == 0:
                        # if this is the maximal match (i==0)
                        # extract prefix and suffix
                        left_context = sig_str.split(p, maxsplit=1)
                        right_context = sig_str.rsplit(p, maxsplit=1)
                        prefixes[
                            left_context[0][:-1]
                        ] += 1  # [:-1] omits the trailing ','
                        if len(right_context) > 1:
                            if right_context[1].startswith(pT):
                                pT_counts[p] += 1
                            suffixes[
                                right_context[1][1:]
                            ] += 1  # [1:] omits the leading ','
                    break

        return partial_counts, prefixes, suffixes, pT_counts


def align_stats(ann, oseq, qmatches, pad=1):
    import spacemake.longread.cache as cache
    from Bio import pairwise2

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
        # print("\n".join(aln_str))
        for r, m, q in zip(aln_str[0], aln_str[1], aln_str[2]):
            # print(r, m, q, i, len(oseq))
            if q == " ":
                # outside of oligo match range
                continue
            if m == " " and i == 0:
                # insertion at position zero is bananas for a local alignment
                # why do we even get this?
                continue
            if m == " " and i >= len(oseq):
                # same as above, we've already aligned all of the oligo sequence
                # what point is there in this "insertion"?
                continue
            if m == "|":
                # print(f"match {i}", matches[i])
                matches[i] += 1
            else:
                edits[i][q + r] += 1
                # print(f"edit {i} {q + r}")
            if q != "-":
                i += 1

    return matches / n, edits


def subsample(qmatches, n=100):
    data = np.array(qmatches, dtype=np.object)
    n_cols, N = data.shape
    I = np.random.choice(N, size=n, replace=False)
    return data[:, I]


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


#     stages[0] += 1
#     # print(sig)
#     l = list(sig)
#     if "SMART" in sig:
#         l = l[sig.index("SMART") + 1 :]
#         stages[1] += 1
#         try:
#             block = l.pop(0)
#             if block == "pT":
#                 pT[0] += 1
#             elif block == "OP1":
#                 stages[2] += 1
#                 block = l.pop(0)
#                 if block == "pT":
#                     pT[1] += 1
#                 elif block == "OP2":
#                     stages[3] += 1
#                     block = l.pop(0)
#                     if block == "pT":
#                         pT[2] += 1
#                     elif block == "OP3":
#                         stages[4] += 1
#                         block = l.pop(0)
#                         if block == "pT":
#                             pT[3] += 1
#         except IndexError:
#             # this just happens when l is empty
#             pass

# sig_intact = ("P5", "SMART", "OP1", "OP2", "OP3", "pT", "N70X")
# S = AnnotatedSequences(
#     sys.argv[1],
#     blocks,
#     path=base_path,
# )
# qintact, qL, qstarts, qends, qsores = S.query_dimensions(sig_intact)


# # estimate extension efficiencies
# # n_reads = 0
# stages = np.zeros(5, dtype=float)
# pT = np.zeros(4)
# for qname, sig in S.signatures.items():
#     stages[0] += 1
#     # print(sig)
#     l = list(sig)
#     if "SMART" in sig:
#         l = l[sig.index("SMART") + 1 :]
#         stages[1] += 1
#         try:
#             block = l.pop(0)
#             if block == "pT":
#                 pT[0] += 1
#             elif block == "OP1":
#                 stages[2] += 1
#                 block = l.pop(0)
#                 if block == "pT":
#                     pT[1] += 1
#                 elif block == "OP2":
#                     stages[3] += 1
#                     block = l.pop(0)
#                     if block == "pT":
#                         pT[2] += 1
#                     elif block == "OP3":
#                         stages[4] += 1
#                         block = l.pop(0)
#                         if block == "pT":
#                             pT[3] += 1
#         except IndexError:
#             # this just happens when l is empty
#             pass
