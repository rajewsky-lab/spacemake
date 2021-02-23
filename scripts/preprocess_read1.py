#!/usr/bin/env python3
import argparse
import logging
import timeit as ti
import pandas as pd
import numpy as np
import multiprocessing as mp
from collections import defaultdict
from more_itertools import grouper
from Bio import pairwise2
from Bio import SeqIO
from collections import namedtuple
atype = namedtuple("atype", "seqA seqB score start end")
NO_CALL="NNNNNNNN"


def read_fq(fname):
    for name, seq, _, qual in grouper(open(fname), 4):
        yield name.rstrip(), seq.rstrip(), qual.rstrip()


def read_source():
    if args.read2 and ("r2" in args.cell or "r2" in args.UMI):
        for (id1, seq1, qual1), (id2, seq2, qual2) in zip(read_fq(args.read1), read_fq(args.read2)):
            assert id1.split()[0] == id2.split()[0]
            yield id1, seq1, seq2
    else:
        for id1, seq1, qual1 in read_fq(args.read1):
            yield id1, seq1, "READ2 IS NOT AVAILABLE"


def hamming(seqA, seqB, costs, match=2):
    score = 0
    for a, b, c in zip(seqA, seqB, costs):
        if a == b:
            score += match
        else:
            if a == 'A' and c == 0:
                # we are looking at the BC2 start position
                # This tweak gives almost ~1% point more matches 72.58% -> 73.43%
                score += 1 
            score -= c
    return score


class BarcodeMatcher:
    def __init__(self, fname, length_specific=True, place="left"):
        self.logger = logging.getLogger("BarcodeMatcher")
        self.length_specific = length_specific
        self.names, self.seqs = self.load_targets(fname)
        self.slen = np.array([len(s) for s in self.seqs])
        self.lmin = self.slen.min()
        self.lmax = self.slen.max()
        self.place = place
        self.costs = {}
        self.slen_masks = {}
        for l in range(self.lmin, self.lmax + 1):
            self.slen_masks[l] = (self.slen == l).nonzero()[0]
            cost = np.ones(l)
            if place == "left":
                cost[-8:] = 2
                cost[-4:] = 3  # the diagnostic 4-mer

            if place == "right":
                cost[0] = 0  # the error-prone first base
                cost[1:8] = 2
                cost[4:8] = 3  # the diagnostic 4-mer

            self.costs[l] = cost

        self.logger.info(f"initialized from {len(self.names)} sequences from lmin={self.lmin} to lmax={self.lmax}")

    def load_targets(self, fname):
        self.logger.info(f"loading target sequences from '{fname}'")
        names = []
        seqs = []
        for rec in SeqIO.parse(fname, "fasta"):
            names.append(rec.id.split()[0])
            seqs.append(str(rec.seq).upper())

        return np.array(names), np.array(seqs)

    def align(self, query, debug=False):
        results = []
        scores = []
        # select set of barcode sequences to align with
        if len(query) < self.lmin:
            return [NO_CALL, ], [NO_CALL, ], [-1, ]

        if self.length_specific:
            lq = min(len(query), self.lmax)
            lmask = self.slen_masks[lq]
            seqs_sel = self.seqs[lmask]
            names_sel = self.names[lmask]
        else:
            seqs_sel = self.seqs
            names_sel = self.names

        costs = self.costs[lq]
        # perform the alignments and keep results
        for seq in seqs_sel:
            res = hamming(query, seq, costs)
            # print(query, seq, res)
            results.append(res)
            scores.append(res)

        # identify best alignments and return tied, best barcode matches
        scores = np.array(scores)
        if debug:
            I = scores.argsort()[::-1]
            print("Q  ", query)
            for i in I:
                seq = seqs_sel[i]
                res = results[i]
                print("*  ", seq, res)

        i = scores.argmax()
        ties = scores == scores[i]
        return names_sel[ties], seqs_sel[ties], scores[ties]


class TieBreaker:
    def __init__(self, fname, place="left"):
        self.logger = logging.getLogger("TieBreaker")
        self.matcher = BarcodeMatcher(fname, place=place)
        self.query_count = defaultdict(float)
        self.bc_count = defaultdict(float)
        self.cache = {}
        self.n_hit = 0
        self.n_align = 0

    def align(self, query, debug=False, w=1):
        self.query_count[query] += w
        if not query in self.cache:
            self.n_align += w

            names, seqs, scores = self.matcher.align(query, debug)
            if debug:
                for n, s, S in zip(names, seqs, scores):
                    print(f"{n}\t{s}\t{S}")

            if (len(names) == 1) or (len(set(names)) == 1):
                # unambiguous best hit
                result = names[0], seqs[0], scores[0]
            else:
                # potentially do more involved resolving?
                result = (NO_CALL, NO_CALL, scores[0])

            self.cache[query] = result

        else:
            self.n_hit += w
            result = self.cache[query]

        self.bc_count[result[0]] += w
        self.bc_count['total'] += w
        return result

    def align_choices(self, queries, debug=False):
        w = 1.0 / len(queries)
        results = [self.align(q, debug=debug, w=w) for q in queries]
        score = np.array([(r[0] != NO_CALL) * r[2] for r in results])
        i = score.argmax()
        ties = (score == score.argmax())
        if ties.sum() > 1:
            return (NO_CALL, NO_CALL, -1), queries[0]
        else:
            return results[i], queries[i]

    def store_cache(self, fname, mincount=2):
        with open(fname, 'w') as f:
            for query in sorted(self.cache.keys(), key=lambda x: x[0] if len(x) else x):
                if self.query_count[query] >= mincount:
                    name, seq, score = self.cache[query]
                    f.write(f"{query}\t{seq}\t{name}\t{score}\t{self.query_count[query]}\n")

    def load_cache(self, fname):
        self.logger.info(f"pre-populating alignment cache from '{fname}'")
        if fname:
            try:
                df = pd.read_csv(fname, sep='\t', index_col=None, names=["query", "seq", "name", "score", "count"])
            except OSError as err:
                self.logger.error(f"error while loading: {err}")
            else:
                for row in df.itertuples():
                    self.cache[row.query] = (row.name, row.seq, row.score)

        self.logger.info(f"loaded {len(self.cache)} queries.")


def output(rname, cell=None, UMI=None, bc1=None, bc2=None, qual="E", r1="", r2=""):
    if bc1 is None:
        bc1 = args.na

    if bc2 is None:
        bc2 = args.na

    # slightly concerned about security here... perhaps replace
    # all () and ; with sth else?
    cell = eval(args.cell)
    UMI = eval(args.UMI)
    seq = args.template.format(**locals())
    # print(bc1, bc2)
    # print(args.template)
    # print(seq)
    print(f"{rname}\n{seq}\n+\n{qual*len(seq)}")


def report_stats(N, prefix=""):
    for k, v in sorted(N.items()):
        logging.info(f"{prefix}{k}\t{v}\t{100.0 * v/N['total']:.2f}")


def match_BC1(bc1_matcher, seq, qstart, tstart, N):
    if tstart == 0:  # the start of opseq primer is intact
        bc1 = seq[:qstart]
        BC1, ref1, score1 = bc1_matcher.align(bc1)
    else:
        # we have a possible deletion, or mutation. Check both options
        bc1_choices = [
            seq[:qstart],  # deletion
            seq[:qstart - tstart]  # mutation
        ]
        (BC1, ref1, score1), bc1 = bc1_matcher.align_choices(bc1_choices)

    N[f'BC1_score_{score1}'] += 1
    if BC1 == NO_CALL:
        N['BC1_ambig'] += 1
    elif score1 < 2*len(BC1) * args.threshold:
        N['BC1_low_score'] += 1
        BC1 = NO_CALL
    else:
        N['BC1_assigned'] += 1

    return bc1, BC1, ref1, score1


def match_BC2(bc2_matcher, seq, qend, tend, N, lopseq=22):
    bc2_choices = [
        seq[qend:],
        seq[qend - 1:],
    ]
    # we have a possible deletion, or mutation. Check both options
    if tend:
        mut = seq[qend + tend:]
        if len(mut) >= 8:
            bc2_choices.append(mut)

    # print(f"bc2 choices {bc2_choices}")
    (BC2, ref2, score2), bc2 = bc2_matcher.align_choices(bc2_choices)

    N[f'BC2_score_{score2}'] += 1
    if BC2 == NO_CALL:
        N['BC2_ambig'] += 1
    elif score2 < 2*len(BC2) * args.threshold:
        N['BC2_low_score'] += 1
        BC2 = NO_CALL
    else:
        N['BC2_assigned'] += 1

    return bc2, BC2, ref2, score2


def opseq_local_align(seq, opseq="GAATCACGATACGTACACCAGT", min_opseq_score=22):
    res = pairwise2.align.localmd(opseq, seq,
                                  2, -1.5, -3, -1, -3, -1,
                                  one_alignment_only=True,
                                  score_only=False)[0]

    res = atype(*res) # convert unpickle-able Alignment class to named tuple
    # print(pairwise2.format_alignment(*res, full_sequences=True))
    # print(res)
    qstart = res.start
    qend = res.end
    tstart = 0
    tend = 0
    L = len(res.seqB)
    # print(f"flags {qstart < 8} {qend > (len(res.seqB) - 8)} {res.score < min_opseq_score}")
    if (qstart < 8) or (qend > (L - 8)) or (res.score < min_opseq_score):
        res = None
    else:
        # check for start/end gaps
        # backtrack from qstart in seqA until we hit '-'
        while res.seqA[qstart - tstart - 1] != '-':
            if qstart - tstart - 1 <= 0:
                tstart = -1
                res = None
                break
            tstart += 1

        if res is not None:
            # and scan forward from qend until we hit '-'
            while res.seqA[qend + tend] != '-':
                if qend + tend >= len(res.seqA) - 8:
                    tend = -1
                    res = None
                    break
                tend += 1

    return res, tstart, tend


def map_opseq(fq):
    fqid, seq1, seq2 = fq
    res, tstart, tend = opseq_local_align(seq1.rstrip(), opseq=args.opseq, min_opseq_score=args.min_opseq_score)
    return fqid, seq1, seq2, res, tstart, tend


def main_combinatorial(args):
    bc1_matcher = TieBreaker(args.bc1_ref, place="left")
    bc2_matcher = TieBreaker(args.bc2_ref, place="right")
    bc1_matcher.load_cache(args.bc1_cache)
    bc2_matcher.load_cache(args.bc2_cache)
    lopseq = len(args.opseq)

    N = defaultdict(int)
    with mp.Pool(args.parallel) as pool:
        # iterate query sequences
        for n, (fqid, r1, r2, res, tstart, tend) in enumerate(pool.imap(map_opseq, read_source(), chunksize=1000)):
            N['total'] += 1
            if res is None:
                N['opseq_broken'] += 1
                output(fqid, r1=r1, r2=r2)
                continue

            bc1, BC1, ref1, score1 = match_BC1(bc1_matcher, res.seqB, res.start, tstart, N)
            bc2, BC2, ref2, score2 = match_BC2(bc2_matcher, res.seqB, res.end, tend, N,
                                               lopseq=lopseq)

            # slo = sQSeq.lower()
            # sout = slo[:qstart] + sQSeq[qstart:qend] + slo[qend:]
            # print(sout, qstart, qend, tstart, tend, bc1, bc2)
            # print(f"bc1: {bc1} -> {ref1} -> {BC1} score={50.0*score1/len(bc1):.1f} %")
            # print(f"bc2: {bc2} -> {ref2} -> {BC2} score={50.0*score2/len(bc2):.1f} %")

            if BC1 != NO_CALL and BC2 != NO_CALL:
                N["called"] += 1

            output(fqid, bc1=BC1, bc2=BC2, r1=r1, r2=r2)
            # print(sQId, sQSeq.rstrip(), "score", score, "target_begin:", tstart, "target_end", tend, "query_begin", qstart, "query_end", qend, sQSeq[qstart:qend], "bc1", bc1, "BC1", BC1, "score=", score1, "bc2", bc2, "BC2", BC2, "score=", score2)

            if n and n % 100000 == 0:
                N['BC1_cache_hit'] = bc1_matcher.n_hit
                N['BC2_cache_hit'] = bc2_matcher.n_hit
                report_stats(N)

    N['BC1_cache_hit'] = bc1_matcher.n_hit
    N['BC2_cache_hit'] = bc2_matcher.n_hit
    report_stats(N)

    if args.update_cache:
        bc1_matcher.store_cache(args.bc1_cache)
        bc2_matcher.store_cache(args.bc2_cache)

    # report_stats(bc1_matcher.bc_count, prefix="BC1\t")
    # report_stats(bc2_matcher.bc_count, prefix="BC2\t")
    return N, bc1_matcher, bc2_matcher


def main_dropseq(args):
    N = defaultdict(int)
    # with mp.Pool(args.parallel) as pool:
    # iterate query sequences
    for n, (fqid, r1, r2) in enumerate(read_source()):
        N['total'] += 1
        output(fqid, r1=r1, r2=r2)

        if n and n % 100000 == 0:
            report_stats(N)

    report_stats(N)
    return N


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert combinatorial barcode read1 sequences to Dropseq pipeline compatible read1 FASTQ')
    parser.add_argument("--read1", default="/dev/stdin", help="source from where to get read1 (FASTQ format)")
    parser.add_argument("--read2", default="", help="source from where to get read2 (FASTQ format)")
    parser.add_argument("--parallel", default=1, type=int, help="how many processes to spawn")
    parser.add_argument("--opseq", default="GAATCACGATACGTACACCAGT", help="opseq primer site sequence")
    parser.add_argument("--min-opseq-score", default=22, type=float, help="minimal score for opseq alignment (default 22 [half of max])")
    parser.add_argument("--cell", default="bc1[-6:][::-1]+bc2[-6:][::-1]")  #"r1[0:12:-1]")
    parser.add_argument("--UMI", default="r1[12:16]")
    parser.add_argument("--template", default="{cell}{UMI}")
    parser.add_argument("--bc1-ref", default="", help="load BC1 reference sequences from this FASTA")
    parser.add_argument("--bc2-ref", default="", help="load BC2 reference sequences from this FASTA")
    parser.add_argument("--bc1-cache", default="", help="load cached BC1 alignments from here")
    parser.add_argument("--bc2-cache", default="", help="load cached BC2 alignments from here")
    parser.add_argument("--update-cache", default=False, action="store_true", help="update cache when run is complete (default=False)")
    parser.add_argument("--threshold", default=0.75, type=float, help="score threshold for calling a match (rel. to max for a given length)")
    parser.add_argument("--na", default="NNNNNNNN", help="code for ambiguous or unaligned barcode")
    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG)
    logging.info(f"starting read1 preprocessing with {sorted(vars(args).items())}")
    t1 = ti.default_timer()
    if args.bc1_ref and args.bc2_ref:
        main_combinatorial(args)
    else:
        main_dropseq(args)

    t2 = ti.default_timer()
    logging.info(f'CPU time: {t2 - t1:.3f} seconds\n')
