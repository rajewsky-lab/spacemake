#!/usr/bin/env python3
__version__ = "0.9"
__author__ = ["Marvin Jens"]
__license__ = "MIT"
__email__ = ["marvin.jens@mdc-berlin.de"]

import argparse
import logging
import time
import os
import sys
from numpy.core.fromnumeric import argsort
import pandas as pd
import numpy as np
import pysam
import multiprocessing as mp
from collections import defaultdict
from Bio import pairwise2
from Bio import SeqIO

from spacemake.parallel import (
    put_or_abort,
    queue_iter,
    join_with_empty_queues,
    chunkify,
    ExceptionLogging,
)
from spacemake.util import read_fq

NO_CALL = "NNNNNNNN"

# TODO:
# * add support for three segments
# * turn Output into a context manager
# * load params from YAML (code from opseq)
# * save params to YAML in run-folder to document
# * get rid of ugly global NO_CALL
# * refactor into multiple modules?


def read_source(args):
    if args.read2:
        for (id1, seq1, qual1), (id2, seq2, qual2) in zip(
            read_fq(args.read1), read_fq(args.read2)
        ):
            # assert id1.split()[0] == id2.split()[0]
            yield id1, seq1, id2, seq2, qual2
    else:
        for id1, seq1, qual1 in read_fq(args.read1):
            yield id1, seq1, id1, "READ2 IS NOT AVAILABLE", "READ2 IS NOT AVAILABLE"


def hamming(seqA, seqB, costs, match=2):
    score = 0
    for a, b, c in zip(seqA, seqB, costs):
        if a == b:
            score += match
        else:
            if a == "A" and c == 0:
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
        if len(self.names) == 0:
            self.logger.warning(f"no references loaded! Disabling matching for {fname}")
            self.align = self.align_na
        else:
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

            self.logger.debug(
                f"initialized from {len(self.names)} sequences from lmin={self.lmin} to lmax={self.lmax}"
            )

    def load_targets(self, fname):
        self.logger.debug(f"loading target sequences from '{fname}'")
        names = []
        seqs = []
        if fname:
            for rec in SeqIO.parse(fname, "fasta"):
                names.append(rec.id.split()[0])
                seqs.append(str(rec.seq).upper())

        return np.array(names), np.array(seqs)

    def align_na(self, query, debug=False):
        return (
            [NO_CALL],
            [NO_CALL],
            [-1],
        )

    def align(self, query, debug=False):
        results = []
        scores = []
        # select set of barcode sequences to align with
        if len(query) < self.lmin:
            return (
                [NO_CALL],
                [NO_CALL],
                [-1],
            )

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
        self.bc_count["total"] += w
        return result

    def align_choices(self, queries, debug=False):
        w = 1.0 / len(queries)
        results = [self.align(q, debug=debug, w=w) for q in queries]
        score = np.array([(r[0] != NO_CALL) * r[2] for r in results])
        i = score.argmax()
        ties = score == score.argmax()
        if ties.sum() > 1:
            return (NO_CALL, NO_CALL, -1), queries[0]
        else:
            return results[i], queries[i]

    def load_cache(self, fname):
        self.logger.debug(f"pre-populating alignment cache from '{fname}'")
        if fname:
            try:
                df = pd.read_csv(
                    fname,
                    sep="\t",
                    index_col=None,
                    names=["query", "seq", "name", "score", "count"],
                )
            except OSError as err:
                self.logger.warning(f"error while loading caches: {err}")
            else:
                for row in df.itertuples():
                    self.cache[row.query] = (row.name, row.seq, row.score)

        self.logger.debug(f"loaded {len(self.cache)} queries.")


def store_cache(fname, cache, query_count, mincount=2):
    with open(fname, "w") as f:
        print(cache.keys())
        for query in sorted(cache.keys()):
            if query_count[query] >= mincount:
                name, seq, score = cache[query]
                f.write(f"{query}\t{seq}\t{name}\t{score}\t{query_count[query]}\n")


def report_stats(N, prefix=""):
    for k, v in sorted(N.items()):
        logging.info(f"{prefix}{k}\t{v}\t{100.0 * v/N['total']:.2f}")


def match_BC1(bc1_matcher, seq, qstart, tstart, N, debug=False, threshold=0.5):
    if tstart == 0:  # the start of opseq primer is intact
        bc1 = seq[qstart - 8 : qstart]
        BC1, ref1, score1 = bc1_matcher.align(bc1, debug=debug)
    else:
        # we have a possible deletion, or mutation. Check both options
        bc1_choices = [
            seq[qstart - 8 : qstart],
            seq[qstart - tstart - 8 : qstart - tstart],
        ]  # deletion  # mutation
        (BC1, ref1, score1), bc1 = bc1_matcher.align_choices(bc1_choices)

    N[f"BC1_score_{score1}"] += 1
    if BC1 == NO_CALL:
        N["BC1_ambig"] += 1
    elif score1 < 2 * len(BC1) * threshold:
        N["BC1_low_score"] += 1
        BC1 = NO_CALL
    else:
        N["BC1_assigned"] += 1

    return bc1, BC1, ref1, score1


def match_BC2(bc2_matcher, seq, qend, tend, N, threshold=0.5):
    bc2_choices = [
        seq[qend:],  # assume local alignment catches entire opseq
        seq[qend - 1 :],  # assume the last opseq base is in fact BC
        # (This indeed happens quite a bit for BC2)
    ]
    # if the end-gap is in fact a mutation, BC would start later
    if tend:
        mut = seq[qend + tend :]
        if len(mut) >= 8:
            bc2_choices.append(mut)

    # print(f"bc2 choices {bc2_choices}")
    (BC2, ref2, score2), bc2 = bc2_matcher.align_choices(bc2_choices)

    N[f"BC2_score_{score2}"] += 1
    if BC2 == NO_CALL:
        N["BC2_ambig"] += 1
    elif score2 < 2 * len(BC2) * threshold:
        N["BC2_low_score"] += 1
        BC2 = NO_CALL
    else:
        N["BC2_assigned"] += 1

    return bc2, BC2, ref2, score2


def opseq_local_align(
    seq,
    opseq="GAATCACGATACGTACACCAGT",
    min_opseq_score=22,
    min_start=8,
    max_end=22 + 12,
    allow_start_gap=False,
    allow_end_gap=False,
):
    results = pairwise2.align.localmd(
        opseq, seq, 2, -1.5, -3, -1, -3, -1, one_alignment_only=True, score_only=False
    )
    tstart = 0
    tend = 0
    if len(results) < 1:
        res = None
    else:
        res = results[0]
        qstart = res.start
        qend = res.end
        L = len(res.seqB)
        # print(res)
        # print(
        #     f"flags {qstart < min_start} {qend > max_end} {res.score < min_opseq_score}"
        # )
        if (
            (min_start and (qstart < min_start))
            or (max_end and (qend > max_end))
            or (res.score < min_opseq_score)
        ):
            res = None
        else:
            # check for start/end gaps
            # backtrack from qstart in seqA until we hit '-'
            if not allow_start_gap:
                while res.seqA[qstart - tstart - 1] != "-":
                    if qstart - tstart - 1 <= 0:
                        tstart = -1
                        res = None
                        break
                    tstart += 1

            if res is not None and not allow_end_gap:
                # and scan forward from qend until we hit '-'
                while res.seqA[qend + tend] != "-":
                    if qend + tend >= len(res.seqA) - 8:
                        tend = -1
                        res = None
                        break
                    tend += 1
    # print("end results", res, tstart, tend)
    return res, tstart, tend


def process_ordered_results(res_queue, args, Qerr, abort_flag):
    with ExceptionLogging("collector", Qerr=Qerr, exc_flag=abort_flag) as el:
        import heapq
        import time

        heap = []
        n_chunk_needed = 0
        t0 = time.time()
        t1 = t0
        n_rec = 0

        logger = el.logger
        out = Output(args)

        for n_chunk, results in queue_iter(res_queue, abort_flag):
            heapq.heappush(heap, (n_chunk, results))

            # as long as the root of the heap is the next needed chunk
            # pass results on to storage
            while heap and (heap[0][0] == n_chunk_needed):
                n_chunk, results = heapq.heappop(heap)  # retrieves heap[0]
                for assigned, record in results:
                    # print("record in process_ordered_results", record)
                    out.write(assigned, record)
                    n_rec += 1

                n_chunk_needed += 1

            # debug output on average throughput
            t2 = time.time()
            if t2 - t1 > 30:
                dT = t2 - t0
                rate = n_rec / dT
                logger.info(
                    "processed {0} reads in {1:.0f} seconds (average {2:.3f} reads/second).".format(
                        n_rec, dT, rate
                    )
                )
                t1 = t2

        out.close()
        # by the time None pops from the queue, all chunks
        # should have been processed!
        if not abort_flag.value:
            assert len(heap) == 0
        else:
            logger.warning(
                f"{len(heap)} chunks remained on the heap due to missing data upon abort."
            )

        dT = time.time() - t0
        logger.info(
            "finished processing {0} reads in {1:.0f} seconds (average {2:.3f} reads/second)".format(
                n_rec, dT, n_rec / dT
            )
        )


def process_fastq(Qfq, args, Qerr, abort_flag):
    """
    reads from two fastq files, groups the input into chunks for
    faster parallel processing, and puts these on a mp.Queue()
    """
    with ExceptionLogging("dispatcher", Qerr=Qerr, exc_flag=abort_flag) as el:
        for chunk in chunkify(read_source(args)):
            logging.debug(f"placing {chunk[0]} {len(chunk[1])} in queue")
            if put_or_abort(Qfq, chunk, abort_flag):
                el.logger.warning("shutdown flag was raised!")
                break


def count_dict_sum(sources):
    dst = defaultdict(float)
    for src in sources:
        for k, v in src.items():
            dst[k] += v

    return dst


def dict_merge(sources):
    dst = {}
    for src in sources:
        dst.update(src)
    return dst


def process_combinatorial(Qfq, Qres, args, Qerr, abort_flag, stat_lists):
    with ExceptionLogging("worker", Qerr=Qerr, exc_flag=abort_flag) as el:
        el.logger.debug(
            f"process_combinatorial starting up with Qfq={Qfq}, Qres={Qres} and args={args}"
        )
        bc1_matcher = TieBreaker(args.bc1_ref, place="left")
        bc2_matcher = TieBreaker(args.bc2_ref, place="right")
        bc1_matcher.load_cache(args.bc1_cache)
        bc2_matcher.load_cache(args.bc2_cache)

        out = Output(args, open_files=False)
        N = defaultdict(int)
        for n_chunk, reads in queue_iter(Qfq, abort_flag):
            el.logger.debug(f"received chunk {n_chunk} of {len(reads)} reads")
            results = []
            for fqid, r1, fqid2, r2, qual2 in reads:
                N["total"] += 1
                out_d = dict(qname=fqid, r1=r1, r2=r2, r2_qual=qual2, r2_qname=fqid2)
                # fallback values for bc1/bc2 so that some BC diversity
                # is maintained for debugging purposes in case we can not
                # assign a decent & unambiguous match
                # lower case bc is for original, uncorrected sequence
                out_d["bc1"] = r1[:12]
                out_d["bc2"] = r2[-12:]
                # upper case BC is for assigned, corrected sequence
                out_d["BC1"] = args.na
                out_d["BC2"] = args.na

                # align opseq sequence to seq of read1
                aln = opseq_local_align(
                    r1.rstrip(),
                    opseq=args.opseq,
                    min_opseq_score=args.min_opseq_score,
                    allow_end_gap=True,  # TODO more permanent fix for this quick'n'dirty hack to get short illumina read to work
                )
                res, tstart, tend = aln
                # print("OPSEQ", res, tstart, tend)
                if res is None:
                    N["opseq_broken"] += 1
                    assigned = False
                else:
                    # identify barcodes
                    bc1, BC1, ref1, score1 = match_BC1(
                        bc1_matcher,
                        res.seqB,
                        res.start,
                        tstart,
                        N,
                        threshold=args.threshold,
                    )
                    # bc2, BC2, ref2, score2 = match_BC2(
                    #     bc2_matcher, res.seqB, res.end, tend, N, threshold=args.threshold
                    # )
                    bc2, BC2, ref2, score2 = "na", "NA", "na", -1

                    # slo = sQSeq.lower()
                    # sout = slo[:qstart] + sQSeq[qstart:qend] + slo[qend:]
                    # print(sout, qstart, qend, tstart, tend, bc1, bc2)
                    # print(f"bc1: {bc1} -> {ref1} -> {BC1} score={50.0*score1/len(bc1):.1f} %")
                    # print(f"bc2: {bc2} -> {ref2} -> {BC2} score={50.0*score2/len(bc2):.1f} %")

                    # best matching pieces of sequence
                    out_d["bc1"] = bc1
                    out_d["bc2"] = bc2
                    # best attempt at assignment
                    out_d["BC1"] = BC1
                    out_d["BC2"] = BC2

                    if BC1 != NO_CALL:  # and BC2 != NO_CALL:
                        N["called"] += 1
                        assigned = True
                    else:
                        assigned = False

                out_d["assigned"] = assigned
                rec = out.make_record(**out_d)
                results.append((assigned, rec))

            Qres.put((n_chunk, results))

        N["BC1_cache_hit"] = bc1_matcher.n_hit
        N["BC2_cache_hit"] = bc2_matcher.n_hit

        # return our counts and observations
        # via these list proxies
        el.logger.debug("synchronizing cache and counts")
        Ns, qcaches1, qcaches2, qcounts1, qcounts2, bccounts1, bccounts2 = stat_lists
        Ns.append(N)
        qcaches1.append(bc1_matcher.cache)
        qcaches2.append(bc2_matcher.cache)
        qcounts1.append(bc1_matcher.query_count)
        qcounts2.append(bc2_matcher.query_count)
        bccounts1.append(bc1_matcher.bc_count)
        bccounts2.append(bc2_matcher.bc_count)


def log_qerr(qerr):
    "helper function for reporting errors in sub processes"
    for name, lines in qerr:
        for line in lines:
            logging.error(f"subprocess {name} exception {line}")


def main_combinatorial(args):

    # queues for communication between processes
    Qfq = mp.Queue(
        args.parallel * 5
    )  # FASTQ reads from process_fastq->process_combinatorial
    Qres = mp.Queue()  # extracted BCs from process_combinatorial->collector
    Qerr = mp.Queue()  # child-processes can report errors back to the main process here

    # Proxy objects to allow workers to report statistics about the run
    manager = mp.Manager()
    abort_flag = mp.Value("b")
    abort_flag.value = False
    qcaches1 = manager.list()
    qcaches2 = manager.list()
    qcounts1 = manager.list()
    qcounts2 = manager.list()
    bccounts1 = manager.list()
    bccounts2 = manager.list()
    Ns = manager.list()
    stat_lists = [Ns, qcaches1, qcaches2, qcounts1, qcounts2, bccounts1, bccounts2]

    with ExceptionLogging("main_combinatorial", exc_flag=abort_flag) as el:

        # read FASTQ in chunks and put them in Qfq
        dispatcher = mp.Process(
            target=process_fastq, name="dispatcher", args=(Qfq, args, Qerr, abort_flag)
        )

        dispatcher.start()
        el.logger.info("Started dispatch")

        # workers consume chunks of FASTQ from Qfq,
        # process them, and put the results in Qres
        workers = []
        for i in range(args.parallel):
            w = mp.Process(
                target=process_combinatorial,
                name=f"worker_{i}",
                args=(Qfq, Qres, args, Qerr, abort_flag, stat_lists),
            )
            w.start()
            workers.append(w)

        el.logger.info("Started workers")
        collector = mp.Process(
            target=process_ordered_results,
            name="output",
            args=(Qres, args, Qerr, abort_flag),
        )
        collector.start()
        el.logger.info("Started collector")
        # wait until all sequences have been thrown onto Qfq
        qfq, qerr = join_with_empty_queues(dispatcher, [Qfq, Qerr], abort_flag)
        el.logger.info("The dispatcher exited")
        if qfq or qerr:
            el.logger.info(f"{len(qfq)} chunks were drained from Qfq upon abort.")
            log_qerr(qerr)

        # signal all workers to finish
        el.logger.info("Signalling all workers to finish")
        for n in range(args.parallel):
            Qfq.put(None)  # each worker consumes exactly one None

        for w in workers:
            # make sure all results are on Qres by waiting for
            # workers to exit. Or, empty queues if aborting.
            qres, qerr = join_with_empty_queues(w, [Qres, Qerr], abort_flag)
            if qres or qerr:
                el.logger.info(f"{len(qres)} chunks were drained from Qres upon abort.")
                log_qerr(qerr)

        el.logger.info(
            "All worker processes have joined. Signalling collector to finish."
        )
        # signal the collector to stop
        Qres.put(None)

        # and wait until all output has been generated
        collector.join()
        el.logger.info("Collector has joined. Merging worker statistics.")

        N = count_dict_sum(Ns)
        if N["total"]:
            el.logger.info(
                f"Run completed. Overall combinatorial barcode assignment "
                f"rate was {100.0 * N['called']/N['total']}"
            )
        else:
            el.logger.error("No reads were processed!")

        if args.update_cache:
            qcount1 = count_dict_sum(qcounts1)
            qcount2 = count_dict_sum(qcounts2)
            cache1 = dict_merge(qcaches1)
            cache2 = dict_merge(qcaches2)
            store_cache(args.bc1_cache, cache1, qcount1)
            store_cache(args.bc2_cache, cache2, qcount2)

        if args.save_stats:
            bccount1 = count_dict_sum(bccounts1)
            bccount2 = count_dict_sum(bccounts2)
            with open(args.save_stats, "w") as f:
                for k, v in sorted(N.items()):
                    f.write(f"freq\t{k}\t{v}\t{100.0 * v/max(N['total'], 1):.2f}\n")

                for k, v in sorted(bccount1.items()):
                    f.write(
                        f"BC1\t{k}\t{v}\t{100.0 * v/max(bccount1['total'], 1):.2f}\n"
                    )

                for k, v in sorted(bccount2.items()):
                    f.write(
                        f"BC2\t{k}\t{v}\t{100.0 * v/max(bccount2['total'], 1):.2f}\n"
                    )


def process_dropseq(Qfq, Qres, args, Qerr, abort_flag, stat_lists):
    with ExceptionLogging("worker", Qerr=Qerr, exc_flag=abort_flag) as el:
        el.logger.debug(
            f"process_dropseq starting up with Qfq={Qfq}, Qres={Qres} and args={args}"
        )
        out = Output(args, open_files=False)
        N = defaultdict(int)
        for n_chunk, reads in queue_iter(Qfq, abort_flag):
            el.logger.debug(f"received chunk {n_chunk} of {len(reads)} reads")
            results = []
            for fqid, r1, fqid2, r2, qual2 in reads:
                N["total"] += 1
                rec = out.make_record(
                    assigned=True,
                    qname=fqid,
                    r1=r1,
                    r2=r2,
                    r2_qual=qual2,
                    r2_qname=fqid2,
                )
                results.append((True, rec))

            Qres.put((n_chunk, results))

        # return our counts and observations
        # via these list proxies
        el.logger.debug("synchronizing cache and counts")
        Ns = stat_lists[0]
        Ns.append(N)


def main_dropseq(args):
    # queues for communication between processes
    Qfq = mp.Queue(
        args.parallel * 5
    )  # FASTQ reads from process_fastq->process_combinatorial
    Qres = mp.Queue()  # extracted BCs from process_combinatorial->collector
    Qerr = mp.Queue()  # child-processes can report errors back to the main process here

    # Proxy objects to allow workers to report statistics about the run
    manager = mp.Manager()
    abort_flag = mp.Value("b")
    abort_flag.value = False
    Ns = manager.list()
    stat_lists = [Ns]

    with ExceptionLogging("main_combinatorial", exc_flag=abort_flag) as el:

        # read FASTQ in chunks and put them in Qfq
        dispatcher = mp.Process(
            target=process_fastq, name="dispatcher", args=(Qfq, args, Qerr, abort_flag)
        )

        dispatcher.start()
        el.logger.info("Started dispatch")

        # workers consume chunks of FASTQ from Qfq,
        # process them, and put the results in Qres
        workers = []
        for i in range(args.parallel):
            w = mp.Process(
                target=process_dropseq,
                name=f"worker_{i}",
                args=(Qfq, Qres, args, Qerr, abort_flag, stat_lists),
            )
            w.start()
            workers.append(w)

        el.logger.info("Started workers")
        collector = mp.Process(
            target=process_ordered_results,
            name="output",
            args=(Qres, args, Qerr, abort_flag),
        )
        collector.start()
        el.logger.info("Started collector")
        # wait until all sequences have been thrown onto Qfq
        qfq, qerr = join_with_empty_queues(dispatcher, [Qfq, Qerr], abort_flag)
        el.logger.info("The dispatcher exited")
        if qfq or qerr:
            el.logger.info(f"{len(qfq)} chunks were drained from Qfq upon abort.")
            log_qerr(qerr)

        # signal all workers to finish
        el.logger.info("Signalling all workers to finish")
        for n in range(args.parallel):
            Qfq.put(None)  # each worker consumes exactly one None

        for w in workers:
            # make sure all results are on Qres by waiting for
            # workers to exit. Or, empty queues if aborting.
            qres, qerr = join_with_empty_queues(w, [Qres, Qerr], abort_flag)
            if qres or qerr:
                el.logger.info(f"{len(qres)} chunks were drained from Qres upon abort.")
                log_qerr(qerr)

        el.logger.info(
            "All worker processes have joined. Signalling collector to finish."
        )
        # signal the collector to stop
        Qres.put(None)

        # and wait until all output has been generated
        collector.join()
        el.logger.info("Collector has joined. Merging worker statistics.")

        N = count_dict_sum(Ns)
        if N["total"]:
            el.logger.info(f"Run completed, {N['total']} reads processed.")
        else:
            el.logger.error("No reads were processed!")

        # N = defaultdict(int)
        # out = Output(args)

        # # iterate query sequences
        # for n, (qname, r1, r2_qname, r2, r2_qual) in enumerate(read_source(args)):
        #     N["total"] += 1
        #     out.write(True, qname=qname, r1=r1, r2_qname=r2_qname, r2=r2, r2_qual=r2_qual)

        #     if n and n % 100000 == 0:
        #         report_stats(N)

        # out.close()

        # report_stats(N)
        if args.save_stats:
            with open(args.save_stats, "w") as f:
                for k, v in sorted(N.items()):
                    f.write(f"freq\t{k}\t{v}\t{100.0 * v/max(N['total'], 1):.2f}\n")

    return N


class Output:
    def __init__(self, args, open_files=True):
        assert Output.safety_check_eval(args.cell_raw)
        assert Output.safety_check_eval(args.cell)
        assert Output.safety_check_eval(args.UMI)
        self.cell_raw = args.cell_raw
        self.cell = args.cell
        self.UMI = args.UMI
        self.na = args.na

        # precompile functions for speed-up
        self.f_cell_raw = compile(self.cell_raw, "<string cell_raw>", "eval")
        self.f_cell = compile(self.cell, "<string cell>", "eval")
        self.f_UMI = compile(self.UMI, "<string UMI>", "eval")

        self.fq_qual = args.fq_qual
        self.bc_na = args.na

        self.tags = []
        for tag in args.bam_tags.split(","):
            self.tags.append(tag.split(":"))

        if args.out_format == "fastq":
            fopen = lambda x: open(x, "w")
            self._write_record = self.write_fastq
            self._make_record = self.make_fastq_record

        elif args.out_format == "bam":
            prog = os.path.basename(__file__)
            header = {
                "HD": {"VN": "1.6"},
                "PG": [
                    {
                        "ID": 0,
                        "PN": prog,
                        "CL": " ".join(sys.argv[1:]),
                        "VN": __version__,
                    },
                ],
                "RG": [
                    {"ID": "A", "SM": args.sample},
                    {"ID": "U", "SM": f"unassigned_{args.sample}"},
                ],
            }
            self.bam_header = pysam.AlignmentHeader.from_dict(header)
            fopen = lambda x: pysam.AlignmentFile(x, "wbu", header=header)
            self._write_record = self.write_bam
            self._make_record = self.make_bam_record
        else:
            raise ValueError(f"unsopported output format '{args.out_format}'")

        if open_files:
            self.out_assigned = fopen(args.out_assigned)
            if args.out_unassigned != args.out_assigned:
                self.out_unassigned = fopen(args.out_unassigned)
            else:
                self.out_unassigned = self.out_assigned

    @staticmethod
    def safety_check_eval(s, danger="();."):
        chars = set(list(s))
        if chars & set(list(danger)):
            return False
        else:
            return True

    def make_bam_record(self, **kw):
        # sys.stderr.write(f"r2_qual={r2_qual}\n")
        a = pysam.AlignedSegment(self.bam_header)
        # STAR does not like spaces in read names so we have to split
        a.query_name = kw["r2_qname"].split()[0]
        a.query_sequence = kw["r2"]
        a.flag = 4
        a.query_qualities = pysam.qualitystring_to_array(kw["r2_qual"])
        a.tags = [(name, templ.format(**kw)) for name, templ in self.tags]
        return a.to_string()

    def make_fastq_record(self, **kw):
        seq = kw["cell"] + kw["UMI"]
        qual = self.fq_qual * len(seq)
        return f"@{kw['qname']}\n{seq}\n+\n{qual}\n"

    def make_record(self, assigned=True, **kw):
        kw["raw"], kw["cell"], kw["UMI"] = self.format(**kw)
        kw["assigned"] = "A" if assigned else "U"
        return self._make_record(**kw)

    def write_bam(self, out, rec):
        out.write(pysam.AlignedSegment.fromstring(rec, self.bam_header))

    def write_fastq(self, out, rec):
        out.write(rec)

    def write(self, assigned, record):
        if assigned:
            self._write_record(self.out_assigned, record)
        else:
            self._write_record(self.out_unassigned, record)

    def format(
        self,
        qname="qname1",
        r2_qname="qname2",
        r2_qual="",
        bc1=None,
        bc2=None,
        BC1=None,
        BC2=None,
        r1="",
        r2="",
        **kw,
    ):
        if BC1 is None:
            BC1 = self.na

        if BC2 is None:
            BC2 = self.na

        # slightly concerned about security here...
        # at least all () and ; raise an assertion in __init__
        cell = eval(self.f_cell)
        raw = eval(self.f_cell_raw)
        UMI = eval(self.f_UMI)
        # cell = eval(self.cell)
        # raw = eval(self.cell_raw)
        # UMI = eval(self.UMI)

        if (cell is None) or (UMI is None):
            raise ValueError(
                f"one of cell=eval('{self.cell}')='{cell}' "
                f"UMI=eval('{self.UMI}')='{UMI}' evaluated to None"
            )

        return raw, cell, UMI

    def close(self):
        self.out_assigned.close()
        self.out_unassigned.close()


# class Process(mp.Process):
#     def __init__(self, *argc, **kw):
#         mp.Process.__init__(self, *argc, **kw)

#     def run(self):


def setup_logging(args):
    FORMAT = f"%(asctime)-20s\t%(levelname)s\t{args.sample}\t%(name)s\t%(message)s"
    # formatter = logging.Formatter(FORMAT)
    lvl = getattr(logging, args.log_level, logging.INFO)
    logging.basicConfig(level=lvl, format=FORMAT)

    fh = logging.FileHandler(filename=args.log_file, mode="a")
    fh.setFormatter(logging.Formatter(FORMAT))
    root = logging.getLogger("")
    root.addHandler(fh)

    logger = logging.getLogger("main")
    logger.info(f"starting raw read preprocessing run")
    for k, v in sorted(vars(args).items()):
        logger.info(f"cmdline arg\t{k}={v}")

    # TODO: store options in YAML for quick
    # re-run with same options
    # {sorted(vars(args).items())}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert raw reads1 and reads2 FASTQ into a single BAM file with cell barcode and UMI as BAM-tags"
    )
    parser.add_argument(
        "--sample",
        default="NA",
        help="source from where to get read1 (FASTQ format)",
    )
    parser.add_argument(
        "--read1",
        default="/dev/stdin",
        help="source from where to get read1 (FASTQ format)",
    )
    parser.add_argument(
        "--read2", default="", help="source from where to get read2 (FASTQ format)"
    )
    parser.add_argument("--cell", default="r1[8:20][::-1]")
    parser.add_argument("--cell-raw", default="None")
    parser.add_argument("--UMI", default="r1[0:8]")
    parser.add_argument(
        "--out-format",
        default="bam",
        choices=["fastq", "bam"],
        help="'bam' for tagged, unaligned BAM, 'fastq' for preprocessed read1 as FASTQ",
    )
    parser.add_argument(
        "--out-assigned",
        default="/dev/stdout",
        help="output for successful assignments (default=/dev/stdout) ",
    )
    parser.add_argument(
        "--out-unassigned",
        default="/dev/stdout",
        help="output for un-successful assignments (default=/dev/stdout) ",
    )
    parser.add_argument(
        "--save-stats",
        default="preprocessing_stats.txt",
        help="store statistics in this file",
    )
    parser.add_argument(
        "--log-file",
        default="preprocessing_run.log",
        help="store all log messages in this file (default=None)",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="set sensitivity of the python logging facility (default=INFO)",
    )
    parser.add_argument(
        "--na", default="NNNNNNNN", help="code for ambiguous or unaligned barcode"
    )
    parser.add_argument(
        "--fq-qual",
        default="E",
        help="phred qual for assigned barcode bases in FASTQ output (default='E')",
    )

    # combinatorial barcode specific options
    parser.add_argument(
        "--parallel", default=1, type=int, help="how many processes to spawn"
    )
    parser.add_argument(
        "--opseq",
        default="GAATCACGATACGTACACCAGT",
        help="opseq primer site sequence",
    )
    parser.add_argument(
        "--min-opseq-score",
        default=22,
        type=float,
        help="minimal score for opseq alignment (default 22 [half of max])",
    )
    parser.add_argument(
        "--threshold",
        default=0.5,
        type=float,
        help="score threshold for calling a match (rel. to max for a given length)",
    )
    parser.add_argument(
        "--bc1-ref", default="", help="load BC1 reference sequences from this FASTA"
    )
    parser.add_argument(
        "--bc2-ref", default="", help="load BC2 reference sequences from this FASTA"
    )
    parser.add_argument(
        "--bc1-cache", default="", help="load cached BC1 alignments from here"
    )
    parser.add_argument(
        "--bc2-cache", default="", help="load cached BC2 alignments from here"
    )
    parser.add_argument(
        "--update-cache",
        default=False,
        action="store_true",
        help="update cache when run is complete (default=False)",
    )
    # parser.add_argument("--template", default="{cell}{UMI}")

    # bam r2 tagging. SAM standard now supports CB=corrected cell barcode, CR=original cell barcode, and MI=molecule identifier/UMI
    parser.add_argument(
        "--bam-tags",
        default="CB:{cell},MI:{UMI},RG:{assigned}",
        help="raw, uncorrected cell barcode",
    )
    args = parser.parse_args()

    return args


def cmdline():
    with ExceptionLogging("main"):
        args = parse_args()
        NO_CALL = args.na
        setup_logging(args)

        if args.out_format == "bam" and not args.read2:
            raise ValueError("bam output format requires --read2 parameter")

        if ("bc1" in args.cell and not args.bc1_ref) or (
            "bc2" in args.cell and not args.bc2_ref
        ):
            raise ValueError(
                "bc1/2 are referenced in --cell or --cell-raw, but no reference barcodes are specified via --bc{{1,2}}-ref"
            )

        if args.bc1_ref or args.bc2_ref:
            main_combinatorial(args)
        else:
            main_dropseq(args)


if __name__ == "__main__":
    cmdline()
