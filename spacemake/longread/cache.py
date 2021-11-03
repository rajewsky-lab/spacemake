__author__ = ["Marvin Jens"]
__license__ = "GPL"
__email__ = ["marvin.jens@mdc-berlin.de"]

import shelve
import os
import logging
import numpy as np
import pandas as pd
from Bio import pairwise2
from spacemake.util import ensure_path, read_fq

"""
A module that uses Smith & Waterman alignment from Biopython to detect
matches of pre-defined sequence blocks (oligos) in PacBio long reads.
Matches are stored using the shelve module for fast access in in-depth
analyses (CachedAlignments class).
Overlapping alignments and co-occurence pattern annotation are handled by
the MultiAlignments class.
"""


def align(ref, query, aln_scores=(2, -1.5, -3, -1, -3, -1), min_score=0.4):
    L = min(len(ref), len(query))
    ms = max(int(min_score * L * aln_scores[0]), 14)  # at least 7 matches
    return [
        r for r in pairwise2.align.localmd(ref, query, *aln_scores) if r.score >= ms
    ]


def print_aln(aln):
    s = pairwise2.format_alignment(*aln, full_sequences=True)
    print(s)
    return s


def non_overlapping_hits(seq, oligo, min_size=7):
    oli_hits = np.zeros(len(seq), dtype=int)
    hit_scores = np.zeros(len(seq))
    #
    for res in align(seq, oligo):
        # print(res.score, res.start, res.end)
        start = res.start
        end = res.end + 1
        scores = hit_scores[start:end]
        hits = oli_hits[start:end]
        oli_hits[start:end] = np.where(scores < res.score, 1, hits)
        hit_scores[start:end] = np.where(scores < res.score, res.score, scores)
    #
    positions = (oli_hits > 0).nonzero()[0]
    if len(positions):
        start = positions[0]
        last = start
        for i in positions[1:]:
            if i > last + 1:
                if last - start > min_size:
                    yield start, last, hit_scores[start:last].max()
                start = i
            last = i
        #
        if last - start > min_size:
            yield start, last, hit_scores[start:last].max()


class CachedAlignments:
    def __init__(self, sample_name, oligo_name, oligo_seq, path="./cache/"):
        """
        Sets up transparent caching of alignments. A shelve DB will be
        created in {path}/{sample_name}__{oligo_name}.<db suffixes>, which allows to
        concurrently align different oligos against the same reads.
        """
        self.logger = logging.getLogger("CachedAlignments")
        self.logger.info(f"initializing for sample={sample_name} oligo={oligo_name}")
        self.sample_name = sample_name
        self.oligo_name = oligo_name
        self.oligo_seq = oligo_seq
        self.oli_lengths = {"NA": 0}
        for oli, seq in zip(self.oligo_name, self.oligo_seq):
            self.oli_lengths[oli] = float(len(seq))

        self.sname = os.path.join(ensure_path(path), f"{sample_name}__{oligo_name}")
        self.shelf = shelve.open(self.sname)
        self.logger.info(f"found {len(self.shelf)} cached entries")

    def align_or_load(self, seq_name, seq):
        if not seq_name in self.shelf:
            hits = list(non_overlapping_hits(seq, self.oligo_seq))
            self.shelf[seq_name] = hits

        return self.shelf[seq_name]

    def close(self):
        self.shelf.sync()
        self.shelf.close()


class MultiAlignments:
    def __init__(self, sample_name, oligo_dict, path="./cache/"):
        """
        Combine alignments of multiple oligos against a set of reads and
        annotate occurence patterns. Caches will be built/opened from
        {path}/{sample_name}__{oligo_name} for all oligos in oligo_dict,
        which maps oligo_name -> oligo_sequence.
        """
        self.logger = logging.getLogger("MultiAlignments")
        self.sample_name = sample_name
        self.oligo_dict = oligo_dict
        self.oligo_names = sorted(oligo_dict.keys())
        self.aln_caches = {}
        for o in self.oligo_names:
            aln = CachedAlignments(sample_name, o, self.oligo_dict[o])
            self.aln_caches[o] = aln

    def annotate(self, seq_name, seq, min_size=7):
        """
        Iterates over all oligo matches for a desired read sequence and yields
        (start, end, label, score) annotations
        """
        labels = [
            "NA",
        ] + self.oligo_names
        L = len(seq)
        oli_hits = np.zeros(L, dtype=int)
        hit_scores = np.zeros(L)

        for i, o in enumerate(self.oligo_names):
            aln = self.aln_caches[o]
            for start, end, score in aln.align_or_load(seq_name, seq):
                scores = hit_scores[start:end]
                hits = oli_hits[start:end]
                oli_hits[start:end] = np.where(scores < score, i + 1, hits)
                hit_scores[start:end] = np.where(scores < score, score, scores)

        positions = (oli_hits > 0).nonzero()[0]
        if len(positions):
            start = positions[0]
            prev = oli_hits[start]
            last = start
            for i in positions[1:]:
                curr = oli_hits[i]
                if curr != prev or i > last + 1:
                    if last - start > min_size:
                        label = labels[prev]
                        score = hit_scores[
                            start:last
                        ].max()  # * (last - start) / self.oli_lengths[label]
                        yield start, last, label, score
                    start = i
                last = i
                prev = curr
            #
            if last - start > min_size:
                label = labels[prev]
                score = hit_scores[
                    start:last
                ].max()  # * (last - start) / self.oli_lengths[label]
                yield start, last, labels[prev], score


def align_one_oligo_one_read(oligo_seq, qname, qseq):
    hits = list(non_overlapping_hits(qseq, oligo_seq))
    return qname, hits


def align_one_oligo(oligo_name, oligo_seq, sample_name, fastq_path, n_proc=None):
    import multiprocessing as mp

    aln = CachedAlignments(sample_name, oligo_name, oligo_seq)
    n = 0
    job_list = []
    for name, seq, qual in read_fq(fastq_path):
        if not name in aln.shelf:
            job_list.append((oligo_seq, name, seq))

        # aln.align_or_load(name, seq)
        n += 1

    if job_list:
        logging.info(
            f"about to align {len(job_list)} reads with oligo {oligo_name} using {n_proc} processes"
        )
        pool = mp.Pool(n_proc)
        for qname, hits in pool.starmap(
            align_one_oligo_one_read, job_list, chunksize=10
        ):
            aln.shelf[qname] = hits
    else:
        logging.info(
            f"every read was already aligned against {oligo_name}. Nothing to do."
        )

    return sample_name, oligo_name, n


def fill_caches(
    fastq_path, sample_name, oligo_dict, path="./cache/", n_proc=None
):
    """
    Utility function to fill the shelves by running alignments of different oligos.
    Not parallelized
    """
    # import multiprocessing as mp
    for oligo_name, oligo_seq in sorted(oligo_dict.items()):
        sample_name, oligo_name, n_reads = align_one_oligo(
            oligo_name, oligo_seq, sample_name, fastq_path, n_proc=n_proc
        )
        print(f"{sample_name}: aligned {n_reads} against {oligo_name}")


def annotate(fastq_path, sample_name, oligo_dict, path="./cache/"):
    multi = MultiAlignments(sample_name, oligo_dict, path=path)
    data = []
    for name, seq, qual in read_fq(fastq_path):
        for start, end, label, score in multi.annotate(name, seq):
            data.append((sample_name, name, len(seq), label, start, end, score))

    df = pd.DataFrame(
        data, columns=["sample_name", "qname", "L", "oligo", "start", "end", "score"]
    )
    return df
