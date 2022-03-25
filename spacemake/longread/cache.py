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


def align_one_oligo_one_read(oligo_seq, qname, qseq):
    # logging.debug(f"{oligo_seq} vs {qname} (len={len(qseq)})")
    hits = list(non_overlapping_hits(qseq, oligo_seq))
    return qname, hits


class CachedAlignments:
    def __init__(self, sample_name, oligo_name, oligo_seq, path="./cache/"):
        """
        Sets up transparent caching of alignments. A DataFrame will be
        created in {path}/{sample_name}__{oligo_name}.tsv.
        """
        self.logger = logging.getLogger("CachedAlignments")
        self.sample_name = sample_name
        self.oligo_name = oligo_name
        self.oligo_seq = oligo_seq
        self.oli_lengths = {"NA": 0}
        for oli, seq in zip(self.oligo_name, self.oligo_seq):
            self.oli_lengths[oli] = float(len(seq))

        self.sname = os.path.join(ensure_path(path), f"{sample_name}__{oligo_name}")
        self.logger.info(
            f"initializing for sample={sample_name} oligo={oligo_name} -> {self.sname}"
        )
        self.df = self.load_df(self.sname + ".tsv")
        self._df_modified = False
        self.logger.info(f"found {len(self.df)} cached entries")
        self.lookup = self.build_lookup(self.df)

    def load_df(self, fname):
        self.logger.info(f"trying to load DataFrame from '{fname}'")
        if not os.access(fname, os.R_OK):
            self.logger.warning(f"'{fname}' not found, looking for shelf-db to convert")
            self.shelf_to_df()

        if not os.access(fname, os.R_OK):
            df = pd.DataFrame(columns=["qname", "start", "end", "score"]).set_index(
                "qname"
            )
            self.logger.info(f"starting '{fname}' from empty dataframe")
        else:
            df = pd.read_csv(fname, sep="\t").set_index("qname")
            self.logger.info(f"loaded {len(df)} rows from '{fname}'")

        return df

    def df_from_hits(self, qname, hits):
        df_new = pd.DataFrame(hits, columns=["start", "end", "score"])
        df_new["qname"] = qname
        return df_new.set_index("qname")

    def build_lookup(self, df):
        from time import time
        from collections import defaultdict

        lkup = defaultdict(list)
        t0 = time()
        for row in df.itertuples():
            # print(row.Index, type(row.Index))
            lkup[row.Index].append((row.start, row.end, row.score))

        dt = time() - t0
        self.logger.info(f"built fast lookup of {len(df)} entries in {dt:.3f} seconds")
        return lkup

    def query(self, qname, qseq):
        if qname in self.lookup:
            hits = self.lookup[qname]
            # # fix pandas stupidity
            # if len(hits.shape) < 2:
            #     hits = hits.reshape((1, 3))
        else:
            hits = []
            # qname, hits = align_one_oligo_one_read(self.oligo_seq, qname, qseq)
            # self.df = self.df.append(self.df_from_hits(qname, hits))
            # self._df_modified = True
            # hits = np.array(hits)
        # print("HITS:", hits)
        return hits

    def align_fastq(self, fastq_path, n_proc=8):
        import multiprocessing as mp

        self.logger.info(
            f"preparing to align {self.oligo_name} against raw reads in '{fastq_path}'"
        )
        n = 0
        job_list = []
        for qname, seq, qual in read_fq(fastq_path):
            if not qname in self.df.index:
                job_list.append((self.oligo_seq, qname, seq))

            n += 1

        if job_list:
            logging.info(
                f"need to align {len(job_list)} reads using {n_proc} processes"
            )
            pool = mp.Pool(n_proc)
            new = [
                self.df_from_hits(qname, hits)
                for qname, hits in pool.starmap(
                    align_one_oligo_one_read, job_list, chunksize=100
                )
            ]

            if len(new):
                self._df_modified = True
                self.df = pd.concat([self.df] + new)

        else:
            logging.info(
                f"every read was already aligned against {self.oligo_name}. Nothing to do."
            )
        self.sync()
        self.lookup = self.build_lookup(self.df)

    def shelf_to_df(self):
        import dbm

        dfname = self.sname + ".tsv"
        data = []
        try:
            shelf = shelve.open(self.sname, flag="r")
        except dbm.error:
            return

        for qname in shelf.keys():
            for start, end, score in shelf[qname]:
                data.append((qname, start, end, score))

        df = pd.DataFrame(data, columns=["qname", "start", "end", "score"]).set_index(
            "qname"
        )
        df.to_csv(dfname, sep="\t")
        return df

    def sync(self):
        dfname = self.sname + ".tsv"
        if self._df_modified:
            self.logger.info(
                f"sync: storing changed dataframe '{dfname}' with {len(self.df)} entries..."
            )
            self.df.to_csv(dfname, sep="\t")
        else:
            self.logger.info(f"sync: dataframe '{dfname}' has not been changed")

        # self.shelf.sync()
        # self.shelf.close()

    @classmethod
    def crawl_and_convert_df(cls, root="./"):
        """
        legacy converter. Convert shelf databases to Pandas DataFrames
        """
        import os

        for (dirpath, dirnames, filenames) in os.walk(root):
            for fname in filenames:
                if fname.endswith(".dat"):
                    print(os.path.join(dirpath, fname))
                    name = os.path.splitext(fname)[0]
                    sample_name, oligo_name = name.split("__")
                    print(sample_name, oligo_name)
                    c = CachedAlignments(sample_name, oligo_name, "NA", path=dirpath)
                    c.to_df()


class MultiAlignments:
    def __init__(
        self, sample_name, oligo_dict, min_length=15, path="./cache/", relevant=[]
    ):
        """
        Combine alignments of multiple oligos against a set of reads and
        annotate occurence patterns. Caches will be built/opened from
        {path}/{sample_name}__{oligo_name} for all oligos in oligo_dict,
        which maps oligo_name -> oligo_sequence.
        """
        self.logger = logging.getLogger("MultiAlignments")
        self.sample_name = sample_name
        self.oligo_dict = oligo_dict
        # sort oligos by length and accept longer oligo match only if score is at least 2 matches better (+4)
        self.relevant = set(relevant)
        for m in relevant:
            #  also include reverse-complement hits of an oligo
            self.relevant.add(f"{m}_RC".replace("_RC_RC", ""))

        self.min_length = min_length

        self.logger.info(
            f"restricting to relevant matches {sorted(self.relevant)} "
            f"and discarding matches under {min_length} bases"
        )

        self.oligo_names = sorted(
            oligo_dict.keys(), key=lambda k: len(self.oligo_dict[k])
        )
        self.oligo_names = [
            o for o in self.oligo_names if o in self.relevant or len(self.relevant) == 0
        ]
        self.aln_caches = {}
        for o in self.oligo_names:
            aln = CachedAlignments(sample_name, o, self.oligo_dict[o], path=path)
            self.aln_caches[o] = aln

    def annotate(self, seq_name, seq, score_gap=3):
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
            # print("INVESTIGATING", o, len(self.oligo_dict[o]))
            for start, end, score in aln.query(seq_name, seq):
                start = int(start)
                end = int(end)
                # print(f"hit: {start}-{end} score={score}")
                scores = hit_scores[start:end]
                # print(f"len(scores)={len(scores)}, end -start={end-start}")
                # print("scores before", scores)

                hits = oli_hits[start:end]
                # print("assignments before", hits)
                # print([self.oligo_names[h] for h in sorted(set(hits))])

                if scores.max() + score_gap < score:
                    oli_hits[start:end] = i + 1
                    hit_scores[start:end] = score
                # oli_hits[start : end + 1] = np.where(
                #     scores + score_gap < score, i + 1, hits
                # )
                # hit_scores[start : end + 1] = np.where(
                #     scores + score_gap < score, score, scores
                # )
                # print("scores after", scores)
                # print("assignments after", hits)
                # print([labels[h] for h in sorted(set(hits))])

        positions = (oli_hits > 0).nonzero()[0]
        if len(positions):
            start = positions[0]
            prev = oli_hits[start]
            last = start
            for i in positions[1:]:
                curr = oli_hits[i]
                if curr != prev or i > last + 1:
                    if last - start > self.min_length:
                        label = labels[prev]
                        score = hit_scores[
                            start:last
                        ].max()  # * (last - start) / self.oli_lengths[label]
                        yield start, last + 1, label, score
                    start = i
                last = i
                prev = curr
            #
            if last - start > self.min_length:
                label = labels[prev]
                score = hit_scores[
                    start:last
                ].max()  # * (last - start) / self.oli_lengths[label]
                yield start, last + 1, labels[prev], score

    def annotate_fastq(self, fastq_path):
        self.logger.info(f"about to annotate reads from '{fastq_path}'")
        data = []
        for i, (name, seq, _) in enumerate(read_fq(fastq_path)):
            # print(i, name)
            for start, end, label, score in self.annotate(name, seq):
                # print("  ", label, start, end, score)
                data.append(
                    (self.sample_name, name, len(seq), label, start, end, score)
                )

            if i and i % 5000 == 0:
                self.logger.info(f"processed {i} reads...")

        df = pd.DataFrame(
            data,
            columns=["sample_name", "qname", "L", "oligo", "start", "end", "score"],
        )
        self.logger.info(f"done. compiled DataFrame with {len(df)} entries.")
        return df


# def align_one_oligo(oligo_name, oligo_seq, sample_name, fastq_path, n_proc=None, **kw):
#     import multiprocessing as mp

#     aln = CachedAlignments(sample_name, oligo_name, oligo_seq, **kw)
#     aln.
#     n = 0
#     job_list = []
#     for name, seq, qual in read_fq(fastq_path):
#         if not name in aln.shelf:
#             job_list.append((oligo_seq, name, seq))

#         # aln.align_or_load(name, seq)
#         n += 1

#     if job_list:
#         logging.info(
#             f"about to align {len(job_list)} reads with oligo {oligo_name} using {n_proc} processes"
#         )
#         pool = mp.Pool(n_proc)
#         for qname, hits in pool.starmap(
#             align_one_oligo_one_read, job_list, chunksize=10
#         ):
#             aln.shelf[qname] = hits
#     else:
#         logging.info(
#             f"every read was already aligned against {oligo_name}. Nothing to do."
#         )

#     return sample_name, oligo_name, n


def rev_comp_name(name):
    return (name + "_RC").replace("_RC_RC", "")


def fill_caches(
    fastq_path, sample_name, oligo_dict, path="./cache/", n_proc=None, relevant=[]
):
    """
    Utility function to run alignments in parallel
    """
    # import multiprocessing as mp
    for oligo_name, oligo_seq in sorted(oligo_dict.items()):
        if relevant:
            if oligo_name not in relevant and rev_comp_name(oligo_name) not in relevant:
                logging.info(
                    f"skipping irrelevant oligo {oligo_name} for this signature"
                )
                continue

        aln = CachedAlignments(sample_name, oligo_name, oligo_seq, path=path)
        aln.align_fastq(fastq_path, n_proc=n_proc)


def annotate(fastq_path, sample_name, oligo_dict, path="./cache/", relevant=[]):
    multi = MultiAlignments(sample_name, oligo_dict, path=path, relevant=relevant)
    return multi.annotate_fastq(fastq_path)
    data = []
    # for name, seq, qual in read_fq(fastq_path):
    #     for start, end, label, score in multi.annotate(name, seq):
    #         data.append((sample_name, name, len(seq), label, start, end, score))

    # df = pd.DataFrame(
    #     data, columns=["sample_name", "qname", "L", "oligo", "start", "end", "score"]
    # )
    # return df


if __name__ == "__main__":
    import spacemake.longread.util as util

    oligo_dict = util.load_oligos("")

    seq = "CGCGGAAGCAGTGGTATCAACGCAGAGTACAAAA"
    print(oligo_dict.keys())
    for oli_name in ["bead_start", "chr_TSO", "chrV2_RT_PRIMER"]:
        oli_seq = oligo_dict[oli_name]
        print(f"oli_name={oli_name} oli_seq={oli_seq}")
        for aln in align(seq, oli_seq):
            print_aln(aln)
            print(aln.start, aln.end)

        print("testing non-overlapping hits")
        for x in non_overlapping_hits(seq, oli_seq):
            print(x)

    m = MultiAlignments("test", oligo_dict)
    for hit in m.annotate(
        "test1",
        seq,
        score_gap=3,
    ):
        print(hit)
        start, end, label, score = hit
        print(seq[start:end])
