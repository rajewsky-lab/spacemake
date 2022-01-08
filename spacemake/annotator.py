import pandas as pd
import numpy as np
import os
import ncls
import argparse
import re
import gzip
import logging
import ncls
import itertools
import pickle
from time import time
from collections import defaultdict, OrderedDict, deque

logging.basicConfig(level=logging.DEBUG)

## reading the GTF
def attr_to_dict(attr_str):
    d = {}
    for match in re.finditer(r"(\w+) \"(\S+)\";", attr_str):
        key, value = match.groups()
        d[key] = value

    return d


def load_GTF(
    src,
    attributes=["gene_id", "gene_type", "gene_name"],
    features=["transcript", "exon", "CDS", "UTR"],
):
    if type(src) is str:
        if src.endswith(".gz"):
            src = gzip.open(src, "rt")
        else:
            src = open(src, "rt")

    fset = set(features)
    data = []
    for line in src:
        parts = line.split("\t")
        if len(parts) != 9:
            # not a GFF formatted line
            continue

        chrom, source, feature, start, end, score, strand, frame, attr_str = parts
        if feature not in fset:
            continue

        attrs = attr_to_dict(attr_str)
        row = [chrom, feature, int(start) - 1, int(end), strand] + [
            attrs[a] for a in attributes
        ]
        data.append(row)

    return pd.DataFrame(
        data, columns=["chrom", "feature", "start", "end", "strand"] + attributes
    ).drop_duplicates()


## Annotation classification matrix
default_map = {
    # CDS   UTR    exon   transcript
    (True, False, True, True): "CDS",
    (True, True, True, True): "CDS",
    (False, True, True, True): "UTR",
    (False, False, True, True): "non-coding",
    (True, False, False, True): "intron",
    (False, True, False, True): "intron",
    (False, False, False, True): "intron",
}


class Classifier:
    def __init__(
        self,
        df,
        hierarchy=["CDS", "UTR", "non-coding", "intron"],
        feature_map=default_map,
    ):
        self.feature_map = feature_map
        self.hierarchy = hierarchy
        self.prio = {}
        for i, h in enumerate(hierarchy):
            self.prio[h] = i

        feat2idx = {
            "CDS": 0,
            "UTR": 1,
            "exon": 2,
            "transcript": 3,
        }
        df["feature_idx"] = df["feature"].apply(lambda x: feat2idx[x])
        self.df = df.values

    def process(self, ids):
        gene_features = defaultdict(lambda: np.zeros(4))
        gene_names = {}
        gene_types = {}
        # iterate over the overlapping GTF features only once.
        # sort by gene_id as we go
        for i in ids:
            row = self.df[i][5:]
            (gene_id, gene_type, gene_name, feature_idx) = row
            gene_names[gene_id] = gene_name
            gene_types[gene_id] = gene_type
            # feature_idx encodes CDS, UTR, exon, transcript records
            # as numbers to quickly construct a correct absence/presence
            # boolean flag-vector to be used for the lookup
            gene_features[gene_id][feature_idx] = True

        top_prio = np.inf
        gn_by_prio = defaultdict(list)
        gf_by_prio = defaultdict(list)
        gt_by_prio = defaultdict(list)

        for gene_id, feature_flags in gene_features.items():
            # here the boolean vector is converted to a tuple
            # which can serve as a key for lookup
            cls = self.feature_map[tuple(feature_flags)]
            # print(tx_id, gene, gene_id, gene_type, features, cls)
            prio = self.prio[cls]
            top_prio = min(top_prio, prio)
            gn_by_prio[prio].append(gene_names[gene_id])
            gf_by_prio[prio].append(cls)
            gt_by_prio[prio].append(gene_types[gene_id])

        return gf_by_prio[top_prio], gn_by_prio[top_prio], gt_by_prio[top_prio]


## Helper functions for working with NCLS
def query(nc, x0, x1):
    """
    report only the target_ids from the overlapping
    (start, end, id) tuples returned by NCLS.find_overlap().
    format is frozenset bc that can be hashed and used as a key.
    """
    return frozenset((o[2] for o in nc.find_overlap(x0, x1)))


def decompose(nc):
    """
    Iterate over the unique start or end coordinates of all the intervals
    in a NCLS, essentially decomposing the original NCLS into non-overlapping
    segments of unique target_id combinations.
    This is used to *compile* an NCLS. Compiled, here means that target_ids
    reference not original GTF features, but enumerated, unique combinations of
    the GTF features as they are positioned and potentially overlap in the
    reference sequence.
    The advantage is that each such combination has to be processed into the
    tagging information only once, and that this processing can be done prior
    to any tagging, reducing the complexity to a simple array lookup with O(1).

    :param nc: The NCLS (nested list) with original genomic features from GTF
    :type nc: NCLS64
    :yield: (start, end, frozenset(target_ids) )
    :rtype: None
    """
    logger = logging.getLogger("compile")
    logger.info(f"compiling strand into non-overlapping and pre-classified annotations")

    starts, ends, ids = np.array(nc.intervals()).T
    logger.debug(f"retrieved {len(starts)} intervals from ncls")

    breakpoints = np.array(sorted(set(starts) | set(ends)))
    logger.debug(f"identified {len(breakpoints)} breakpoints")

    if not len(starts):
        return

    # since we sorted, the first must be a start position
    last_pos = breakpoints[0]
    last_key = query(nc, last_pos, last_pos + 1)
    for bp in breakpoints[1:]:
        # depending on wether this is another start or
        # end position, the transition point can be off by one nt
        # in either direction. Most robust fix was to scan all 3
        # options until the key changes (while ensuring that end > start)
        for x in [-1, 0, +1]:
            if bp + x <= last_pos:
                continue

            key = query(nc, bp + x, bp + x + 1)
            if key != last_key:
                if len(last_key):
                    # ensure we do not yield the empty set
                    yield last_pos, bp + x, last_key
                last_key = key
                last_pos = bp + x
                break


## NCLS-based lookup of genome annotation
class GenomeAnnotation:
    def __init__(self, df, processor):
        self.logger = logging.getLogger("GenomeAnnotation")
        self.df = df
        self.processor = processor

        # find all unique combinations of chrom + strand
        strands = df[["chrom", "strand"]].drop_duplicates().values
        self.strand_keys = [tuple(s) for s in strands]
        self.strand_map = {}

        t0 = time()
        for strand_key in sorted(self.strand_keys):
            chrom, strand = strand_key

            d = self.df.query(f"chrom == '{chrom}' and strand == '{strand}'")
            nested_list = ncls.NCLS(d["start"], d["end"], d.index)
            self.strand_map[strand_key] = nested_list

        dt = time() - t0
        self.logger.info(
            f"constructed nested lists of {len(df)} features on {len(self.strand_keys)} strands in {dt:.3f}s"
        )

    def query_idx(self, chrom, start, end, strand):
        strand_key = (chrom, strand)
        if not strand_key in self.strand_map:
            return []

        nested_list = self.strand_map[strand_key]
        return query(nested_list, start, end)

    def query_idx_blocks(self, chrom, strand, blocks):
        idx = set()
        for start, end in blocks:
            idx |= self.query_idx(chrom, start, end, strand)
        return frozenset(idx)

    def query(self, chrom, start, end, strand):
        idx = self.query_idx(chrom, start, end, strand)
        return self.processor(idx)

    def query_blocks(self, chrom, strand, blocks):
        idx = self.query_idx_blocks(chrom, strand, blocks)
        return self.processor(idx)

    def compile(self):
        chroms = []
        strands = []
        cstarts = []
        cends = []
        cids = []
        cidx = []
        cid_lkup = {}

        for strand_key, nc in self.strand_map.items():
            chrom, strand = strand_key
            self.logger.info(f"decomposing {chrom} {strand}")
            for start, end, idx in decompose(nc):
                # print(f"start={start} end={end} idx={idx}")
                chroms.append(chrom)
                strands.append(strand)
                cstarts.append(start)
                cends.append(end)
                if idx not in cid_lkup:
                    cid_lkup[idx] = len(cid_lkup)
                    cidx.append(idx)

                cids.append(cid_lkup[idx])

        self.logger.info("done")
        cdf = pd.DataFrame(
            dict(chrom=chroms, strand=strands, start=cstarts, end=cends, cid=cids)
        )
        return cdf, cidx


# The only piece of work left is merging multiple pre-classified annotations
def joiner(idx):
    if len(idx) == 1:
        # fast path
        return classifications[tuple(idx)]

    # we want to report each combination only once
    # Note: This code path is executed rarely
    seen = set()
    gf = []
    gn = []
    gt = []
    for i in idx:
        for ann in zip(*classifications[i]):
            if not ann in seen:
                seen.add(ann)
                gf.append(ann[0])
                gn.append(ann[1])
                gt.append(ann[2])

    return gf, gn, gt


def annotate_BAM(ga, src, out):
    import pysam

    bam = pysam.AlignmentFile(src)
    out = pysam.AlignmentFile(out, "wbu", template=bam)
    t0 = time()
    T = 5
    for n, read in enumerate(bam.fetch(until_eof=True)):
        if not read.is_unmapped:
            chrom = bam.get_reference_name(read.tid)
            strand = "-" if read.is_reverse else "+"
            gf, gn, gt = ga.query_blocks(chrom, strand, read.get_blocks())
            if len(gf):
                read.tags += [
                    ("gf", ",".join(gf)),
                    ("gn", ",".join(gn)),
                    ("gt", ",".join(gt)),
                ]
            else:
                read.tags += [("gf", "INTERGENIC")]

        out.write(read)
        dt = time() - t0
        if dt > T:
            print(f"processed {n} reads in {dt:.2f} seconds ({n/dt:.2f} reads/second)")
            T += 5


def has_compiled_annotation(args):
    cdf_path = os.path.join(args.compiled, "non_overlapping.csv")
    cclass_path = os.path.join(args.compiled, "classifications.npy")

    if os.access(cdf_path, os.R_OK) and os.access(cclass_path, os.R_OK):
        found = True
    else:
        found = False

    return found, cdf_path, cclass_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gtf",
        required=True,
        help="path to the original annotation (e.g. gencodev38.gtf.gz)",
    )
    parser.add_argument(
        "--compiled",
        required=True,
        help="path to keep a compiled version of the GTF in (used as cache)",
    )
    parser.add_argument("--bam-in", help="path for the input BAM to be tagged")
    parser.add_argument(
        "--bam-out",
        help="path for the tagged BAM output",
    )
    args = parser.parse_args()

    found, cdf_path, cclass_path = has_compiled_annotation(args)
    if found:
        ## re-read the compiled index
        t0 = time()
        cdf = pd.read_csv(cdf_path, sep="\t")
        classifications = np.load(cclass_path, allow_pickle=True)
        dt = time() - t0
        print(f"re-read compiled annotation index in {dt:.3f} seconds")

    else:
        # load GTF the first time. Need to build compiled annotation
        t0 = time()
        df = load_GTF(args.gtf)
        dt = time() - t0
        print(f"loaded {len(df)} records in {dt:.3f} seconds")

        ## Build NCLS with original GTF features
        cl = Classifier(df)
        ga = GenomeAnnotation(df, lambda idx: cl.process(idx))
        t0 = time()
        ## Compile into non-overlapping, unique combinations of GTF features
        cdf, idx_combinations = ga.compile()
        dt = time() - t0
        print(
            f"decomposed complete original GTF annotation overlaps in {dt:.3f} seconds"
        )
        # print(cdf, len(idx_combinations))

        ## pre-classify each unique combination of GTF features
        classifications = []
        t0 = time()
        for n, idx in enumerate(idx_combinations):
            res = ga.processor(idx)
            classifications.append(res)

        ## Turn into np.array to turn lookup into super-fast array index operation
        classifications = np.array(classifications, dtype=object)
        dt = time() - t0
        print(f"processed {n} feature combinations in {dt:.3f} seconds")

        ## Store the compiled DataFrame and pre-classifications
        t0 = time()
        cdf.to_csv(cdf_path, sep="\t", index=False)
        np.save(cclass_path, classifications)
        dt = time() - t0
        print(f"stored compiled annotation index in {dt:.3f} seconds")

    ## Create a secondary Annotator which uses the non-overlapping combinations
    ## and the pre-classified annotations for the actual tagging
    gc = GenomeAnnotation(cdf, joiner)

    t0 = time()
    annotate_BAM(gc, args.bam_in, args.bam_out)
    dt = time() - t0
    print(f"{dt:.3f} seconds. ")
