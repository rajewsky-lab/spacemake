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


def attr_to_dict(attr_str):
    """
    Helper for reading the GTF. Parsers the attr_str found in the GTF
    attribute column (col 9). Returns a dictionary from the key, value
    pairs.
    """
    d = {}
    for match in re.finditer(r"(\w+) \"(\S+)\";", attr_str):
        key, value = match.groups()
        d[key] = value

    return d


def load_GTF(
    src,
    attributes=["gene_id", "gene_type", "gene_name"],
    features=["exon", "CDS", "UTR"],
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
    # CDS   UTR    exon -> ga gf tag values
    (True, False, True): "CDS",  # CODING
    (True, True, True): "CDS",  # CODING
    (False, True, True): "UTR",  # UTR
    (False, False, True): "non-coding",  # CODING
    (True, False, False): "intron",  # INTRONIC
    (False, True, False): "intron",  # INTRONIC
    (False, False, False): "intron",  # INTRONIC
}


class GTFClassifier:
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
            # "transcript": 3,
        }
        df["feature_idx"] = df["feature"].apply(lambda x: feat2idx[x])
        self.df = df.values

    def process(self, ids):
        gene_features = defaultdict(lambda: np.zeros(3, dtype=bool))
        gene_names = {}
        gene_types = {}
        # iterate over the overlapping GTF features only once.
        # sort by gene_id as we go
        for i in ids:
            row = self.df[i][5:]
            # print(row)
            (strand, gene_id, gene_type, gene_name, feature_idx) = row
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
        gs_by_prio = defaultdict(list)

        for gene_id, feature_flags in gene_features.items():
            # here the boolean vector is converted to a tuple
            # which can serve as a key for lookup
            if not feature_flags[-1]:
                print(gene_id)
                for i in ids:
                    print(self.df[i])

            cls = self.feature_map[tuple(feature_flags)]
            # print(tx_id, gene, gene_id, gene_type, features, cls)
            prio = self.prio[cls]
            top_prio = min(top_prio, prio)
            gn_by_prio[prio].append(gene_names[gene_id])
            gf_by_prio[prio].append(cls)
            gt_by_prio[prio].append(gene_types[gene_id])
            gs_by_prio[prio].append(strand)

        return (
            gf_by_prio[top_prio],
            gn_by_prio[top_prio],
            gs_by_prio[top_prio],
            gt_by_prio[top_prio],
        )


class CompiledClassifier:
    def __init__(self, df, classifications):
        self.cid_table = df["cid"].values
        self.classifications = np.array(classifications, dtype=object)

    @staticmethod
    def get_filenames(path):
        cdf_path = os.path.join(path, "non_overlapping.csv")
        cclass_path = os.path.join(path, "classifications.npy")

        return cdf_path, cclass_path

    @staticmethod
    def files_exist(path):
        cdf_path, cclass_path = CompiledClassifier.get_filenames(path)
        return os.access(cdf_path, os.R_OK) and os.access(cclass_path, os.R_OK)

    # The only piece of work left is merging multiple pre-classified annotations
    def joiner(self, cidx):
        # we want to report each combination only once
        # Note: This code path is executed rarely
        seen = set()
        gf = []
        gn = []
        gs = []
        gt = []
        for i in tuple(cidx):
            cid = self.cid_table[i]
            for ann in zip(*self.classifications[cid]):
                if not ann in seen:
                    seen.add(ann)
                    gf.append(ann[0])
                    gn.append(ann[1])
                    gs.append(ann[2])
                    gt.append(ann[3])

        return gf, gn, gs, gt

    def process(self, cidx):
        # fast path
        if len(cidx) == 1:
            # fast path
            idx = tuple(cidx)[0]
            return self.classifications[self.cid_table[idx]]

        return self.joiner(cidx)


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


class GenomeAnnotation:
    """
    NCLS-based lookup of genome annotation
    """

    logger = logging.getLogger("GenomeAnnotation")

    def __init__(self, df, processor, is_compiled=False):
        """
        [summary]

        :param df: Annotation data of either raw GTF features or pre-classified combinations
        :type df: DataFrame with at least the following columns: chrom, strand, start, end
        :param processor: evaluates indices into df and returns annotation
            the indices point into the dataframe to all
            overlapping features
        :type processor: function that takes frozenset(indices) as sole argument
        """
        # self.df = df
        self.processor = processor

        # find all unique combinations of chrom + strand
        strands = df[["chrom", "strand"]].drop_duplicates().values
        self.strand_keys = [tuple(s) for s in strands]
        self.strand_map = {}
        self.empty = frozenset([])

        t0 = time()
        for strand_key in sorted(self.strand_keys):
            chrom, strand = strand_key

            d = df.query(f"chrom == '{chrom}' and strand == '{strand}'")
            nested_list = ncls.NCLS(d["start"], d["end"], d.index)
            self.strand_map[strand_key] = nested_list

        dt = time() - t0
        self.logger.info(
            f"constructed nested lists of {len(df)} features on {len(self.strand_keys)} strands in {dt:.3f}s"
        )
        self.is_compiled = is_compiled

    @classmethod
    def from_compiled_index(cls, path):
        cdf_path, cclass_path = CompiledClassifier.get_filenames(path)
        t0 = time()
        cdf = pd.read_csv(cdf_path, sep="\t")
        classifications = np.load(cclass_path, allow_pickle=True)
        dt = time() - t0
        cls.logger.info(
            f"re-read compiled annotation index with {len(cdf)} original GTF feature combinations and {len(classifications)} in {dt:.3f} seconds"
        )
        ## Create a secondary Annotator which uses the non-overlapping combinations
        ## and the pre-classified annotations for the actual tagging
        cl = CompiledClassifier(cdf, classifications)
        gc = cls(cdf, lambda idx: cl.process(idx), is_compiled=True)
        return gc

    @classmethod
    def from_GTF(cls, gtf, df_cache=""):
        # load GTF the first time. Need to build compiled annotation
        t0 = time()
        df = load_GTF(gtf)
        dt = time() - t0
        cls.logger.info(f"loaded {len(df)} GTF records in {dt:.3f} seconds")
        if df_cache:
            df.to_csv(df_cache, sep="\t")

        ## Build NCLS with original GTF features
        cl = GTFClassifier(df)
        ga = cls(df, lambda idx: cl.process(idx))
        return ga

    @classmethod
    def from_uncompiled_df(cls, path):
        t0 = time()
        df = pd.read_csv(path, sep="\t")
        dt = time() - t0
        cls.logger.info(f"loaded {len(df)} tabular records in {dt:.3f} seconds")

        ## Build NCLS with original GTF features
        cl = GTFClassifier(df)
        ga = cls(df, lambda idx: cl.process(idx))
        return ga

    def sanity_check(self, df):
        for strand, nc in self.strand_map.items():
            starts, ends, idx = np.array(nc.intervals()).T
            self.logger.debug(
                f"checking {strand} with starts from {starts.min()}-{starts.max()} and idx from {idx.min()}-{idx.max()}"
            )

            assert idx.max() < len(df)

    def query_idx(self, chrom, start, end, strand):
        strand_key = (chrom, strand)
        if not strand_key in self.strand_map:
            return self.empty

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

    def compile(self, path=""):
        chroms = []
        strands = []
        cstarts = []
        cends = []
        cids = []
        cidx = []
        cid_lkup = {}

        t0 = time()
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
        dt = time() - t0
        self.logger.debug(f"decomposed original GTF annotation in {dt:.3f} seconds")

        # pre-classify all unique feature combinations
        classifications = []
        t0 = time()
        for n, idx in enumerate(cidx):
            res = self.processor(idx)
            classifications.append(res)

        ## Turn into np.array to turn lookup into super-fast array index operation
        classifications = np.array(classifications, dtype=object)
        dt = time() - t0
        self.logger.debug(
            f"pre-classified {n} feature combinations in {dt:.3f} seconds"
        )

        if path:
            cdf_path, cclass_path = CompiledClassifier.get_filenames(path)
            ## Store the compiled DataFrame and pre-classifications
            t0 = time()
            cdf.to_csv(cdf_path, sep="\t", index=False)
            np.save(cclass_path, classifications)
            dt = time() - t0
            self.logger.debug(f"stored compiled annotation index in {dt:.3f} seconds")

        ## Create a secondary Annotator which uses the non-overlapping combinations
        ## and the pre-classified annotations for the actual tagging
        cl = CompiledClassifier(cdf, classifications)
        gc = GenomeAnnotation(cdf, lambda idx: cl.process(idx), is_compiled=True)

        return gc

    def annotate_BAM(self, src, out, antisense=False, interval=5):
        import pysam

        as_strand = {"+": "-", "-": "+"}
        self.logger.info(
            f"beginning BAM annotation: {src} -> {out}. is_compiled={self.is_compiled}"
        )
        bam = pysam.AlignmentFile(src)
        out = pysam.AlignmentFile(out, "wbu", template=bam)
        t0 = time()
        T = interval
        for n, read in enumerate(bam.fetch(until_eof=True)):
            if not read.is_unmapped:
                chrom = bam.get_reference_name(read.tid)
                strand = "-" if read.is_reverse else "+"

                gf, gn, gs, gt = self.query_blocks(chrom, strand, read.get_blocks())
                if antisense:
                    gf_as, gn_as, gs_as, gt_as = self.query_blocks(
                        chrom, as_strand[strand], read.get_blocks()
                    )
                    gf += gf_as
                    gn += gn_as
                    gs += gs_as
                    gt += gt_as

                if len(gf):
                    read.tags += [
                        ("gF", ",".join(gf)),
                        ("gN", ",".join(gn)),
                        ("gS", ",".join(gs)),
                        ("gT", ",".join(gt)),
                    ]
                else:
                    read.tags += [("gF", "INTERGENIC")]

            out.write(read)
            dt = time() - t0
            if dt > T:
                self.logger.info(
                    f"processed {n} alignments in {dt:.2f} seconds ({n/dt:.2f} reads/second)"
                )
                T += interval

        self.logger.info(
            f"processed {n} alignments in {dt:.2f} seconds ({n/dt:.2f} reads/second)"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gtf",
        required=True,
        help="path to the original annotation (e.g. gencodev38.gtf.gz)",
    )
    parser.add_argument(
        "--tabular",
        default="",
        help="path to tabular version of the relevant features only (e.g. gencodev38.tsv)",
    )
    parser.add_argument(
        "--compiled",
        help="path to keep a compiled version of the GTF in (used as cache)",
    )
    parser.add_argument(
        "--use-compiled",
        default=False,
        action="store_true",
        help="enable using the compiled version of GenomeAnnotation (faster)",
    )
    parser.add_argument(
        "--antisense",
        default=False,
        action="store_true",
        help="enable annotating against the opposite strand (antisense to the alignment) as well",
    )
    parser.add_argument(
        "--dropseqtools",
        default=False,
        action="store_true",
        help="emulate dropseqtools behavior",
    )

    parser.add_argument("--bam-in", help="path for the input BAM to be tagged")
    parser.add_argument(
        "--bam-out",
        help="path for the tagged BAM output",
    )
    args = parser.parse_args()

    if CompiledClassifier.files_exist(args.compiled) and args.use_compiled:
        ga = GenomeAnnotation.from_compiled_index(args.compiled)

    elif args.tabular and os.access(args.tabular, os.R_OK):
        ga = GenomeAnnotation.from_uncompiled_df(args.tabular)

    else:
        ga = GenomeAnnotation.from_GTF(args.gtf, df_cache=args.tabular)

    # perform compilation if that's what we want
    if not ga.is_compiled and args.use_compiled:
        ga = ga.compile(args.compiled)

    ga.annotate_BAM(args.bam_in, args.bam_out, antisense=args.antisense)
