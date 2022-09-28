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
from collections import defaultdict, OrderedDict, deque, namedtuple
from dataclasses import dataclass
import typing
import cProfile

logging.basicConfig(level=logging.DEBUG)


def attr_to_dict(attr_str: str) -> dict:
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
    src: typing.Union[str, typing.TextIO],
    attributes: list = ["transcript_id", "transcript_type", "gene_name"],
    features: list = ["exon", "CDS", "UTR", "transcript"],
) -> pd.DataFrame:
    """
    Loads GTF-formatted annotation data (such as GENCODE) and converts into a tabular DataFrame. Coordinates will be converted
    to zero-based, half-open interval (C, python, BED,...).

    Args:
        src (filename or file-like object): Source of GTF annotation records. Can be gzip compressed (filename ending in .gz)
        attributes (list, optional): Which attributes to extract from the attribute column. Defaults to ["transcript_id", "transcript_type", "gene_name"].
        features (list, optional): which GTF features to care about. Defaults to ["exon", "CDS", "UTR", "transcript"].

    Returns:
        pandas.DataFrame: all GTF annotation data converted to a DataFrame with columns: chrom, feature, start (0-based), end, strand, [attributes]
    """
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
            # not a GTF/GFF formatted line
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


feat2flag = {
    "transcript": 1,
    "exon": 2,
    "UTR": 4,
    "CDS": 8,
}

# idea: encode as bitmask. Sense (low nibble) and antisense (high nibble) sense would both fit into one byte!
default_map = {
    # CDS UTR exon transcript-> gm, gf tag values
    0b1011: ("exon", "cds"),
    # This should not arise! Unless, isoforms-are merged
    0b1111: ("exon", "cds_utr"),
    0b0111: ("exon", "utr"),
    0b0011: ("exon", "nc"),
    0b1001: ("intron", "cds"),
    0b0101: ("intron", "utr"),
    0b0001: ("intron", "nc"),
    0b0000: ("intergenic", "n/a"),
}
default_lookup = np.array(
    [default_map.get(i, ("undefined", "n/a")) for i in np.arange(256)], dtype=object
)


@dataclass
class IsoformOverlap:
    gene_name: str = "no_name"
    transcript_type: str = "no_type"
    flags: int = 0


class GTFClassifier:
    """
    _summary_
    A simple classifier. Initizialized with a DataFrame of GTF features, the process() function
    will combine GTF features for each transcript to produce condensed annotations for the
    annotation bam tags.
    To that end, the presence/absence of CDS, UTR, exon, transcript features in the set of DataFrame
    indices passed into process() is compiled into a boolean vector. This vector is translated into
    values for gene-model and function values using a feature_map. No hierarchy is enforced at this
    level. All annotations are reported based on the aggregate isoform info retrieved via the feature
    ids.
    """

    def __init__(
        self,
        df,
        hierarchy=["CDS", "UTR", "non-coding", "intron"],
        flags_lookup=default_lookup,
    ):
        self.flags_lookup = flags_lookup
        self.hierarchy = hierarchy
        self.prio = {}
        for i, h in enumerate(hierarchy):
            self.prio[h] = i

        df["feature_flag"] = df["feature"].apply(lambda x: feat2flag[x])
        self.df_columns = df.columns
        self.df = df.values
        self.df_core = df[
            ["transcript_id", "transcript_type", "gene_name", "feature_flag"]
        ].values
        # print(self.df)

    def process(self, ids: typing.Iterable[int]) -> tuple:
        """
        Gather annotation records identified by DataFrame (array) indices
        by isoform: gene_name, gene_type and overlap-flags

        Args:
            ids (Iterable of ints): indices into the array-converted DataFrame of GTF records

        Returns:
            (gn, gm, gf, gt): a tuple with tag values encoding the overlapping GTF features
        """
        # print(f">>> process({ids})")
        isoform_overlaps = defaultdict(IsoformOverlap)
        # iterate over the overlapping GTF features group by transcript_id as we go
        i: int
        for i in ids:
            # print(f"df entry for id={i}: {self.df_core[i]}")
            (transcript_id, transcript_type, gene_name, flag) = self.df_core[i]
            info = isoform_overlaps[transcript_id]
            info.gene_name = gene_name
            info.transcript_type = transcript_type
            info.flags |= flag

        return overlaps_to_tags(isoform_overlaps)  # isoform_overlaps


def overlaps_to_tags(isoform_overlaps: dict, flags_lookup=default_lookup) -> tuple:
    tags = set()
    for transcript_id, info in isoform_overlaps.items():
        _gm, _gf = flags_lookup[info.flags]
        # print(f">> {transcript_id} => {info.flags:04b} => {_gm} {_gf}")
        tags.add((info.gene_name, _gm, _gf, info.transcript_type))

    gn = []
    gm = []
    gf = []
    gt = []
    for n, m, f, t in tags:
        gn.append(n)
        gm.append(m)
        gf.append(f)
        gt.append(t)

    return gn, gm, gf, gt


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
        cls.logger.info(f"loading compiled index from '{path}'")

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
        cls.logger.info(f"loading GTF from '{gtf}'")
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
        cls.logger.info(f"loading from uncompiled df '{path}'")
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
        self.logger.info(f"compiling to '{path}'...")

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
        self.logger.debug(
            f"decomposed original GTF annotation in {dt:.3f} seconds. Pre-classifying all disjunct combinations..."
        )

        # pre-classify all unique feature combinations
        classifications = []
        t0 = time()
        for n, idx in enumerate(cidx):
            res = self.processor(idx)
            # classifications.append(overlaps_to_tags(res))
            classifications.append(res)

        ## Turn into np.array to make lookup a super-fast array index operation
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

    def annotate_BAM(self, src, out, antisense=False, interval=5, repeat=5):
        import pysam

        as_strand = {"+": "-", "-": "+"}
        self.logger.info(
            f"beginning BAM annotation: {src} -> {out}. is_compiled={self.is_compiled}"
        )
        bam = pysam.AlignmentFile(src)
        out = pysam.AlignmentFile(out, "wbu", template=bam, threads=8)
        t0 = time()
        T = interval
        n: str = 0
        for read in bam.fetch(until_eof=True):
            for i in range(repeat):
                n += 1
                if not read.is_unmapped:
                    chrom = bam.get_reference_name(read.tid)
                    strand = "-" if read.is_reverse else "+"

                    # isoform_overlaps = self.query_blocks(
                    #     chrom, strand, read.get_blocks()
                    # )
                    # gn, gm, gf, gt = overlaps_to_tags(isoform_overlaps)
                    gn, gm, gf, gt = self.query_blocks(chrom, strand, read.get_blocks())
                    if antisense:
                        gf_as, gn_as, gm_as, gt_as = self.query_blocks(
                            chrom, as_strand[strand], read.get_blocks()
                        )
                        gf += gf_as
                        gn += gn_as
                        gm += gm_as
                        gt += gt_as

                    if len(gf):
                        read.set_tag("gF", ",".join(gf))
                        read.set_tag("gN", ",".join(gn))
                        read.set_tag("gM", ",".join(gm))
                        read.set_tag("gT", ",".join(gt))
                    else:
                        read.set_tag("gF", "INTERGENIC")

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


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gtf",
        required=True,
        help="path to the original annotation (e.g. gencodev38.gtf.gz)",
    )
    parser.add_argument(
        "--repeat",
        type=int,
        default=1,
        help="repeat each read n times (for benchmarking purposes)",
    )
    parser.add_argument(
        "--tabular",
        default="",
        help="path to tabular version of the relevant features only (e.g. gencodev38.tsv)",
    )
    parser.add_argument(
        "--compiled",
        default="",
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

    parser.add_argument("--bam-in", default="", help="path to the input BAM")
    parser.add_argument("--bam-out", default="", help="path for the tagged BAM output")
    args = parser.parse_args()

    return args


def main(args):
    logger = logging.getLogger("annotator")
    if CompiledClassifier.files_exist(args.compiled) and args.use_compiled:
        ga = GenomeAnnotation.from_compiled_index(args.compiled)

    elif args.tabular and os.access(args.tabular, os.R_OK):
        ga = GenomeAnnotation.from_uncompiled_df(args.tabular)

    else:
        ga = GenomeAnnotation.from_GTF(args.gtf, df_cache=args.tabular)
    # perform compilation if that's what we want
    if not ga.is_compiled and args.use_compiled:
        ga = ga.compile(args.compiled)

    if args.bam_in:
        ga.annotate_BAM(
            args.bam_in, args.bam_out, antisense=args.antisense, repeat=args.repeat
        )


def cmdline():
    args = parse_args()
    main(args)


if __name__ == "__main__":
    cProfile.run("cmdline()", "prof_stats")
