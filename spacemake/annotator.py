import pandas as pd
import numpy as np
import os
import ncls
import argparse
import re
import gzip
import logging
import ncls
from time import time
from collections import defaultdict
from dataclasses import dataclass
import typing
import cProfile
import pysam
import multiprocessing as mp

from spacemake.parallel import (
    join_with_empty_queues,
    order_results,
    ExceptionLogging,
    log_qerr,
)
import spacemake.util as util

## GTF I/O
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

    df = pd.DataFrame(
        data, columns=["chrom", "feature", "start", "end", "strand"] + attributes
    ).drop_duplicates()

    df["feature_flag"] = df["feature"].apply(lambda x: feat2bit[x])
    df["transcript_type"] = df["transcript_type"].apply(
        lambda t: abbreviations.get(t, t)
    )

    return df


## Efficient representation of GTF features as bit-mask
feat2bit = {
    "transcript": 1,
    "exon": 2,
    "UTR": 4,
    "CDS": 8,
}

# IDEA: plus (low nibble) and minus strand (high nibble) sense would both fit into one byte!
# extracting only specific-strand features:
#
#   f_plus = mask & 0b1111
#   f_minus = mask >> 4

default_map = {
    # CDS UTR exon transcript-> gm, gf tag values
    0b1011: "C",  # CDS exon
    # This should not arise! Unless, isoforms-are merged
    0b1111: "CU",  # CDS and UTR exon
    0b0111: "U",  # UTR exon
    0b0011: "N",  # exon of a non-coding transcript
    0b1001: "I",  # intron
    0b0101: "I",
    0b0001: "I",
    0b0000: "-",  # intergenic
    ## TODO: Antisense is shifted left by 4 bits and uses lower-case letters
    0b10110000: "c",  # CDS exon (antisense)
    # This should not arise! Unless, isoforms-are merged
    0b11110000: "cu",  # CDS and UTR exon (antisense)
    0b01110000: "u",  # UTR exon (antisense)
    0b00110000: "n",  # exon of a non-coding transcript (antisense)
    0b10010000: "i",  # intron (antisense)
    0b01010000: "i",  # (antisense)
    0b00010000: "i",  # (antisense)
}
default_lookup = np.array(
    [default_map.get(i, "?") for i in np.arange(256)], dtype=object
)
abbreviations = {
    "protein_coding": "c",
    "processed_transcript": "t",
    "processed_pseudogene": "p",
    "lncRNA": "l",
    "nonsense_mediated_decay": "n",
    "retained_intron": "i",
}


## Space-efficient container of aggregating overlap information (isoform resolution)
## TODO: find a better, more compact/faster intermediate layer
@dataclass
class IsoformOverlap:
    gene_name: str = "no_name"
    transcript_type: str = "no_type"
    flags: int = 0


def overlaps_to_tags(isoform_overlaps: dict, flags_lookup=default_lookup) -> tuple:
    """_summary_
    Take a dictionary of transcript_id -> IsoformOverlap as input and populate annotation tuples

    Args:
        isoform_overlaps (dict): _description_
        flags_lookup (_type_, optional): _description_. Defaults to default_lookup.

    Returns:
        tuple: gn, gf, gt (each is a list of strings n: gene names f: function t: transcript type)
    """
    tags = set()
    for transcript_id, info in isoform_overlaps.items():
        _gf = flags_lookup[info.flags]
        # print(f">> {transcript_id} => {info.flags:04b} => {_gm} {_gf}")
        tags.add((info.gene_name, _gf, info.transcript_type))

    gn = []
    gf = []
    gt = []
    for n, f, t in tags:
        gn.append(n)
        gf.append(f)
        gt.append(t)

    return gn, gf, gt


# First implementation: uses the GTF features directly. Not very optimized. Mainly used to
# build compiled annotation (see below).
# A Classifier instance needs to be passed to GenomeAnnotation as a parameter to __init__ (see below)
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

        self.df_columns = df.columns
        self.df = df.values
        self.df_core = df[
            ["transcript_id", "transcript_type", "gene_name", "feature_flag"]
        ].values

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


## Second iteration. Here, we load already pre-compiled, unique combinations of original GTF
# features as they occur in the genome. The beauty is that each of these "compiled" features now
# has zero overlap with other compiled features and each of them already has the gn, gf, gt values
# pre-evaluated ("compiled"). The only issue arises if a read overlaps multiple of these pre-compiled
# "meta"-features. For this we need the joiner() function which aims to combine two sets of gn, gf, gt
# values in a meaningful way. This is slow, but rare enough to leave a substantial speed boost
# on average, when comparing to using raw GTF features directly.
# of note, this still plugs into the same GenomeAnnotation facility (see below). It just changes the
# meaning and number of the underlying features (from individual GTF records to combinatorial overlaps
# these GTF records as they occur in the genome), and it also changes the processing because we no
# longer deal with DataFrame entries but already with gn, gf, gt tuples.
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
        gt = []
        for i in tuple(cidx):
            cid = self.cid_table[i]
            for ann in zip(*self.classifications[cid]):
                if not ann in seen:
                    seen.add(ann)
                    gf.append(ann[0])
                    gn.append(ann[1])
                    gt.append(ann[2])

        return gf, gn, gt

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
    TODO: directly hook into NCLSIterator functionality to get rid of this overhead
    """
    # return frozenset((o[2] for o in nc.find_overlap(x0, x1)))
    return [o[2] for o in nc.find_overlap(x0, x1)]


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
    last_key = frozenset(query(nc, last_pos, last_pos + 1))
    for bp in breakpoints[1:]:
        # depending on wether this is another start or
        # end position, the transition point can be off by one nt
        # in either direction. Most robust fix was to scan all 3
        # options until the key changes (while ensuring that end > start)
        for x in [-1, 0, +1]:
            if bp + x <= last_pos:
                continue

            key = frozenset(query(nc, bp + x, bp + x + 1))
            if key != last_key:
                if len(last_key):
                    # ensure we do not yield the empty set
                    yield last_pos, bp + x, last_key
                last_key = key
                last_pos = bp + x
                break


class GenomeAnnotation:
    """
    NCLS-based lookup of genome annotation. This is the top-level interface for the query of genomic intervals.
    GenomeAnnotation is actually a pretty thin wrapper around NCLS for fast 1d-interval queries. The raw
    query results are just sets of numerical ids. Interpretation of these sets is deferred to a "processor", an
    instance of either GTFClassifier or CompiledClassifier. These take on the task of processing the raw ids
    into tuples of gn, gf, gt.
    Because these two layers are tightly integrated, there are convenience classmethods to construct a fully
    functional GenomeAnnotation instance for the most common scenarios (from GTF, from compiled data, ...).
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
        self.strand_keys = []
        self.strand_map = {}
        self.empty = []  # frozenset([])

        t0 = time()
        for chrom, strand in df[["chrom", "strand"]].drop_duplicates().values:
            strand_key = chrom + strand

            d = df.query(f"chrom == '{chrom}' and strand == '{strand}'")
            nested_list = ncls.NCLS(d["start"], d["end"], d.index)
            self.strand_map[strand_key] = nested_list
            self.strand_keys.append(strand_key)

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
            f"loaded compiled annotation index with {len(cdf)} original GTF feature combinations and {len(classifications)} in {dt:.3f} seconds"
        )
        ## Create a classifier which uses the pre-compiled, non-overlapping combinations
        ## and corresponding pre-computed annotation tags
        cl = CompiledClassifier(cdf, classifications)
        ga = cls(cdf, lambda idx: cl.process(idx), is_compiled=True)
        return ga

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
        strand_key = chrom + strand
        nested_list = self.strand_map.get(strand_key, [])
        if nested_list:
            return query(nested_list, start, end)
        else:
            return self.empty

    def query_idx_blocks(self, chrom, strand, blocks):
        idx = []
        for start, end in blocks:
            idx.extend(self.query_idx(chrom, start, end, strand))

        return idx

    def query(self, chrom, start, end, strand):
        idx = self.query_idx(chrom, start, end, strand)
        return self.processor(idx)

    def query_blocks(self, chrom, strand, blocks):
        idx = self.query_idx_blocks(chrom, strand, blocks)
        return self.processor(idx)

    def compile(self, path=""):
        if self.is_compiled:
            raise ValueError(
                "attempting to compile an already compiled annotation. Why?? :-( "
            )

        self.logger.info(f"compiling to '{path}'...")
        if path:
            path = util.ensure_path(path + "/")

        chroms = []
        strands = []
        cstarts = []
        cends = []
        cids = []
        cidx = []
        cid_lkup = {}

        t0 = time()
        for strand_key, nc in self.strand_map.items():
            chrom = strand_key[:-1]
            strand = strand_key[-1]
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
            # classifications.append(to_tags(res))
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

    def get_annotation_tags(self, chrom: str, strand: str, blocks, antisense=False):
        def collapse_tag_values(values):
            if len(set(values)) == 1:
                return [
                    values[0],
                ]
            else:
                return values

        gn, gf, gt = self.query_blocks(chrom, strand, blocks)
        if antisense:
            gn_as, gf_as, gt_as = self.query_blocks(
                chrom, as_strand[strand], read.get_blocks()
            )
            gn += gn_as
            gf += gf_as
            gt += gt_as

        return collapse_tag_values(gn), gf, collapse_tag_values(gt)


def merge_annotations_from_queue_with_BAM(
    bam_in, bam_out, bam_mode, Qres, Qerr, abort_flag
):
    """_summary_
    The output generating sub-process. Iterate over the input BAM and gather the annotation tags generated by
    parallel worker sub-processes simultaneously. Then assign the tags to the BAM records and write them to the
    output BAM.

    Args:
        bam_in (str): file name of input BAM file
        bam_out (str): fila name of output BAM file
        bam_mode (str): mode for output BAM. Can be used to set compression level (e.g. bu or b6)
        Qres (mp.Queue): source for chunks of annotation data
        Qerr (mp.Queue): in case of error, dump exception information onto this queue so that it can be
            extracted from the main process
        abort_flag (mp.Value): inter-process flag to indicate abort conditions (i.e sth went wrong elsewhere)
    """
    with ExceptionLogging(
        "spacemake.annotator.order_results", Qerr=Qerr, exc_flag=abort_flag
    ) as el:
        logger = el.logger
        logger.info(f"beginning BAM annotation: {bam_in} -> {bam_out}.")
        bam = pysam.AlignmentFile(bam_in, threads=2)
        out = pysam.AlignmentFile(bam_out, f"w{bam_mode}", template=bam, threads=8)

        for read, annotation in zip(
            bam.fetch(until_eof=True), order_results(Qres, abort_flag, logger=el.logger)
        ):
            # set BAM tags and write
            # (n_read, qname, gn_val, gf_val, gt_val) = annotation
            (gn_val, gf_val, gt_val) = annotation
            # assert qname == read.query_name
            read.set_tag("gn", gn_val)
            read.set_tag("gf", gf_val)
            read.set_tag("gt", gt_val)
            read.set_tag("XF", None)
            read.set_tag("gs", None)
            out.write(read)


def annotation_chunks_from_BAM(
    bam_in, compiled_annotation, rr_n=0, rr_i=0, chunk_size=100
):
    """_summary_
    The main work of the parallel workers is actually implemented in this generator.
    It can also be run in a non-parallel, single-process fashion as well (set rr_n=0 and rr_i=0).

    Args:
        bam_in (_type_): _description_
        compiled_annotation (_type_): _description_
        rr_n (int, optional): _description_. Defaults to 0.
        rr_i (int, optional): _description_. Defaults to 0.
        chunk_size (int, optional): _description_. Defaults to 100.

    Returns:
        None

    Yields:
        n_chunk, chunk: list of up to chunk_size tuples of (gn_val, gf_val, gt_val) strings. n_chunk is consistent
        across parallel workers and increments linearly. Allows to sort chunks for correct output ordering and merging
        with the original BAM records from the input.
    """
    import pysam

    # as_strand = {"+": "-", "-": "+"}
    logger = logging.getLogger("interval_chunks_from_BAM")
    ga = GenomeAnnotation.from_compiled_index(compiled_annotation)

    t0 = time()
    T = 5
    tid_cache = {}

    bam = pysam.AlignmentFile(bam_in, threads=2)

    def get_reference_name(tid):
        if not tid in tid_cache:
            tid_cache[tid] = bam.get_reference_name(tid)
        return tid_cache[tid]

    n_read = 0
    n_chunk = 0
    n_in_chunk = 0
    chunk = []
    for read in bam.fetch(until_eof=True):
        if rr_n:
            rr = n_chunk % rr_n
        else:
            rr = 0
        if rr == rr_i:
            gn_val = None
            gf_val = "-"
            gt_val = None
            if not read.is_unmapped:
                chrom = get_reference_name(read.tid)
                strand = "-" if read.is_reverse else "+"
                gn, gf, gt = ga.get_annotation_tags(chrom, strand, read.get_blocks())
                if len(gf):
                    gn_val = ",".join(gn)
                    gf_val = ",".join(gf)
                    gt_val = ",".join(gt)

            chunk.append((gn_val, gf_val, gt_val))

        n_read += 1
        n_in_chunk += 1
        if n_in_chunk >= chunk_size:
            if rr == rr_i:
                yield n_chunk, chunk

            chunk = []
            n_in_chunk = 0
            n_chunk += 1

        dt = time() - t0
        if dt > T:
            logger.info(
                f"processed {n_read} alignments in {dt:.2f} seconds ({n_read/dt:.2f} reads/second)"
            )
            T += 5

    if chunk:
        yield n_chunk, chunk

    logger.info(
        f"processed {n_read} alignments in {dt:.2f} seconds ({n_read/dt:.2f} reads/second)"
    )


def round_robin_worker(
    bam_in, compiled_annotation, Qres, Qerr, rr_n, rr_i, abort_flag, chunk_size
):
    """_summary_
    Worker sub-process. Iterate over the BAM file independently, skip everything that is not in our chunk
    (every rr_i'th chunk) if chunk_size BAM records. Annotate our chunk using compiled GenomeAnnotation
    and place complete gn, gf, gt string values on a results Queue.

    Args:
        bam_in (str): file name of input BAM
        compiled_annotation (str): path to directory with compiled annotation information
        Qres (mp.Quere): Queue to place chunks of complete annotation tags in
        Qerr (mp.Queue): see merge_annotations_from_queue_with_BAM
        rr_n (int): how many round-robin workers are there in total
        rr_i (int): which round-robin worker are we? 0..rr_n-1 . Determines which
            chunks this worker feels responsible for.
        abort_flag (mp.Value): see merge_annotations_from_queue_with_BAM
        chunk_size (int): number of consecutive BAM records to process
    """
    with ExceptionLogging(f"worker_{rr_i}", Qerr=Qerr, exc_flag=abort_flag) as el:
        el.logger.debug(
            f"round_robin_worker {rr_i}/{rr_n} starting up with Qres={Qres}"
        )

        for n_chunk, data in annotation_chunks_from_BAM(
            bam_in, compiled_annotation, rr_n=rr_n, rr_i=rr_i, chunk_size=chunk_size
        ):
            Qres.put((n_chunk, data))


def annotate_BAM_parallel(args):
    """_summary_
    Main function of the 'annotate' command. Create the plumbing for parallel worker processes and a single
    collector/writer process.

    Args:
        args (namespace): the command-line arguments reported from the parser
    """
    Qres = (
        mp.Queue()
    )  # annotation tags are placed by the worker processes in chunks onto this queue
    Qerr = mp.Queue()  # child-processes can report errors back to the main process here

    # Proxy objects to allow workers to report statistics about the run
    abort_flag = mp.Value("b")
    abort_flag.value = False

    with ExceptionLogging("annotate_BAM_parallel", exc_flag=abort_flag) as el:
        workers = []
        for i in range(args.parallel):
            w = mp.Process(
                target=round_robin_worker,
                name=f"worker_{i}",
                args=(
                    args.bam_in,
                    args.compiled,
                    Qres,
                    Qerr,
                    args.parallel,
                    i,
                    abort_flag,
                    args.chunk_size,
                ),
            )
            w.start()
            workers.append(w)

        el.logger.info("Started workers")
        collector = mp.Process(
            target=merge_annotations_from_queue_with_BAM,
            name="output",
            args=(args.bam_in, args.bam_out, args.bam_mode, Qres, Qerr, abort_flag),
        )
        collector.start()
        el.logger.info("Started collector")

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


def annotate_BAM_linear(bam, ga, out, repeat=1, interval=5):
    # import pysam

    # as_strand = {"+": "-", "-": "+"}
    logger = logging.getLogger("chunks_from_BAM")

    t0 = time()
    T = interval
    n: str = 0

    tid_cache = {}

    def get_reference_name(tid):
        if not tid in tid_cache:
            tid_cache[tid] = bam.get_reference_name(tid)
        return tid_cache[tid]

    for read in bam.fetch(until_eof=True):
        n += 1
        gn_val = None
        gf_val = "-"
        gt_val = None
        if not read.is_unmapped:
            chrom = get_reference_name(read.tid)
            strand = "-" if read.is_reverse else "+"
            gn, gf, gt = ga.get_annotation_tags(chrom, strand, read.get_blocks())
            if len(gf):
                gn_val = ",".join(gn)
                gf_val = ",".join(gf)
                gt_val = ",".join(gt)

        read.set_tag("gn", gn_val)
        read.set_tag("gf", gf_val)
        read.set_tag("gt", gt_val)
        read.set_tag("XF", None)
        read.set_tag("gs", None)
        out.write(read)

        dt = time() - t0
        if dt > T:
            logger.info(
                f"processed {n} alignments in {dt:.2f} seconds ({n/dt:.2f} reads/second)"
            )
            T += interval

    logger.info(
        f"processed {n} alignments in {dt:.2f} seconds ({n/dt:.2f} reads/second)"
    )


def build_compiled_annotation(args):
    logger = logging.getLogger("spacemake.annotator.build_compiled_annotation")
    if CompiledClassifier.files_exist(args.compiled):
        logger.warning(
            "already found a compiled annotation. use --force-overwrite to overwrite"
        )
        # ga = GenomeAnnotation.from_compiled_index(args.compiled)
        if not args.force_overwrite:
            return

    if args.tabular and os.access(args.tabular, os.R_OK):
        ga = GenomeAnnotation.from_uncompiled_df(args.tabular)
    else:
        ga = GenomeAnnotation.from_GTF(args.gtf, df_cache=args.tabular)

    ga = ga.compile(args.compiled)
    return ga


def parse_args():

    parser = argparse.ArgumentParser()

    def usage(args):
        parser.print_help()

    parser.set_defaults(func=usage)

    subparsers = parser.add_subparsers()

    build_parser = subparsers.add_parser("build")
    build_parser.set_defaults(func=build_compiled_annotation)

    tag_parser = subparsers.add_parser("tag")
    tag_parser.set_defaults(func=annotate_BAM_parallel)

    build_parser.add_argument(
        "--gtf",
        default=None,
        required=True,
        help="path to the original annotation (e.g. gencodev38.gtf.gz)",
    )
    build_parser.add_argument(
        "--compiled",
        default=None,
        help="path to a directoy in which a compiled version of the GTF is stored",
    )
    build_parser.add_argument(
        "--tabular",
        default="",
        help="path to a cache of the tabular version of the relevant GTF features (optional)",
    )
    build_parser.add_argument(
        "--force-overwrite",
        default=False,
        action="store_true",
        help="re-compile GTF and overwrite the pre-existing compiled annotation",
    )

    tag_parser.add_argument(
        "--compiled",
        default=None,
        help="path to a directoy in which a compiled version of the GTF is stored",
    )
    tag_parser.add_argument("--bam-in", default="", help="path to the input BAM")
    tag_parser.add_argument(
        "--bam-out", default="", help="path for the tagged BAM output"
    )
    tag_parser.add_argument(
        "--bam-mode", default="b", help="mode of the output BAM file (default=b)"
    )
    tag_parser.add_argument(
        "--parallel",
        type=int,
        default=4,
        help="how many parallel annotation processes to create (default=4)",
    )
    tag_parser.add_argument(
        "--chunk-size",
        type=int,
        default=10000,
        help="how many BAM-records form a chunk for parallel processing (default=10000)",
    )
    tag_parser.add_argument(
        "--antisense",
        default=False,
        action="store_true",
        help="enable annotating against the opposite strand (antisense to the alignment) as well",
    )

    return parser.parse_args()


# def main(args):
#     logger = logging.getLogger("annotator")
#     if not args.bam_in:
#         return

#     logger.info(f"beginning BAM annotation: {args.bam_in} -> {args.bam_out}.")
#     bam = pysam.AlignmentFile(args.bam_in, threads=2)
#     out = pysam.AlignmentFile(args.bam_out, "wb", template=bam, threads=8)

#     for chunk in interval_chunks_from_BAM(args.bam_in, args.compiled, 0, 0):
#         pass
#         # for data in chunk:
#         #     print(data)

# annotate_BAM_linear(bam, ga, out, repeat=args.repeat)

# for read_chunk, block_chunk in chunks_from_BAM(bam):
#     for read, aln in zip(read_chunk, block_chunk):
#         read.set_tag("XF", None)
#         read.set_tag("gs", None)
#         if aln:
#             gn, gf, gt = ga.get_annotation_tags(*aln)
#             if len(gf):
#                 read.set_tag("gn", ",".join(gn))
#                 read.set_tag("gf", ",".join(gf))
#                 read.set_tag("gt", ",".join(gt))
#             else:
#                 read.set_tag("gn", None)
#                 read.set_tag("gf", "-")
#                 read.set_tag("gt", None)

#         out.write(read)


def cmdline():
    args = parse_args()
    args.func(args)
    # main(args)
    # annotate_BAM_parallel(args)

    # if args.parallel > 1:
    #     annotate_BAM_parallel(args)
    # else:
    #     annotate_BAM_linear(args)


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    cmdline()
    # cProfile.run("cmdline()", "prof_stats")
