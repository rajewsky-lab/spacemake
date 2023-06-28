import pandas as pd
import numpy as np
import os
import ncls
import re
import gzip
import logging
from time import time
from collections import defaultdict
from dataclasses import dataclass
from ctypes import c_int
import typing
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

    if not "gene_name" in d:
        d["gene_name"] = d.get("gene_id", d.get("transcript_id", "na"))

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

    if type(src) is str:
        if src.endswith(".gz"):
            src = gzip.open(src, "rt")
        else:
            src = open(src, "rt")

    fset = set(features)
    data = []
    for line in src:
        if line.startswith('#'):
            continue

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

    def bitflags(columns):
        strand, feature = columns
        flag = feat2bit[feature]
        if strand == "-":
            # features on minus strand use the high-nibble of the bit-flags
            flag = flag << 4

        return flag

    df["feature_flag"] = df[["strand", "feature"]].apply(bitflags, axis=1)
    df["transcript_type"] = df["transcript_type"].apply(
        lambda t: abbreviations.get(t, t)
    )
    return df

def make_lookup(d):
    return np.array(
        [d.get(i, "?") for i in np.arange(256)], dtype=object
    )

default_lookup = make_lookup({
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
    0b10110000: "c",  # CDS exon
    0b11110000: "cu",  # CDS and UTR exon
    0b01110000: "u",  # UTR exon
    0b00110000: "n",  # exon of a non-coding transcript
    0b10010000: "i",  # intron
    0b01010000: "i",
    0b00010000: "i",
    0b00000000: "-",  # intergenic
})
default_strand_translation = str.maketrans("CUNIcuni", "cuniCUNI", "")

# short-hands for common transcript types
abbreviations = {
    "protein_coding": "C",
    "processed_transcript": "N",
    "processed_pseudogene": "P",
    "transcribed_unprocessed_pseudogene": "P",
    "transcribed_processed_pseudogene": "P",
    "transcribed_unitary_pseudogene": "P",
    "lncRNA": "L",
    "lincRNA": "I",
    "nonsense_mediated_decay": "D",
    "retained_intron": "R",
    "miRNA": "M",
    "tRNA": "t",
    "snoRNA": "s",
    "rRNA": "r",
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
        
        tuple: (gn, gf, gt)
            a tuple with tag values encoding the overlapping
            GTF features. Each is a list of strings :

                gn: gene names 
                gf: function
                gt: transcript type

                gf uses upper case letters for features on the '+' strand and
                lower case letters for features on the '-' strand.

                For a query against the minus strand, upper and lower case
                need to be swapped in order to get features that align in the 
                "sense" direction in upper case and "antisense" features in lower-case.

                For a strand-agnostic query, .upper() is called. 
                These case-mangling operations are carried out in get_annotation_tags()
        
    """
    # PROBLEM:
    # order of tag values is given by order of annotation records, not anything
    # reproducible, such as alphabetical order. This is makes testing and counting
    # un-necessarily hard. Let's try to collect by gene and sort gf,gt on gf for same gene

    # by_gn = defaultdict(set)
    # for transcript_id, info in isoform_overlaps.items():
    #     _gf = default_lookup[info.flags]

    #     # print(f">> {transcript_id} => {info.flags:04b} => {_gm} {_gf}")
    #     by_gn[info.gene_name].add((_gf, info.transcript_type))

    # gn = []
    # gf = []
    # gt = []
    # for n, tags in by_gn.items():
    #     for f,t in sorted(tags):
    #         gn.append(n)
    #         gf.append(f)
    #         gt.append(t)

    tags = set()
    for transcript_id, info in isoform_overlaps.items():
        _gf = default_lookup[info.flags]

        # print(f">> {transcript_id} => {info.flags:04b} => {_gm} {_gf}")
        tags.add((info.gene_name, _gf, info.transcript_type))

    gn = []
    gf = []
    gt = []
    for n, f, t in tags:
        gn.append(n)
        gf.append(f)
        gt.append(t)

    # print(f"overlaps_to_tags {gn} {gf} {gt}")
    return (gn, gf, gt)


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

    def process(self, ids: typing.Iterable[c_int]) -> dict:
        """
        Gather annotation records identified by DataFrame (array) indices
        by isoform: gene_name, gene_type and overlap-flags

        Args:
            ids (Iterable of ints): indices into the array-converted DataFrame of GTF records

        Returns:
            output of overlaps_to_tags
        """
        # print(f">>> process({ids})")
        isoform_overlaps = defaultdict(IsoformOverlap)
        # iterate over the overlapping GTF features group by transcript_id as we go
        i: c_int
        for i in ids:
            (transcript_id, transcript_type, gene_name, flag) = self.df_core[i]
            info = isoform_overlaps[transcript_id]
            info.gene_name = gene_name
            info.transcript_type = transcript_type
            info.flags |= flag 
            # print(f"df entry for id={i}: {self.df_core[i]} flags={info.flags:08b}")

        return overlaps_to_tags(isoform_overlaps)  # isoform_overlaps

    # def preprocess(self, ids: typing.Iterable[c_int]) -> tuple:
    #     """
    #     pre-process a list of features into full annotation for either strand, 
    #     as well as intermediate data objects that allow correct merging.
    #     These results will be stored when annotation is "compiled" in order to offer
    #     much faster lookups.

    #     Args:
    #         idx (_type_): _description_
    #     """

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
    # TODO: ensure deterministic ordering of joined annotation tags
    def joiner(self, cidx):
        # we want to report each combination only once
        # Note: This code path is executed rarely

        seen = set()
        gn = []
        gf = []
        gt = []
        for idx in tuple(cidx):
            ann_tags = self.classifications[self.cid_table[idx]]
            for ann in zip(*ann_tags):
                if not ann in seen:
                    seen.add(ann)
                    gn.append(ann[0])
                    gf.append(ann[1])
                    gt.append(ann[2])

        return (gn, gf, gt)

    def process(self, cidx, check_contained=False):
        # fast path
        gene_contains = defaultdict(int)

        if check_contained:
            cidx, start_included, end_included = cidx
            for idx, si, ei in zip(cidx, start_included, end_included):
                c = self.classifications[self.cid_table[idx]]
                for gene in c[0]:
                    gene_contains[gene] |= si
                    gene_contains[gene] |= 2*ei

        if len(cidx) == 1:
            # fast path
            idx = cidx.pop()
            ann = self.classifications[self.cid_table[idx]]
        else:
            ann = self.joiner(cidx)

        if check_contained:
            for gene, fcontains in gene_contains.items():
                if fcontains < 3:
                    ann[0].append(gene)
                    ann[1].append('-')
                    ann[2].append('-')

        return ann


## Helper functions for working with NCLS
def query_ncls(nc, x0, x1, check_contained=False):
    """
    report only the target_ids from the overlapping
    (start, end, id) tuples returned by NCLS.find_overlap().
    format is frozenset bc that can be hashed and used as a key.
    TODO: directly hook into NCLSIterator functionality to get rid of this overhead
    """
    
    if check_contained:
        keys = []
        start_contained = []
        end_contained = []
        for start, end, key in nc.find_overlap(x0, x1):
            start_contained.append(x0 >= start)
            end_contained.append(x1 <= end)
            keys.append(key)

        # print(f"query ({x0} - {x1}) strand ={strand} -> {res}")
        return keys, start_contained, end_contained
    
    else:
        res = [o[2] for o in nc.find_overlap(x0, x1)]
        return res

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
    logger = logging.getLogger("spacemake.annotator.compile")
    logger.debug(
        f"compiling strand into non-overlapping and pre-classified annotations"
    )

    starts, ends, ids = np.array(nc.intervals()).T
    logger.debug(f"retrieved {len(starts)} intervals from ncls")

    breakpoints = np.array(sorted(set(starts) | set(ends)))
    logger.debug(f"identified {len(breakpoints)} breakpoints")

    if not len(starts):
        return

    # since we sorted, the first must be a start position
    last_pos = breakpoints[0]
    last_key = frozenset(query_ncls(nc, last_pos, last_pos + 1))
    for bp in breakpoints[1:]:
        # depending on wether this is another start or
        # end position, the transition point can be off by one nt
        # in either direction. Most robust fix was to scan all 3
        # options until the key changes (while ensuring that end > start)
        for x in [-1, 0, +1]:
            if bp + x <= last_pos:
                continue

            key = frozenset(query_ncls(nc, bp + x, bp + x + 1))
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

    logger = logging.getLogger("spacemake.annotator.GenomeAnnotation")

    def __init__(self, df, classifier, is_compiled=False):
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
        self.classifier = classifier

        # find all unique combinations of chrom + strand
        # strands = df["chrom"].drop_duplicates().values
        self.chroms = []
        self.chrom_ncls_map = {}
        self.empty = []
        # df['sign'] = df['strand'].apply(lambda x: -1 if x == '-' else +1)

        t0 = time()
        for chrom in df["chrom"].drop_duplicates().values:
            d = df.query(f"chrom == '{chrom}'")
            nested_list = ncls.NCLS(d["start"], d["end"], d.index)
            self.logger.debug(
                f"constructed nested_list for {chrom} with {len(d)} features"
            )

            self.chrom_ncls_map[chrom] = nested_list
            self.chroms.append(chrom)

        dt = time() - t0
        self.logger.debug(
            f"constructed nested lists of {len(df)} features on {len(self.chroms)} strands in {dt:.3f}s"
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
        cls.logger.debug(
            f"loaded compiled annotation index with {len(cdf)} original GTF feature combinations and {len(classifications)} pre-compiled classifications in {dt:.3f} seconds"
        )
        ## Create a classifier which uses the pre-compiled, non-overlapping combinations
        ## and corresponding pre-computed annotation tags
        ga = cls(cdf, CompiledClassifier(cdf, classifications), is_compiled=True)
        return ga

    @classmethod
    def from_GTF(cls, gtf, df_cache=""):
        cls.logger.info(f"loading GTF from '{gtf}'")
        # load GTF the first time. Need to build compiled annotation
        t0 = time()
        df = load_GTF(gtf)
        dt = time() - t0
        cls.logger.debug(f"loaded {len(df)} GTF records in {dt:.3f} seconds")
        if df_cache:
            df.to_csv(df_cache, sep="\t")

        ## Build NCLS with original GTF features
        ga = cls(df, GTFClassifier(df))
        return ga

    @classmethod
    def from_uncompiled_df(cls, path):
        cls.logger.info(f"loading from uncompiled df '{path}'")
        t0 = time()
        df = pd.read_csv(path, sep="\t")
        dt = time() - t0
        cls.logger.debug(f"loaded {len(df)} tabular records in {dt:.3f} seconds")

        ## Build NCLS with original GTF features
        ga = cls(df, GTFClassifier(df))
        return ga

    # def sanity_check(self, df):
    #     for chrom, nc in self.chrom_ncls_map.items():
    #         starts, ends, idx = np.array(nc.intervals()).T
    #         self.logger.debug(
    #             f"checking {chrom} with starts from {starts.min()}-{starts.max()} and idx from {idx.min()}-{idx.max()}"
    #         )

    #         assert idx.max() < len(df)

    def query_idx(self, chrom, start, end, **kw):
        nested_list = self.chrom_ncls_map.get(chrom, [])
        # print(f"nested list for {chrom} -> {nested_list}")
        if nested_list:
            return query_ncls(nested_list, start, end, **kw)
        else:
            return self.empty

    def query_idx_blocks(self, chrom, blocks, **kw):
        idx = []
        for start, end in blocks:
            idx.extend(self.query_idx(chrom, start, end, **kw))

        return idx

    # TODO: do we ever need query_idx_blocks separate from query?
    def query(self, chrom, start, end, **kw):
        idx = self.query_idx(chrom, start, end, **kw)
        return self.classifier.process(idx, **kw)

    def query_blocks(self, chrom, blocks, **kw):
        idx = self.query_idx_blocks(chrom, blocks, **kw)
        return self.classifier.process(idx, **kw)

    def compile(self, path=""):
        if self.is_compiled:
            raise ValueError(
                "attempting to compile an already compiled annotation. Why?? %-/ "
            )

        self.logger.info(f"compiling to '{path}'...")
        if path:
            path = util.ensure_path(path + "/")

        chroms = []
        # strands = []
        cstarts = []
        cends = []
        cids = []
        cidx = []
        cid_lkup = {}

        t0 = time()
        for chrom, nc in self.chrom_ncls_map.items():
            # chrom = strand_key[:-1]
            # strand = strand_key[-1]
            self.logger.debug(f"decomposing {chrom}")
            # decompose will report combinations of real indices into the df
            # for both strands + and - . The task of properly assigning these is
            # shifted to process(), merge() and later queries.
            # In consequence a single compiled annotation element can include features 
            # from both strands! 
            for start, end, idx in decompose(nc):
                # print(f"start={start} end={end} idx={idx}")
                chroms.append(chrom)
                # strands.append(strand)
                cstarts.append(start)
                cends.append(end)
                if idx not in cid_lkup:
                    cid_lkup[idx] = len(cid_lkup)
                    cidx.append(idx)

                cids.append(cid_lkup[idx])

        self.logger.debug("done")
        cdf = pd.DataFrame(
            dict(chrom=chroms, start=cstarts, end=cends, cid=cids) # strand=strands, 
        )
        dt = time() - t0
        self.logger.debug(
            f"decomposed original GTF annotation in {dt:.3f} seconds. Pre-classifying all disjunct combinations..."
        )

        # pre-classify all unique feature combinations
        classifications = []
        t0 = time()
        for n, idx in enumerate(cidx):
            res = postprocess_tags(*self.classifier.process(idx)) # <- ((gn, gf, gt), isoform_overlaps dict)
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
        gc = GenomeAnnotation(cdf, CompiledClassifier(cdf, classifications), is_compiled=True)

        return gc

    def get_annotation_tags(self, chrom: str, strand: str, blocks):
        def collapse_tag_values(values):
            if len(set(values)) == 1:
                return [
                    values[0],
                ]
            else:
                return values

        gn, gf, gt = self.query_blocks(chrom, blocks, check_contained=False)  # True is experimental and breaks stuff
        if len(gf):
            gn = ",".join(collapse_tag_values(gn))
            gf = ",".join(gf)
            gt = ",".join(collapse_tag_values(gt))
        else:
            gn = "-"
            gf = "-"
            gt = "-"

        if strand == '-':
            gf = gf.translate(default_strand_translation)
        elif strand == '*':
            gf = gf.upper()

        return gn, gf, gt
    

def postprocess_tags(gn, gf, gt):
    """
    post-process multiple functional annotations for a single gene,
    as can arise if we have multiple isoforms in a region.
    Assuming that this is called during compilation, the query region
    is a minimal combination of uniquely overlapping features. 
    Since we postprocess annotation for each gene separately, some 
    combinations can then ONLY arise by isoform differences. This is
    not the same scenario as overlapping distinct functional elements of
    all isoforms. We want:
        'C,I' -> "overlaps exonic regions *and* intronic regions in all known isoforms"
        'C|I' -> "overlaps a region that is coding in some and intronic in other isoforms"
        "C,C|I,I" -> "overlaps a region that is coding in all isoforms, a region that is coding in
        some and intronic in others, and a region that is intronic in all isoforms"
    """
    # print("input", gn, gf, gt)
    from collections import defaultdict
    by_gn = defaultdict(set)
    for n, f, t in zip(gn, gf, gt):
        by_gn[n].add((f, t))
    
    # print(gf_by_gn)
    # print(gt_by_gn)
    gn = sorted(by_gn.keys())
    # print(gn)
    # print("|".join(sorted(gf_by_gn['B'])))
    gf = []
    gt = []
    for n in gn:
        _gf = set()
        _gt = set()
        for f, t in sorted(by_gn[n]):
            _gf.add(f)
            _gt.add(t)

        gf.append("|".join(sorted(_gf)))
        gt.append("|".join(sorted(_gt)))

    # gf = ["|".join(sorted(by_gn[n])) for n in gn]
    # gt = ["|".join(sorted(gt_by_gn[n])) for n in gn]
    return gn, gf, gt




