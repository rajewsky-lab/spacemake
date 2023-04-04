import pandas as pd
import numpy as np
import sys
import re
import logging
import pysam
import argparse
from collections import defaultdict, OrderedDict
import spacemake.util as util
from spacemake.contrib import __version__, __author__, __license__, __email__
import datetime
import os
import multiprocessing as mp
from spacemake.parallel import (
    put_or_abort,
    queue_iter,
    join_with_empty_queues,
    chunkify,
    ExceptionLogging,
    log_qerr,
)
import spacemake.util as util

"""
This module implements the functionality needed to aggregate annotated alignments from a BAM file
into a digital gene expression (DGE) matrix.
The matrix has gene names as columns and cell barcodes as rows. It is stored as a scipy.sparse.csr_matrix
to save RAM. Depending on the kinds of features a read can be annotated with, counting goes through
two stages of disambiguation (see below) and then a determination of "how to count" the read. A read can be 
counted into more than one count-matrix which will be stored as 'layers' in the resulting AnnData object
(see below).

Configuration of the counting rules is done via the config.yaml. All parameters below have a 
corresponding entry there. The section "quant" holds is subdivided by flavors (for example "default") 
and each flavor can set the values of any of the options below. If you choose to customize further,
you can even set a different (sub-)class to be used for counting and potentially define additional paramters
to configure it in the YAML as well.

## BAM ingestion

BAM records are first grouped into 'bundles' of alternative alignments for the same read. Each bundle
consists of a cell barcode, a UMI sequence and a list of one or more possible alignments. For the purpose
of counting gene expression, the different alignments are represented by the values of the gn, gf, and gt 
BAM tags.

## alignment disambiguation

If we have more than one alignment in a bundle, the next step is to try and disambiguate between them.
For this purpose, the alignments are assigned priorites from the 'alignment_priorities' dictionary.
If there is a tie with multiple alignments having the same priority, the disambiguation failed and
the read can not be counted (for example aligning to two paralogous genes with exactly the same alignment
score). However, if there is one alignment with the highest priority, disambiguation has succeeded 
(for example an exonic region of one specific gene being preferred over alternative intergenic 
alignments).

## gene disambiguation

Even if a single alignment can be selected, there can still be conflicting annotations for a read: it 
may align to a region that is annotated for more than one overlapping genes. In this case another round
of disambiguation is attempted, this time using the 'gene_priorities' dictionary.
Here, for example, a coding exon annotation in the sense direction of the read would take precedence
over an intronic annotation to another gene on the opposite strand.

## channel determination

Lastly, if a read has been found to be assignable to one gene, we may still count it into different
channels (note, within each channel, the sum of all counts will never exceed the number of reads
Reads are never counted multiple times into the same channel. That includes the main adata.X channel).

For example, if a read aligns in an intronic region, it can be counted into 'intronic_reads'. If it is
furthermore the first time we observe this combination of cell barcode and UMI, it would additionally be 
counted as 'intronic_counts' ('counts' meaning "unique molecules").
If 'intronic_counts' are included in the "default_X_counts" parameters, it would also count to the main
DGE matrix stored as adata.X .

## A note on parallelization

In order to decide if a read has been seen more than once already (copies of the same BC & UMI sequences)
each DefaultCounter instance maintains a set() of hash values. To ensure correctness but still allow
parallelization, the input process will place bundles into distinct queues based on the cell barcode
sequence.

For example, if we want to run 4 counting-processes in parallel, we would map the first nucleotide to the 
processes such that barcodes starting with 'A'-> process 1, 'C' -> process 2, 'G' -> process 3, 
'T' -> process 3. In this manner the separate sets of hash values are guaranteed to never overlap and the
test for UMI uniqueness is not an issue.


"""
default_channels = ["counts", "reads", "exonic_counts", "exonic_reads", "intronic_counts", "intronic_reads"]
default_X_counts = ["exonic_counts", "intronic_counts"] # what should be counted into adata.X matrix
default_X_reads = ["exonic_reads", "intronic_reads"] # what should be counted into adata.layers["reads"] matrix (correspond with adata.X but for reads not UMIs)
default_alignment_priorities = {
    'C': 101, # coding exon
    'c': 100, # coding exon (lower case == antisense)
    'U': 51,  # UTR exon
    'u': 50,
    'CU': 51, # overlaps both, CDS+UTR (should in fact never occur as 'CU')
    'cu': 50,
    'N': 21, # exon of non-coding transcript
    'n': 20,
    'I': 11, # intronic region
    'i': 10,
    '-': 0,
}
default_gene_priorities = {
    'C': 101, # coding exon
    'c': 100, # coding exon (lower case == antisense)
    'U': 51,  # UTR exon
    'u': 50,
    'CU': 51, # overlaps both, CDS+UTR (should in fact never occur as 'CU')
    'cu': 50,
    'N': 21, # exon of non-coding transcript
    'n': 20,
    'I': 11, # intronic region
    'i': 10,
    '-': 0,
}
default_exonic_tags = ["C", "U", "CU", "N", "c", "u", "cu", "n"]
default_intronic_tags = ["I", "i"]
exon_intron_disambiguation_functions = {
    "exon_wins": lambda channels: channels - set(["intronic_reads", "intronic_counts"]),
    "intron_wins": lambda channels: channels - set(["exonic_reads", "exonic_counts"]),
    "count_both": lambda channels: channels,
}


class DGE:
    logger = logging.getLogger("spacemake.quant.DGE")

    def __init__(self, channels=["count"], cell_bc_allowlist=""):
        self.channels = channels
        self.counts = defaultdict(lambda: defaultdict(int))
        self.DGE_cells = set()
        self.DGE_genes = set()
        self.allowed_bcs = self.load_list(cell_bc_allowlist)

    @staticmethod
    def load_list(fname):
        if fname:
            allowed_bcs = set([line.rstrip() for line in open(fname)])
            DGE.logger.debug(f"restricting to {len(allowed_bcs)} allowed cell barcodes")
            return allowed_bcs
        else:
            DGE.logger.debug("all barcodes allowed")

    def add_read(self, gene, cell, channels=["count"]):
        if channels:
            if self.allowed_bcs and not (cell in self.allowed_bcs):
                # print(f"discarding '{cell}' as NA", cell in self.allowed_bcs)
                cell = "NA"

            self.DGE_cells.add(cell)
            self.DGE_genes.add(gene)

            for channel in channels:
                self.counts[(gene, cell)][channel] += 1

    def make_sparse_arrays(self):
        """_summary_
        Convert counts across all channels into a dictionary of sparse arrays

        Returns:
            dict: channel is key, value is scipt.sparse.csr_array
        """
        import scipy.sparse
        import anndata

        ind_of_cell = {}
        obs = []
        for cell in sorted(self.DGE_cells):
            ind_of_cell[cell] = len(obs)
            obs.append(cell)

        ind_of_gene = {}
        var = []
        for gene in sorted(self.DGE_genes):
            ind_of_gene[gene] = len(var)
            var.append(gene)

        counts_by_channel = defaultdict(list)
        row_ind = []
        col_ind = []
        for (gene, cell), channel_dict in self.counts.items():
            for channel in self.channels:
                counts_by_channel[channel].append(channel_dict[channel])

            row_ind.append(ind_of_cell[cell])
            col_ind.append(ind_of_gene[gene])

        if (not len(row_ind)) or (not len(col_ind)):
            self.logger.warning(f"empty indices len(row_ind)={len(row_ind)} len(col_ind)={len(col_ind)}")
            return {}, [], []

        sparse_channels = OrderedDict()
        for channel in self.channels:
            counts = counts_by_channel[channel]
            self.logger.debug(f"constructing sparse-matrix for channel '{channel}' from n={len(counts)} non-zero entries (sum={np.array(counts).sum()})")
            m = scipy.sparse.csr_matrix(
                (counts, (row_ind, col_ind)),
                dtype=np.float32
            )
            self.logger.debug(f"resulting sparse matrix shape={m.shape} len(obs)={len(obs)} len(var)={len(var)}")
            sparse_channels[channel] = m

        return sparse_channels, obs, var

    @staticmethod
    def sparse_arrays_to_adata(channel_d, obs, var, main_channel='counts'):
        """_summary_
        Creates an AnnData object from sparse matrices in the channel_d dictionary. The main_channel
        will be assigned to adata.X, the others to adata.layers[channel].

            obs (cell barcodes) -> rows
            var (gene names) -> columns

        Args:
            channel_d (_type_): _description_
            obs (_type_): _description_
            var (_type_): _description_
            main_channel (_type_): _description_

        Returns:
            AnnData object: a proper, scanpy-compatible AnnData object with (optional) layers.
        """
        import anndata

        X = channel_d[main_channel]
        adata = anndata.AnnData(X, dtype=X.dtype)
        adata.obs[f'n_counts'] = sparse_summation(adata.X, axis=1)
        adata.obs_names = obs
        adata.var_names = var

        for channel, sparse in channel_d.items():
            if channel != main_channel:
                adata.layers[channel] = sparse
                # add marginal counts
                adata.obs[f'n_{channel}'] = sparse_summation(adata.layers[channel], axis=1)

        # number of genes detected in this cell
        adata.obs["n_genes"] = sparse_summation(adata.X > 0, axis=1)
        return adata


class CountingStatistics:
    def __init__(self):
        self.stats_by_ref = defaultdict(defaultdict(int).copy)

    def count(self, ref, name):
        self.stats_by_ref[ref][name] += 1

    def add_other_stats(self, other):
        for ref, stats in other.items():
            for k, v in stats.items():
                self.stats_by_ref[ref][k] += v

    def get_stats_df(self):
        dfs = []
        for ref in sorted(self.stats_by_ref.keys()):
            data = {}
            data['ref'] = ref
            for k, v in self.stats_by_ref[ref].items():
                data[k] = [v]

            df = pd.DataFrame(data)
            dfs.append(df[sorted(df.columns)])

        if not dfs:
            dfs = [pd.DataFrame({})]

        return pd.concat(dfs)

    def save_stats(self, path, sep='\t', **kw):
        self.get_stats_df().to_csv(path, sep=sep, **kw)


class DefaultCounter:

    logger = logging.getLogger("spacemake.quant.DefaultCounter")

    def __init__(self, 
        channels = default_channels,
        alignment_priorities = default_alignment_priorities,
        gene_priorities = default_gene_priorities,
        exonic_tags = default_exonic_tags,
        intronic_tags = default_intronic_tags,
        X_counts=default_X_counts,
        X_reads=default_X_reads,
        exon_intron_disambiguation = "exon_wins",
        alignment_selection="priority",
        gene_selection="priority",
        uniq=set(),
        stats=CountingStatistics(),
        **kw):

        self.kw = kw
        self.channels = channels
        self.alignment_priorities = alignment_priorities
        self.gene_priorities = gene_priorities
        self.exonic_set = set(exonic_tags)
        self.intronic_set = set(intronic_tags)

        # how to choose among multiple alignments
        self.select_alignment = {
            'priority': self.select_alignment_by_priority,
            'take_first': self.select_first_alignment,
            'take_first_plus': self.select_first_plus_alignment,
        }[alignment_selection]

        # how to choose among multiple genes
        self.select_gene = {
            'priority': self.select_gene_by_priority,
            'chrom': self.select_chrom_as_gene,
        }[gene_selection]

        # which channels shall contribute to the adata.X (main channel)
        self.count_X_channels = set(X_counts)
        self.read_X_channels = set(X_reads)

        # how to handle reads that align to both intron and exon features
        self.exon_intron_disambiguation_func = exon_intron_disambiguation_functions[exon_intron_disambiguation]
    
        self.stats = stats
        self.uniq = uniq

    ## Alignment selection strategies
    def select_first_alignment(self, bundle):
        return bundle[0]

    def select_first_plus_alignment(self, bundle):
        # only consider alignments on the + strand (for custom indices)
        plus = [b for b in bundle if b[1] == '+']
        if plus:
            return plus[0]


    def select_alignment_by_priority(self, bundle):
        """
        If multiple alignments are reported for one fragment/read try to make a reasonable choice:
        If one of the alignments is to a coding exon, but the others are not, choose the exonic one.
        If multiple alternative alignments point to exons of different coding genes, do not count anything
        because we can not disambiguate.
        """        
        from collections import defaultdict

        def get_prio(gf):
            prios = [self.alignment_priorities.get(f, 0) for f in gf]
            return max(prios)

        top_prio = 0
        n_top = 0
        kept = None
          
        for chrom, strand, gn, gf, score in bundle:
            p = get_prio(gf)
            if p > top_prio:
                top_prio = p
                kept = (chrom, strand, gn, gf, score)
                n_top = 1
            elif p == top_prio:
                n_top += 1
                kept = (chrom, strand, gn, gf, score)

        if n_top == 1:
            return kept

    ## Gene selection strategies
    def select_chrom_as_gene(self, chrom, strand, gn, gf, score):
        return chrom, ['N']

    def select_gene_by_priority(self, chrom, strand, gn, gf, score):
        # let's see if we can prioritize which gene we are interested in
        gene_prio = defaultdict(int)
        gene_gf = defaultdict(list)
        max_prio = 0
        for n, f in zip(gn, gf):
            p = self.gene_priorities.get(f, 0)
            gene_prio[n] = max(gene_prio[n], p)
            max_prio = max([max_prio, p])
            gene_gf[n].append(f)
        
        # restrict to genes tied for highest priority hits
        gn_new = [n for n in set(gn) if gene_prio[n] == max_prio]
        # print(len(gn_new))
        if len(gn_new) == 1:
            # YES we can salvage this read
            gene = gn_new[0]
            return gene, gene_gf[gene]
        else:
            # NO still ambiguous
            return None, None

    ## Count-channel determination
    def determine_channels(self, gf, uniq):
        channels = set()
        exon = False
        intron = False
        for f in gf:
            if f in self.exonic_set:
                channels.add("exonic_reads")
                exon = True
                if uniq:
                    channels.add("exonic_counts")

            elif f in self.intronic_set:
                channels.add("intronic_reads")
                intron = True
                if uniq:
                    channels.add("intronic_counts")

        # post-process: if both exon & intron annotations 
        # (but from different isoforms) are present
        # decide how to count
        if exon and intron:
            channels = self.exon_intron_disambiguation_func(channels)

        if channels & self.count_X_channels:
            channels.add("reads")
            if uniq:
                channels.add("counts")
        
        return channels

    ## main function: alignment bundle -> counting channels
    def process_bam_bundle(self, cell, umi, bundle):
        gene = None
        selected = None
        channels = set()
        # self.set_reference(reference_name)
        # print(f"bundle={bundle}")
        if len(bundle) == 1:
            self.stats['N_aln_unique'] += 1
            selected = bundle[0]
        else:
            self.stats['N_aln_multi'] += 1
            selected = self.select_alignment(bundle)
            if selected:
                self.stats['N_aln_selected'] += 1
            else:
                self.stats['N_aln_selection_failed'] += 1

        if selected:
            self.stats['N_aln_countable'] += 1
            chrom, strand, gn, gf, score = selected
            # self.stats[f'N_aln_{chrom}'] += 1
            self.stats[f'N_aln_{strand}'] += 1
            if len(gn) == 1:
                gene = gn[0]
                if gene is None:
                    self.stats['N_gene_none'] += 1    
                else:
                    self.stats['N_gene_unique'] += 1
                # gf = gf
                # print(f"uniq gene: {gene} {gf}")
            else:
                self.stats['N_gene_multi'] += 1
                gene, gf = self.select_gene(*selected)
                if gene:
                    self.stats['N_gene_selected'] += 1
                else:
                    self.stats['N_gene_selection_failed'] += 1
                    
            # count the alignment every way that the counter prescribes
            if gene:
                # print(f"gene={gene} gf={gf}")
                self.stats['N_aln_counted'] += 1
                if gf[0].islower():
                    self.stats['N_aln_antisense'] += 1
                else:
                    self.stats['N_aln_sense'] += 1

                # handle the whole uniqueness in a way that can parallelize
                # maybe split by UMI[:2]? This way, distributed uniq() sets would
                # never overlap
                key = hash((cell, umi))
                uniq = not (key in self.uniq)
                if uniq:
                    self.uniq.add(key)

                channels = self.determine_channels(gf, uniq)
                for c in channels:
                    self.stats[f'N_channel_{c}'] += 1

                if not channels:
                    self.stats[f'N_channel_NONE'] += 1

        return gene, channels


def parse_cmdline():
    parser = util.make_minimal_parser(prog="quant.py", description="quantify per-cell gene expression from BAM files by counting into a (sparse) DGE matrix")
    parser.add_argument("--config", default="config.yaml", help="path to config-file")

    parser.add_argument(
        "bam_in",
        help="bam input (default=stdin)",
        default=["/dev/stdin"],
        nargs="+",
    )
    parser.add_argument(
        "--skim",
        help="skim through the BAM by investigating only every <skim>-th record (default=1 off)",
        default=1,
        type=int
        # nargs="+",
    )
    parser.add_argument(
        "--parallel",
        help="number of parallel worker processes. Needs to be divisible by 2 (default=4)",
        default=4,
        type=int
    )
    parser.add_argument(
        "--buffer-size",
        help="number of bundles per buffer (default=10000)",
        default=10000,
        type=int
    )
    parser.add_argument(
        "--flavor",
        help="name of the adapter flavor used to retrieve sequences and parameters from the config.yaml. Can also be a mapping choosing specific flavor for each reference (example: 'first@miRNA,default@genome,chrom@rRNA')",
        default="default",
    )
    parser.add_argument(
        "--cell-bc-allowlist",
        help="[OPTIONAL] a text file with cell barcodes. All barcodes not in the list are ignored.",
        default="",
    )
    parser.add_argument(
        "--output",
        help="directory to store the output h5ad and statistics/marginal counts",
        default="dge",
    )
    parser.add_argument(
        "--out-dge",
        help="filename for the output h5ad",
        default="{args.output}/{args.sample}.h5ad",
    )
    parser.add_argument(
        "--out-summary",
        help="filename for the output summary (sum over all vars/genes)",
        default="{args.output}/{args.sample}.summary.tsv",
    )
    parser.add_argument(
        "--out-bulk",
        help="filename for the output summary (sum over all vars/genes)",
        default="{args.output}/{args.sample}.pseudo_bulk.tsv",
    )
    parser.add_argument(
        "--out-stats",
        help="filename for the statistics/counts output",
        default="{args.output}/{args.sample}.stats.tsv",
    )

    return util.load_config_with_fallbacks(parser.parse_args())


def get_counter_class(classpath):
    mod, cls = classpath.rsplit(".", 1)

    if mod == "spacemake.quant":
        return globals()[cls]
    else:
        import importlib

        m = importlib.import_module(mod)
        return getattr(m, cls)


def get_config_for_refs(args):
    """
    Parse flavor definitions, filling a dictionary with reference name as key and configuration
    as values.
    The values are kwargs dictionaries, including the counter_class argument already replaced
    with the actual class object ready for instantiation.

    The key '*' is used for default settings.

    For convenience, we also collect the names of all channels known by all counter classes, so
    that the DGE object can be created before any data are parsed.
    """
    flavors = args.config['quant']
    # replace str of counter class name/path with actual class object
    for name, config in flavors.items():
        config['counter_class'] = get_counter_class(config.get('counter_class', "spacemake.quant.DefaultCounter"))
        config['name'] = name

    ref_d = {}
    default = 'default'
    for f in args.flavor.split(','):
        if '@' in f:
            ref, name = f.split('@')
            ref_d[ref] = flavors[name]
        else:
            default = f

    ref_d['*'] = flavors[default]

    # collect all channel names that can be expected to be generated
    channels = []
    for ref, config in ref_d.items():
        channels.extend(config.get('channels', default_channels))

    return ref_d, sorted(set(channels))


def sparse_summation(X, axis=0):
    # for csr_array this would be fine
    # return np.array(X.sum(axis=axis))
    # unfortunately, for the time being scanpy needs csr_matrix
    # due to h5py lagging somewhat
    if axis == 0:
        return np.array(X.sum(axis=0))[0]
    elif axis == 1:
        return np.array(X.sum(axis=1))[:, 0]



## This will become the reader process
def bam_iter_bundles(bam_src, logger, args, stats):
    """
    Generator: gathers successive BAM records into bundles as long as they belong to the same read/fragment/mate pair
    (by qname). Recquired for chunking input data w/o breaking up this information.
    """

    def extract(rec):
        return (
            bam_src.header.get_reference_name(rec.tid), # chrom
            '-' if rec.is_reverse else '+', # strand
            rec.get_tag('gn').split(',') if rec.has_tag('gn') else '-', # gene names
            rec.get_tag('gf').split(',') if rec.has_tag('gf') else '-', # gene feature overlap encodings
            rec.get_tag('AS') # alignment score
        )

    last_qname = None
    bundle = []
    bundle_ref = None
    CB = 'NN'
    MI = 'NN'

    # for ref, rec in util.timed_loop(bam_src, logger, skim=args.skim):
    for ref, rec in bam_src:
        ref_stats = stats.stats_by_ref[ref]
        ref_stats['N_records'] += 1
        if rec.is_paired and rec.is_read2:
            continue # TODO: proper sanity checks and/or disambiguation

        if rec.is_unmapped:
            ref_stats['N_unmapped'] += 1
            continue

        if rec.query_name != last_qname:
            ref_stats['N_frags'] += 1

            if bundle:
                yield bundle_ref, CB, MI, bundle

            CB = rec.get_tag('CB')
            MI = rec.get_tag('MI')
            bundle = [extract(rec),]
            last_qname = rec.query_name
            bundle_ref = ref
        else:
            bundle.append(extract(rec))

    if bundle:
        yield ref, CB, MI, bundle
        

def main(args):
    logger = util.setup_logging(args, "spacemake.quant.main")
    util.ensure_path(args.output + "/")
    # load detailed counter configurations for each reference name
    # plus a default setting ('*' key)
    conf_d, channels = get_config_for_refs(args)
    
    # prepare the sparse-matrix data collection (for multiple channels in parallel)
    dge = DGE(channels=channels, cell_bc_allowlist=args.cell_bc_allowlist)

    # keep track of which reference contained which gene names
    gene_source = {}
    
    # these objects can be (re-)used by successive counter_class instances
    stats = CountingStatistics() # statistics on disambiguation and counting
    uniq = set() # keep track of (CB, UMI) tuples we already encountered

    last_ref = None
    counter = None
    ref_stats = None
    # iterate over all annotated BAMs from the input
    for bam_name in args.bam_in:
        reference_name = os.path.basename(bam_name).split(".")[0]
        if last_ref != reference_name:
            # create a new counter-class instance with the configuration
            # for this reference name (genome, miRNA, rRNA, ...)
            last_ref = reference_name
            ref_stats = stats.stats_by_ref[reference_name]
            config = conf_d.get(reference_name, conf_d['*'])
            logger.info(f"processing alignments to reference '{reference_name}'. Building counter with '{config['name']}' rules")
            counter = config['counter_class'](stats=ref_stats, uniq=uniq, **config)

        # iterate over BAM and count
        for cell, umi, bundle in bam_iter_bundles(bam_name, logger, args, stats=ref_stats):
            # all the work done in the counter instance should be delegated to workers
            # uniq() is tricky. Perhaps only one worker for now...
            gene, channels = counter.process_bam_bundle(cell, umi, bundle)
            if gene:
                # print(f"add_read gene={gene} cell={cell} channels={channels}")
                dge.add_read(gene=gene, cell=cell, channels=channels)

            # keep track of the BAM origin of every gene (in case we combine multiple BAMs)
            gene_source[gene] = reference_name

            if ref_stats['N_records'] % 1000000 == 0:
                print(stats.get_stats_df().T)

    sparse_d, obs, var = dge.make_sparse_arrays()
    adata = dge.sparse_arrays_to_adata(
        sparse_d, obs, var
    )
    if args.out_stats:
        stats.save_stats(args.out_stats.format(**locals()))

    adata.var["reference"] = [gene_source[gene] for gene in var]
    adata.uns["sample_name"] = args.sample
    adata.uns[
        "DGE_info"
    ] = f"created with spacemake.quant version={__version__} on {datetime.datetime.today().isoformat()}"
    adata.uns["DGE_cmdline"] = sys.argv
    adata.write(args.out_dge.format(args=args))

    ## write out marginal counts across both axes:
    #   summing over obs/cells -> pseudo bulk
    #   summing over genes -> UMI/read distribution
    pname = args.out_bulk.format(args=args)

    data = OrderedDict()

    data["sample_name"] = args.sample
    data["reference"] = [gene_source[gene] for gene in var]
    data["gene"] = var
    for channel, M in sparse_d.items():
        data[channel] = sparse_summation(M)

    import pandas as pd

    df = pd.DataFrame(data).sort_values('counts')
    df.to_csv(pname, sep="\t")

    tname = args.out_summary.format(args=args)

    data = OrderedDict()

    # TODO: add back in once we dropped the Rmd QC sheet scripts
    # data["sample_name"] = args.sample
    data["cell_bc"] = obs
    for channel, M in sparse_d.items():
        data[channel] = sparse_summation(M, axis=1)

    df = pd.DataFrame(data)
    df.to_csv(tname, sep="\t")


class SplitBAMbyUMI:
    def __init__(self, Qin, buf_size=10000, k=4):
        self.logger = logging.getLogger("spacemake.quant.SplitBAMbyUMI")
        self.Qin = Qin
        self.n = len(Qin)
        self.k = k
        self.regexp = re.compile(f'MI:Z:({"." * self.k})')

        self.buffers = [list() for n in range(self.n)]
        self.buffer_map = {}

        from spacemake.util import generate_kmers
        for i, kmer in enumerate(generate_kmers(k, nts='ACGTN')):
            j = i % self.n
            self.buffer_map[kmer] = self.buffers[j]
            self.logger.debug(f"assigning {kmer} -> worker {j}")

        self.logger.debug(f"splitting reads by first 4 bases of the UMI to assign to {self.n} queues")
        self.max_buf_size = buf_size
        self.n_chunks = 0

    def push_to_queues(self, ref):
        sizes = []
        # print(f"push_to_queues, buffers={self.buffers}")
        for buf, Q in zip(self.buffers, self.Qin):
            sizes.append(len(buf))
            if buf:
                # print("pushing buffer", buf)
                # MJ: just spent the better part of an hour debugging 
                # to realize there is some async black magic happening 
                # in Q.put() or whatever.
                # cause if you send the actual buf and then clear() it
                # right after put() the receiver gets truncated or no data.
                # Yeah. For. real.
                # so buf.copy() it is... :-/
                Q.put((self.n_chunks, ref, buf.copy()))
                buf.clear()
                self.n_chunks += 1

        return np.array(sizes)
  
    def broadcast(self, rec):
        # print(f"broadcast, buffers={self.buffers}")
        for buf in self.buffers:
            buf.append(rec)

    def shutdown(self):
        self.logger.info("sending None to all queues to signal end of input")
        for Q in self.Qin:
            Q.put(None)

    def iter_bam(self, bam_name, ref='na'):
        import subprocess
        import re

        proc = subprocess.Popen(['samtools', 'view', '-h', '--no-PG', '--threads=4', bam_name], stdout=subprocess.PIPE, text=True)
        while line := proc.stdout.readline():
            if line.startswith('@'):
                self.broadcast(line)
            else:
                M = re.search(self.regexp, line)
                buf = self.buffer_map[M.groups()[0]]
                # print(f"appending to buffer of L={len(buf)}")
                buf.append(line)

                if len(buf) > self.max_buf_size:
                    yield self.push_to_queues(ref)

        yield self.push_to_queues(ref)


class AlignedSegmentsFromQueue:
    def __init__(self, Qin, abort_flag):
        self.Qin = Qin
        self.abort_flag = abort_flag
        self.header = None
        self.header_lines = []
        self.last_ref = None

    def __iter__(self):
        for n_chunk, ref, sam_lines in queue_iter(self.Qin, self.abort_flag):
            if ref != self.last_ref:
                self.header = None
                self.header_lines = []
                self.last_ref = ref

            # print(f"received chunk {n_chunk} with {len(sam_lines)} lines: {sam_lines}")
            for line in sam_lines:
                # print(line)
                if line.startswith('@'):
                    self.header_lines.append(line)
                else:
                    if not self.header:
                        # print("about to construct header:")
                        # print(self.header_lines)
                        self.header = pysam.AlignmentHeader.from_text("".join(self.header_lines))
                        # print(self.header)
                    else:
                        sam = pysam.AlignedSegment.fromstring(line, self.header)
                        # print(sam)
                        yield ref, sam
            
            self.Qin.task_done()

# def AlignedSegmentsFromFile(bam_name):
#     bam = util.quiet_bam_open(bam_name, "rb", check_sq=False, threads=4)
#     for rec in bam.fetch(until_eof=True):
#         yield rec


def BAM_reader(Qin, args, Qerr, abort_flag, stat_list):
    with ExceptionLogging(
        "spacemake.quant.BAM_reader", Qerr=Qerr, exc_flag=abort_flag
    ) as el:
        from time import time
        stats = CountingStatistics() # statistics on disambiguation and counting
        ref_stats = None
        last_ref = None
        dispatch = SplitBAMbyUMI(Qin, buf_size=args.buffer_size)

        # iterate over all BAMs with countable alignments
        N = 0
        T0 = time()
        for bam_name in args.bam_in:
            reference_name = os.path.basename(bam_name).split(".")[0]
            if last_ref != reference_name:
                # create a new counter-class instance with the configuration
                # for this reference name (genome, miRNA, rRNA, ...)
                last_ref = reference_name
                # ref_stats = stats.stats_by_ref[reference_name]
                t0 = time()
                t1 = time()
                N_ref = 0
                el.logger.info(f"reading alignments to reference '{reference_name}'.")

            # split-by-UMI dispatch
            for n_pushed in dispatch.iter_bam(bam_name, ref=reference_name):
                n = n_pushed.sum()
                N += n
                N_ref += n
                # ref_stats['N_records'] = N_ref
                if time() - t1 > 2:
                    t = time() - t0
                    el.logger.debug(f"worker load distribution: {n_pushed / float(n)}")
                    el.logger.info(f"ingested {N_ref} BAM records in {t:.1f} seconds ({0.001 * N_ref/t:.2f} k/sec).")
                    t1 = time()

        t = time() - T0
        el.logger.debug(f"worker load distribution: {n_pushed / float(n)}")
        el.logger.info(f"ingested {N} BAM records in {t:.1f} seconds ({0.001 * N/t:.2f} k/sec).")
        t1 = time()

        el.logger.debug("done. closing queues.")
        dispatch.shutdown()
        
        el.logger.debug("syncing stats...")
        stat_list.append(stats.stats_by_ref)
        
        el.logger.debug("shutting down...")


def bundle_processor(Qin, Qout, args, Qerr, abort_flag, stat_list):
    with ExceptionLogging(
        "spacemake.quant.bundle_processor", Qerr=Qerr, exc_flag=abort_flag
    ) as el:
        # load detailed counter configurations for each reference name
        # plus a default setting ('*' key)
        conf_d, channels = get_config_for_refs(args)

        # these objects can be (re-)used by successive counter_class instances
        stats = CountingStatistics() # statistics on disambiguation and counting
        uniq = set() # keep track of (CB, UMI) tuples we already encountered

        last_ref = None
        counter = None
        ref_stats = None
        n_chunk = 0
        count_data_chunk = []

        bam_src = AlignedSegmentsFromQueue(Qin, abort_flag)
        for ref, cell, umi, bundle in bam_iter_bundles(bam_src, el.logger, args, stats=stats):
            if last_ref != ref:
                # create a new counter-class instance with the configuration
                # for this reference name (genome, miRNA, rRNA, ...)
                last_ref = ref
                ref_stats = stats.stats_by_ref[ref]
                config = conf_d.get(ref, conf_d['*'])
                el.logger.info(f"processing alignments to reference '{ref}'. Building counter with '{config['name']}' rules")
                counter = config['counter_class'](stats=ref_stats, uniq=uniq, **config)

            # try to select alignment and gene and determine into which
            # channels we want to count
            gene, channels = counter.process_bam_bundle(cell, umi, bundle)
            if gene:
                count_data_chunk.append((gene, cell, channels))
                # print(f"add_read gene={gene} cell={cell} channels={channels}")
                # dge.add_read(gene=gene, cell=cell, channels=channels)

            if len(count_data_chunk) >= args.buffer_size:
                Qout.put((n_chunk, ref, count_data_chunk))
                n_chunk += 1
                count_data_chunk = []
        
        if count_data_chunk:
            Qout.put((n_chunk, ref, count_data_chunk))
            n_chunk += 1
            count_data_chunk = []

        el.logger.debug("done. syncing stats and gene_source")
        stat_list.append(stats.stats_by_ref)

        el.logger.debug("shutting down")


def DGE_counter(Qin, args, Qerr, abort_flag, stat_list):
    with ExceptionLogging(
        "spacemake.quant.DGE_counter", Qerr=Qerr, exc_flag=abort_flag
    ) as el:
        # load detailed counter configurations for each reference name
        # plus a default setting ('*' key)
        conf_d, channels = get_config_for_refs(args)
        
        # prepare the sparse-matrix data collection (for multiple channels in parallel)
        dge = DGE(channels=channels, cell_bc_allowlist=args.cell_bc_allowlist)

        # keep track of which reference contained which gene names
        gene_source = {}
    
        # iterate over chunks of count_data prepared by the worker processes
        for n_chunk, ref, count_data in queue_iter(Qin, abort_flag):
            for gene, cell, channels in count_data:
                dge.add_read(gene=gene, cell=cell, channels=channels)
                # keep track of the mapping reference origin of every gene 
                # (in case we combine multiple BAMs)
                gene_source[gene] = ref

            Qin.task_done()

        el.logger.debug("completed counting. Writing output.")

        sparse_d, obs, var = dge.make_sparse_arrays()
        adata = dge.sparse_arrays_to_adata(
            sparse_d, obs, var
        )
        adata.var["reference"] = [gene_source[gene] for gene in var]
        adata.uns["sample_name"] = args.sample
        adata.uns[
            "DGE_info"
        ] = f"created with spacemake.quant version={__version__} on {datetime.datetime.today().isoformat()}"
        adata.uns["DGE_cmdline"] = sys.argv
        adata.write(args.out_dge.format(args=args))

        ## write out marginal counts across both axes:
        #   summing over obs/cells -> pseudo bulk
        #   summing over genes -> UMI/read distribution
        pname = args.out_bulk.format(args=args)

        data = OrderedDict()

        data["sample_name"] = args.sample
        data["reference"] = [gene_source[gene] for gene in var]
        data["gene"] = var
        for channel, M in sparse_d.items():
            data[channel] = sparse_summation(M)

        import pandas as pd

        df = pd.DataFrame(data).sort_values('counts')
        df.to_csv(pname, sep="\t")

        tname = args.out_summary.format(args=args)

        data = OrderedDict()

        # TODO: add back in once we dropped the Rmd QC sheet scripts
        # data["sample_name"] = args.sample
        data["cell_bc"] = obs
        for channel, M in sparse_d.items():
            data[channel] = sparse_summation(M, axis=1)

        df = pd.DataFrame(data)
        df.to_csv(tname, sep="\t")
        el.logger.debug("shutdown")

def main_parallel(args):
    logger = util.setup_logging(args, "spacemake.quant.main_parallel")
    util.ensure_path(args.output + "/")

    # queues for communication between processes
    Qin = [mp.JoinableQueue(args.parallel * 10) for i in range(args.parallel)]
    Qout = mp.JoinableQueue()
    Qerr = mp.Queue()  # child-processes can report errors back to the main process here

    # Proxy objects to allow workers to report statistics about the run
    with mp.Manager() as manager:
        abort_flag = mp.Value("b")
        abort_flag.value = False

        stat_list = manager.list()
        with ExceptionLogging(
            "spacemake.quant.main_parallel", exc_flag=abort_flag
        ) as el:

            # read BAM in chunks and put them in Qsam
            dispatcher = mp.Process(
                target=BAM_reader,
                name="BAM_reader",
                args=(Qin, args, Qerr, abort_flag, stat_list),
            )

            dispatcher.start()
            el.logger.debug("Started dispatch")

            # workers consume chunks of BAM from Qsam
            # process them, and put the results in Qres
            workers = []
            for i in range(args.parallel):
                w = mp.Process(
                    target=bundle_processor,
                    name=f"worker_{i}",
                    args=(Qin[i], Qout, args, Qerr, abort_flag, stat_list),
                )
                w.start()
                workers.append(w)

            el.logger.debug("Started workers")

            collector = mp.Process(
                target=DGE_counter,
                name="output",
                args=(Qout, args, Qerr, abort_flag, stat_list),
            )
            collector.start()
            el.logger.debug("Started collector")
            # wait until all sequences have been thrown onto Qfq
            qerr = join_with_empty_queues(dispatcher, Qin + [Qerr], abort_flag, logger=el.logger)
            el.logger.debug("The dispatcher exited")

            if qerr[-1]:
                el.logger.info(f"{len(qerr)} chunks were drained from Qfq upon abort.")
                log_qerr(qerr[-1])

            # wait until the workers have completed every task and placed
            # every output chunk onto Qres.

            import time
            el.logger.debug("Waiting for workers to process all data")
            for Q, w in zip(Qin, workers):
                while Q.qsize() > 0:
                    time.sleep(.1)
                # el.logger.debug(f"Q-size: {Q.qsize()}")
                # w.join()
                # el.logger.debug(f"Q-size: {Q.qsize()}")
                # el.logger.debug(".")

            el.logger.debug("Waiting for writer to accumulate all data")
            # Qout.join()
            while Qout.qsize():
                time.sleep(.1)
                # print(Qout.qsize())

            el.logger.debug("Telling writer to stop.")
            Qout.put(None)


            collector.join()
            el.logger.debug("Collector has joined.")

            # el.logger.debug("Joining input queues")
            # [Q.join() for Q in Qin]

            el.logger.debug("Waiting for workers to exit")
            for i, w in enumerate(workers):
                w.join()

            el.logger.debug(
                "All worker processes have joined. Gathering statistics."
            )

            
            stats = CountingStatistics()
            for s in stat_list:
                stats.add_other_stats(s)
                
            if args.out_stats:
                stats.save_stats(args.out_stats.format(**locals()))

        if el.exception:
            return -1



if __name__ == "__main__":
    args = parse_cmdline()
    if args.parallel > 0:
        main_parallel(args)
    else:
        main(args)
