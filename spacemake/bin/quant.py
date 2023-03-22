import pandas as pd
import numpy as np
import sys
import logging
import argparse
from collections import defaultdict, OrderedDict
import spacemake.util as util
from spacemake.contrib import __version__, __author__, __license__, __email__
import datetime
import os

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
default_channels = ["counts", "exonic_counts", "exonic_reads", "intronic_counts", "intronic_reads"]
default_X_counts = ["exonic_counts", "intronic_counts"] # what should be counted into adata.X matrix
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

        sparse_channels = OrderedDict()
        for channel in self.channels:
            counts = counts_by_channel[channel]
            self.logger.debug(f"making channel {channel} from counts len={len(counts)} sum={np.array(counts).sum()}")
            m = scipy.sparse.csr_matrix(
                (counts, (row_ind, col_ind))
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
        adata.obs_names = obs
        adata.var_names = var

        for channel, sparse in channel_d.items():
            if channel != main_channel:
                adata.layers[channel] = sparse

        return adata


class CountingStatistics:
    def __init__(self):
        self.stats_by_ref = defaultdict(lambda : defaultdict(int))

    def count(self, ref, name):
        self.stats_by_ref[ref][name] += 1

    def get_stats_df(self):
        dfs = []
        for ref in sorted(self.stats_by_ref.keys()):
            data = {}
            data['ref'] = ref
            for k, v in self.stats_by_ref[ref].items():
                data[k] = [v]

            dfs.append(pd.DataFrame(data))

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
        exon_intron_disambiguation = "exon_wins",
        stats=None,
        **kw):

        self.kw = kw
        self.channels = channels
        self.alignment_priorities = alignment_priorities
        self.gene_priorities = gene_priorities
        self.exonic_set = set(exonic_tags)
        self.intronic_set = set(intronic_tags)

        # which channels shall contribute to the adata.X (main channel)
        self.count_X_channels = set(X_counts)

        # how to handle reads that align to both intron and exon features
        self.exon_intron_disambiguation_func = exon_intron_disambiguation_functions[exon_intron_disambiguation]
    
        self.uniq = set()
        self.last_ref = None
        if stats is None:
            stats = CountingStatistics()

        self.stats = stats

    def set_reference(self, ref):
        if ref != self.last_ref:
            self.logger.info(f"switching to reference '{ref}'")
            # fresh counter
            self.last_ref = ref
            ## TODO extract new paramters for process_bam_bundle..
            #self.select_single_alignment_kw = lookup(ref)
            #self.select_gene_kw = lookup(ref)

    def count_stats(self, name):
        self.stats.count(self.last_ref, name)

    def select_single_alignment(self, bundle, by_score=False):
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
          
        for gn, gf, score in bundle:
            if by_score:
                p = score
            else:
                p = get_prio(gf)
            if p > top_prio:
                top_prio = p
                kept = (gn, gf)
                n_top = 1
            elif p == top_prio:
                n_top += 1
                kept = (gn, gf)

        if n_top == 1:
            return kept
        else:
            self.count_stats('N_frags_ambiguous')

    def select_gene(self, tags):
        gn, gf = tags

        gene = None
        if len(gn) > 1:
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
                gn = gn_new
                gf = gene_gf[gn[0]]
                self.count_stats('N_gene_disambiguated')
            else:
                self.count_stats('N_gene_disambiguation_failed')
                # print(gene_prio)
                # print(gene_gf)

        if len(gn) == 1:
            gene = gn[0]
            # we have one gene!

        return gene, gf

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
            channels.add("counts")
        
        for c in channels:
            self.count_stats(f'N_channel_{c}')

        if not channels:
            self.count_stats(f'N_channel_NONE')

        return channels

    def process_bam_bundle(self, cell, umi, bundle, reference_name):
        gene = None
        channels = set()
        self.set_reference(reference_name)

        if len(bundle) == 1:
            self.count_stats('N_frags_unique_mapped')
            tags = bundle[0][:2] # ignore score
        else:
            self.count_stats('N_frags_multimapped')
            # TODO: make by_score configurable
            tags = self.select_single_alignment(bundle, by_score=(reference_name == "miRNA"))
            
        if tags:
            self.count_stats('N_frags_countable')
            # count the alignment every way that the counter prescribes
            gene, gf = self.select_gene(tags)
            if gene:
                self.count_stats('N_gene_determined')

                # handle the whole uniqueness in a way that can parallelize
                # maybe split by CB[:2]? This way, distributed uniq() sets would
                # never overlap
                key = hash((cell, gene, umi))
                uniq = not (key in self.uniq)
                if uniq:
                    self.uniq.add(key)

                channels = self.determine_channels(gf, uniq)
            else:
                self.count_stats('N_gene_ambig')

        return gene, channels

    def postprocess_adata(self, adata, sparse_d):
        # number of genes detected in this cell
        adata.obs["n_genes"] = sparse_summation(adata.X > 0, axis=1)

        # add marginal counts:
        for c in self.channels:
            adata.obs[f"n_{c}"] = sparse_summation(sparse_d[c], axis=1)

        adata.obs["n_reads"] = adata.obs["n_exonic_reads"] + adata.obs["n_intronic_reads"]
        adata.obs["n_counts"] = adata.obs["n_exonic_counts"] + adata.obs["n_intronic_counts"]
        return adata


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
        "--flavor",
        help="name of the adapter flavor used to retrieve sequences and parameters from the config.yaml",
        default="default",
    )

    # parser.add_argument(
    #     "--count-class",
    #     help="full name of a class initialized with a BAM header, that needs to provide a parse_bam_record() function (see documentation). default=spacemake.quant.DefaultCounter",
    #     default="spacemake.quant.DefaultCounter",
    # )
    # parser.add_argument(
    #     "--count-class-params",
    #     help="OPTIONAL: path to a YAML file which will be parsed and the dictionary passed on to the --count-class constructor",
    #     default="",
    # )
    # parser.add_argument(
    #     "--main-channel",
    #     help="name of the channel to be stored as AnnData.X (default='counts')",
    #     default="counts",
    # )
    # parser.add_argument(
    #     "--count-X",
    #     help="names of the channels to be stored as AnnData.X (default='exonic_UMI')",
    #     default=["exonic_UMI"],
    #     nargs="+"
    # )
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
def bam_iter_bundles(bam_name, logger, args, stats):
    """
    Generator: gathers successive BAM records into bundles as long as they belong to the same read/fragment/mate pair
    (by qname). Recquired for chunking input data w/o breaking up this information.
    """

    def extract(rec):
        return (rec.get_tag('gn').split(','), rec.get_tag('gf').split(','), rec.get_tag('AS'))

    bam = util.quiet_bam_open(bam_name, "rb", check_sq=False, threads=4)
    last_qname = None
    bundle = []
    CB = 'NN'
    MI = 'NN'

    for rec in util.timed_loop(bam.fetch(until_eof=True), logger, skim=args.skim):
        stats['N_records'] += 1
        if rec.is_paired and rec.is_read2:
            continue # TODO: proper sanity checks and/or disambiguation

        if rec.is_unmapped:
            stats['N_unmapped'] += 1
            continue

        if rec.query_name != last_qname:
            stats['N_frags'] += 1

            if bundle:
                yield CB, MI, bundle

            CB = rec.get_tag('CB')
            MI = rec.get_tag('MI')
            bundle = [extract(rec),]
            last_qname = rec.query_name
        else:
            bundle.append(extract(rec))

    if bundle:
        yield CB, MI, bundle
        

def main(args):
    logger = util.setup_logging(args, "spacemake.quant.main")
    util.ensure_path(args.output + "/")

    # prepare counter instance
    # count_kw = vars(args)

    # from_yaml = util.load_yaml(args.count_class_params)
    # count_kw.update(from_yaml)

    stats = CountingStatistics()
    config = args.config['quant'][args.flavor]
    count_class = get_counter_class(config['counter_class'])
    counter = count_class(stats=stats, **config)

    # prepare the sparse-matrix data collection (for multiple channels in parallel)
    dge = DGE(channels=counter.channels, cell_bc_allowlist=args.cell_bc_allowlist)

    # iterate over all annotated BAMs from the input
    gene_source = {}
    for bam_name in args.bam_in:
        reference_name = os.path.basename(bam_name).split(".")[0]
        counter.set_reference(reference_name)

        # iterate over BAM and count
        for cell, umi, bundle in bam_iter_bundles(bam_name, logger, args, stats.stats_by_ref[reference_name]):
            # all the work done in the counter instance should be delegated to workers
            # uniq() is tricky. Perhaps only one worker for now...
            gene, channels = counter.process_bam_bundle(cell, umi, bundle, reference_name)
            if gene:
                # print(f"add_read gene={gene} cell={cell} channels={channels}")
                dge.add_read(gene=gene, cell=cell, channels=channels)

            # keep track of the BAM origin of every gene (in case we combine multiple BAMs)
            gene_source[gene] = reference_name

            if stats.stats_by_ref[reference_name]['N_records'] % 1000000 == 0:
                print(stats.get_stats_df().T)

    sparse_d, obs, var = dge.make_sparse_arrays()
    adata = dge.sparse_arrays_to_adata(
        sparse_d, obs, var
    )
    adata = counter.postprocess_adata(adata, sparse_d)
    if args.out_stats:
        stats.save_stats(args.out_stats.format(**locals()))

    # if args.cell_rank_cutoff > 0:
    #     counts = adata.obs[count_class.rank_cells_by].sort_values(ascending=False)
    #     keep = counts.iloc[:args.cell_rank_cutoff].index
    #     adata = adata[keep, :]

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

    import pandas as pd

    df = pd.DataFrame(data)
    df.to_csv(tname, sep="\t")


if __name__ == "__main__":
    args = parse_cmdline()
    main(args)
