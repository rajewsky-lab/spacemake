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
            if self.allowed_bcs and not cell in self.allowed_bcs:
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
    def sparse_arrays_to_adata(channel_d, obs, var, main_channel):
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


class DefaultCounter:

    channels = ["counts", "exonic_counts", "exonic_reads", "intronic_counts", "intronic_reads"]
    rank_cells_by = "n_reads"
    logger = logging.getLogger("spacemake.quant.DefaultCounter")
    priorities={
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
    }
    exonic_set = set(["C", "U", "CU", "N", "c", "u", "cu", "n"])
    intronic_set = set(["I", "i"])

    def __init__(self, X_channels=["exonic_counts", "intronic_counts"], **kw):
        self.kw = kw
#         if not stranded:
#             self.exonic_set |= set([x.lower() for x in self.exonic_set])
#  #           self.coding_set |= set([x.lower() for x in self.coding_set])
#             self.intronic_set |= set([x.lower() for x in self.intronic_set])

        # which channels shall contribute to the adata.X (main channel)
        self.count_X_channels = set(X_channels)
    
        self.uniq = set()
        self.stats = defaultdict(int)

    def select_single_alignment(self, bundle):
        """
        If multiple alignments are reported for one fragment/read try to make a reasonable choice:
        If one of the alignments is to a coding exon, but the others are not, choose the exonic one.
        If multiple alternative alignments point to exons of different coding genes, do not count anything
        because we can not disambiguate.
        """        
        from collections import defaultdict

        def get_prio(gf):
            prios = [self.priorities.get(f, 0) for f in gf]
            return max(prios)

        top_prio = 0
        n_top = 0
        kept = None
          
        for gn, gf in bundle:
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
            self.stats['N_frags_ambiguous'] += 1


    def select_gene(self, tags):
        gn, gf = tags

        gene = None
        if len(gn) > 1:
            # let's see if we can prioritize which gene we are interested in
            gene_prio = defaultdict(int)
            gene_gf = defaultdict(list)
            max_prio = 0
            for n, f in zip(gn, gf):
                p = self.priorities.get(f, 0)
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
                self.stats['N_gene_disambiguated'] += 1
            else:
                self.stats['N_gene_disambiguation_failed'] += 1
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
        # count the read as exonic
        if exon and intron:
            channels -= set(["intronic_reads", "intronic_counts"])

        if channels & self.count_X_channels:
            channels.add("counts")
        
        for c in channels:
            self.stats[f'N_channel_{c}'] += 1

        if not channels:
            self.stats[f'N_channel_NONE'] += 1

        return channels


    def process_bam_bundle(self, cell, umi, bundle):
        gene = None
        channels = set()

        if len(bundle) == 1:
            self.stats['N_frags_unique_mapped'] += 1
            tags = bundle[1]
        else:
            self.stats['N_frags_multimapped'] += 1
            tags = self.select_single_alignment(bundle)
            
        if tags:
            self.stats['N_frags_countable'] += 1
            # count the alignment every way that the counter prescribes
            gene, gf = self.select_gene(tags)
            if gene:
                self.stats['N_gene_determined'] += 1

                # handle the whole uniqueness in a way that can parallelize
                # maybe split by CB[:2]? This way, distributed uniq() sets would
                # never overlap
                key = hash(cell, gene, umi)
                uniq = not (key in self.uniq)
                if uniq:
                    self.uniq.add(key)

                channels = self.determine_channels(gf, uniq)
            else:
                self.stats['N_gene_ambig'] += 1

        return gene, cell, channels


    def get_stats_df(self):
        data = {}
        for k, v in self.stats.items():
            data[k] = [v]

        return pd.DataFrame(data, index=['count']).T


    def save_stats(self, path, sep='\t', **kw):
        self.get_stats_df().to_csv(path, sep=sep, **kw)

    @staticmethod
    def postprocess_adata(adata, sparse_d):
        # number of genes detected in this cell
        adata.obs["n_genes"] = sparse_summation(adata.X > 0, axis=1)

        # add marginal counts:
        for c in DefaultCounter.channels:
            adata.obs[f"n_{c}"] = sparse_summation(sparse_d[c], axis=1)

        adata.obs["n_reads"] = adata.obs["n_exonic_reads"] + adata.obs["n_intronic_reads"]
        adata.obs["n_counts"] = adata.obs["n_exonic_counts"] + adata.obs["n_intronic_counts"]
        return adata


def parse_cmdline():
    parser = util.make_minimal_parser(prog="quant.py", description="quantify per-cell gene expression from BAM files by counting into a (sparse) DGE matrix")

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
        "--count-class",
        help="full name of a class initialized with a BAM header, that needs to provide a parse_bam_record() function (see documentation). default=spacemake.quant.DefaultCounter",
        default="spacemake.quant.DefaultCounter",
    )
    parser.add_argument(
        "--count-class-params",
        help="OPTIONAL: path to a YAML file which will be parsed and the dictionary passed on to the --count-class constructor",
        default="",
    )
    parser.add_argument(
        "--main-channel",
        help="name of the channel to be stored as AnnData.X (default='counts')",
        default="counts",
    )
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


    return parser.parse_args()


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
        return (rec.get_tag('gn').split(','), rec.get_tag('gf').split(','))

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

            CB = rec.get_tag('CB'),
            MI = rec.get_tag('MI'),
            bundle = [extract(rec),]
            last_qname = rec.query_name
        else:
            bundle.append(extract(rec))

        if stats['N_records'] % 1000000 == 0:
            for k, v in sorted(stats.items()):
                print(f"{k}\t{v/1000:.1f}k")

    if bundle:
        yield CB, MI, bundle
        

def main(args):
    logger = util.setup_logging(args, "spacemake.quant.main")
    util.ensure_path(args.output + "/")

    # prepare counter instance
    count_kw = vars(args)
    from_yaml = util.load_yaml(args.count_class_params)
    count_kw.update(from_yaml)
    count_class = get_counter_class(args.count_class)
    counter = count_class(**count_kw)

    # prepare the sparse-matrix data collection (for multiple channels in parallel)
    dge = DGE(channels=count_class.channels, cell_bc_allowlist=args.cell_bc_allowlist)

    # iterate over all annotated BAMs from the input
    gene_source = {}
    for bam_name in args.bam_in:
        reference_name = os.path.basename(bam_name).split(".")[0]
        # counter.set_reference(reference_name)

        # iterate over BAM and count
        for cell, umi, bundle in bam_iter_bundles(bam_name, logger, args, counter.stats):
            # all the work done in the counter instance should be delegated to workers
            # uniq() is tricky. Perhaps only one worker for now...
            gene, channels = counter.process_bam_bundle(cell, umi, bundle)
            if gene:
                dge.add_read(gene=gene, cell=cell, channels=channels)

            # keep track of the BAM origin of every gene (in case we combine multiple BAMs)
            gene_source[gene] = reference_name

    sparse_d, obs, var = dge.make_sparse_arrays()
    adata = dge.sparse_arrays_to_adata(
        sparse_d, obs, var, main_channel=args.main_channel
    )
    adata = count_class.postprocess_adata(adata, sparse_d)
    if args.out_stats:
        count_class.save_stats(args.out_stats)

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

    df = pd.DataFrame(data).sort_values(args.main_channel)
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
