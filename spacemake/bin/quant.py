import pysam
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
                return

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
    uniq = set()

    def __init__(self, bam, X_channels=["exonic_counts", "intronic_counts"], **kw):
        self.kw = kw
        self.bam = bam
        self.exonic_set = set(["C", "U", "CU", "N"])
        self.intronic_set = set(["I"])

        # which channels shall contribute to the adata.X (main channel)
        self.count_X_channels = set(X_channels)

    def parse_bam_record(self, aln, tags):
        channels = set()
        gn = tags.get("gn", "").split(",")
        gf = tags.get("gf", "").split(",")
        cell = tags.get("CB", "NN")
        umi = tags.get("MI", "NN")

        gene = "NA"
        if len(gn) == 1:
            gene = gn[0]
            # we have one gene!
            key = (cell, gene, umi)
            uniq = not (key in self.uniq)
            if uniq:
                self.uniq.add(key)

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

        return gene, cell, umi, channels


    def set_reference(self, ref):
        self.logger.info(f"switching to new reference '{ref}'")


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
        help="name of the channel to be stored as AnnData.X (default='exonic_UMI')",
        default="counts",
    )
    parser.add_argument(
        "--count-X",
        help="names of the channels to be stored as AnnData.X (default='exonic_UMI')",
        default=["exonic_UMI"],
        nargs="+"
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
    # parser.add_argument(
    #     "--cell-rank-cutoff",
    #     default=100000,
    #     type=int,
    #     help="if set to positive value, keep only the top n cells, ranked by counts as specified by the counter class (see Documentation). default=100000"
    # )
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


# def read_bundles(src):
#     last_qname = None
#     batch = []
#     for rec in src:
#         if rec.query_name != last_qname:
#             if batch:
#                 yield batch
#             batch = [rec,]
#             last_qname = rec.query_name
#         else:
#             batch.append(rec)
    
#     if batch:
#         yield batch
    

# def select_alignment(alignments):

#     if len(alignments) == 1:
#         return alignments[0]


#     exonic_alignments = [aln for aln in alignments if is_exonic(aln)]
#     if len(exonic_alignments) == 1:
#         return exonic_alignments[0]
#     else:
#         return None
        

def select_alignment(src, countable_regions=set(['N', 'C', 'U', 'CU'])):
    
    def is_countable(aln):
        return set(aln.get_tag('gf').split(',')) & countable_regions

    last_qname = None
    batch = []
    for rec in src:
        if rec.is_unmapped:
            continue

        if rec.query_name != last_qname:
            if len(batch) == 1:
                yield batch[0]

            last_qname = rec.query_name
            batch = []

        if is_countable(rec):
            batch.append(rec)
    
    if len(batch) == 1:
        yield batch[0]


def main(args):
    logger = util.setup_logging(args, "spacemake.quant.main")
    util.ensure_path(args.output + "/")

    # prepare counter instance
    count_kw = vars(args)
    from_yaml = util.load_yaml(args.count_class_params)
    count_kw.update(from_yaml)
    count_class = get_counter_class(args.count_class)

    # prepare the sparse-matrix data collection (for multiple channels in parallel)
    dge = DGE(channels=count_class.channels, cell_bc_allowlist=args.cell_bc_allowlist)

    # iterate over all annotated BAMs from the input
    gene_source = {}
    for bam_name in args.bam_in:
        bam = util.quiet_bam_open(bam_name, "rb", check_sq=False, threads=4)
        
        reference_name = os.path.basename(bam_name).split(".")[0]
        counter = count_class(bam, **count_kw)
        # counter.set_reference(reference_name)

        # iterate over BAM and count
        for aln in select_alignment(util.timed_loop(bam.fetch(until_eof=True), logger, skim=args.skim)):

            tags = dict(aln.get_tags())
            # count the alignment every way that the counter prescribes
            gene, cell, umi, channels = counter.parse_bam_record(aln, tags)
            dge.add_read(gene=gene, cell=cell, channels=channels)

            # keep track of the BAM origin of every gene (in case we combine multiple BAMs)
            gene_source[gene] = reference_name

    sparse_d, obs, var = dge.make_sparse_arrays()
    adata = dge.sparse_arrays_to_adata(
        sparse_d, obs, var, main_channel=args.main_channel
    )
    adata = count_class.postprocess_adata(adata, sparse_d)
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
