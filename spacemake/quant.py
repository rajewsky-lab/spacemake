import pysam
import numpy as np
import sys
import logging
import argparse
from collections import defaultdict, OrderedDict
import spacemake.util as util
import os


class DGE:
    def __init__(self, channels=["count"]):
        self.channels = channels
        self.counts = defaultdict(lambda: defaultdict(int))
        self.DGE_cells = set()
        self.DGE_genes = set()

    def add_read(self, gene, cell, channels=["count"]):
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
            sparse_channels[channel] = scipy.sparse.csr_array(
                (counts, (row_ind, col_ind))
            )
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

        adata = anndata.AnnData(channel_d[main_channel])
        adata.obs_names = obs
        adata.var_names = var

        for channel, sparse in channel_d.items():
            if channel != main_channel:
                adata.layers[channel] = sparse

        return adata


class DefaultCounter:
    def __init__(self, bam, **kw):
        self.kw = kw
        self.bam = bam
        self.channels = ["exonic_UMI", "exonic_read", "intronic_UMI", "intronic_read"]
        self.exonic_set = set(["C", "U", "CU", "N"])
        self.intronic_set = set(["I"])
        self.uniq = set()

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

            for f in gf:
                if f in self.exonic_set:
                    channels.add("exonic_read")
                    if uniq:
                        channels.add("exonic_UMI")
                elif f in self.intronic_set:
                    channels.add("intronic_read")
                    if uniq:
                        channels.add("intronic_UMI")

        return gene, cell, umi, channels


def parse_cmdline():
    parser = argparse.ArgumentParser(
        description="quantify per-cell gene expression from BAM files by counting into a (sparse) DGE matrix"
    )

    parser.add_argument(
        "bam_in",
        help="bam input (default=stdin)",
        default="/dev/stdin",
        # nargs="+",
    )
    parser.add_argument(
        "--sample-name",
        help="sample identifier (default='NA')",
        default="NA",
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
        default="exonic_UMI",
    )

    parser.add_argument(
        "--output",
        help="directory to store the output h5ad and statistics/marginal counts",
        default="dge",
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
    return np.array(X.sum(axis=axis))


def main(args):
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger("spacemake.quant.main")
    # prepare input and output
    bam = pysam.AlignmentFile(args.bam_in, "rb", check_sq=False, threads=4)
    util.ensure_path(args.output + "/")

    # prepare counter instance
    count_kw = util.load_yaml(args.count_class_params)
    count_cls = get_counter_class(args.count_class)
    counter = count_cls(bam, **count_kw)

    # prepare the sparse-matrix data collection (for multiple channels in parallel)
    dge = DGE(channels=counter.channels)
    # iterate over BAM and count
    for aln in util.timed_loop(bam.fetch(until_eof=True), logger, skim=args.skim):
        if aln.is_unmapped:
            continue

        tags = dict(aln.get_tags())
        gene, cell, umi, channels = counter.parse_bam_record(aln, tags)
        dge.add_read(gene=gene, cell=cell, channels=channels)

    sparse_d, obs, var = dge.make_sparse_arrays()
    adata = dge.sparse_arrays_to_adata(
        sparse_d, obs, var, main_channel=args.main_channel
    )
    adata.write(os.path.join(args.output, f"{args.sample_name}.h5ad"))

    # write out marginal counts across both axes:
    #   summing over obs/cells -> pseudo bulk
    #   summing over genes -> UMI/read distribution
    pname = os.path.join(args.output, f"{args.sample_name}.pseudo_bulk.tsv")

    data = OrderedDict()

    data["sample_name"] = args.sample_name
    data["gene"] = var
    for channel, M in sparse_d.items():
        data[channel] = sparse_summation(M)

    import pandas as pd

    df = pd.DataFrame(data)
    df.to_csv(pname, sep="\t")

    tname = os.path.join(args.output, f"{args.sample_name}.totals.tsv")

    data = OrderedDict()

    data["sample_name"] = args.sample_name
    data["cell"] = obs
    for channel, M in sparse_d.items():
        data[channel] = sparse_summation(M, axis=1)

    import pandas as pd

    df = pd.DataFrame(data)
    df.to_csv(tname, sep="\t")


if __name__ == "__main__":
    args = parse_cmdline()
    main(args)
