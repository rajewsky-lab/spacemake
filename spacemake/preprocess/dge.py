import logging

logger_name = "spacemake.preprocess.dge"
logger = logging.getLogger(logger_name)

def calculate_adata_metrics(adata, dge_summary_path=None, n_reads=None):
    import scanpy as sc
    import pandas as pd

    # calculate mitochondrial gene percentage
    adata.var["mt"] = (
        adata.var_names.str.startswith("Mt-")
        | adata.var_names.str.startswith("mt-")
        | adata.var_names.str.startswith("MT-")
    )

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    add_reads = False
    if dge_summary_path is not None:
        dge_summary = pd.read_csv(
            dge_summary_path,
            skiprows=7,
            sep="\t",
            index_col="cell_bc",
            names=["cell_bc", "n_reads", "n_umi", "n_genes"],
        )

        adata.obs = pd.merge(
            adata.obs, dge_summary[["n_reads"]], left_index=True, right_index=True
        )

        add_reads = True

    if n_reads is not None:
        adata.obs["n_reads"] = n_reads
        add_reads = True

    if add_reads:
        adata.obs["reads_per_counts"] = adata.obs.n_reads / adata.obs.total_counts


def calculate_shannon_entropy_scompression(adata):
    import math
    import itertools
    import numpy as np
    from collections import Counter

    def compute_shannon_entropy(barcode):
        prob, length = Counter(barcode), float(len(barcode))
        return -sum(
            count / length * math.log(count / length, 2) for count in prob.values()
        )

    def compute_string_compression(barcode):
        compressed_barcode = "".join(
            letter + str(len(list(group)))
            for letter, group in itertools.groupby(barcode)
        )

        return len(compressed_barcode)

    bc = adata.obs.index.to_numpy()
    bc_len = len(bc[0])
    theoretical_barcodes = np.random.choice(
        ["A", "C", "T", "G"], size=(bc.shape[0], bc_len)
    )

    adata.obs["exact_entropy"] = np.round(
        np.array([compute_shannon_entropy(cell_bc) for cell_bc in bc]), 2
    )
    adata.obs["theoretical_entropy"] = np.round(
        np.array(
            [compute_shannon_entropy(cell_bc) for cell_bc in theoretical_barcodes]
        ),
        2,
    )
    adata.obs["exact_compression"] = np.round(
        np.array([compute_string_compression(cell_bc) for cell_bc in bc]), 2
    )
    adata.obs["theoretical_compression"] = np.round(
        np.array(
            [compute_string_compression(cell_bc) for cell_bc in theoretical_barcodes]
        ),
        2,
    )


def dge_to_sparse_adata(dge_path, dge_summary_path):
    import anndata
    import numpy as np
    import gzip
    import pandas as pd
    from scipy.sparse import coo_matrix, hstack

    gene_names = []

    with gzip.open(dge_path, "rt") as dge:
        first_line = dge.readline().strip().split("\t")
        has_mt = False
        barcodes = first_line[1:]
        N_bc = len(barcodes)
        X = None

        # read DGE line by line
        # first row: contains CELL BARCODEs
        # each next row contains one gene name, and the counts of that gene
        for line in dge:
            vals = line.strip()
            _idx_tab = vals.index("\t")
            _gene_name = vals[:_idx_tab]
            gene_names.append(_gene_name)

            if _gene_name.lower().startswith("mt-"):
                has_mt = True
                
            # store counts as np.array
            _vals = np.fromstring(vals[_idx_tab:], dtype=np.int32, count=N_bc, sep='\t').flatten()
            _idx_nonzero = np.argwhere(_vals != 0).flatten()

            if len(_idx_nonzero) > 0:
                gene_sp = coo_matrix((_vals[_idx_nonzero].astype(np.int32), (_idx_nonzero, np.zeros(len(_idx_nonzero)))), shape=(N_bc, 1), dtype=np.int32)
            else:
                gene_sp = coo_matrix((N_bc, 1), dtype=np.int32)

            if X is None:
                 X = gene_sp
            else:
                 X = hstack([X, gene_sp])

        if X is None:
            X = coo_matrix((len(barcodes), 0), dtype=np.int32)
    
        if not has_mt:
            # ensure we have an entry for mitochondrial transcripts even if it's just all zeros
            print(
                "need to add mt-missing because no mitochondrial stuff was among the genes for annotation"
            )
            gene_names.append("mt-missing")
            X = hstack([X, np.zeros(X.shape[0])[:, None]])

        X = X.tocsr()
        X = X.astype(np.float32)
        adata = anndata.AnnData(
            X, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=gene_names)
        )

        # name the index
        adata.obs.index.name = "cell_bc"

        # attach metrics such as: total_counts, pct_mt_counts, etc
        # also attach n_genes, and calculate pcr
        calculate_adata_metrics(adata, dge_summary_path)

        # calculate per shannon_entropy and string_compression per bead
        calculate_shannon_entropy_scompression(adata)

        if adata.X.sum() == 0:
            logger.warn(f"The DGE from {dge_path} is empty")

        return adata


def load_external_dge(dge_path):
    import scanpy as sc

    from scanpy._utils import check_nonnegative_integers
    from scipy.sparse import issparse, csc_matrix
    from spacemake.errors import SpacemakeError

    adata = sc.read(dge_path)

    if not check_nonnegative_integers(adata.X):
        raise SpacemakeError(
            f"External dge seems to contain values "
            + "which are already normalised. Raw-count matrix expected."
        )

    if not issparse(adata.X):
        adata.X = csc_matrix(adata.X)

    # name the index
    adata.obs.index.name = "cell_bc"

    # attach metrics such as: total_counts, pct_mt_counts, etc
    # also attach n_genes, and calculate pcr
    calculate_adata_metrics(adata)

    return adata


def parse_barcode_file(barcode_file):
    import pandas as pd

    bc = pd.read_csv(barcode_file, sep="[,|\t]", engine='python')

    # rename columns
    bc = (
        bc.rename(
            columns={
                "xcoord": "x_pos",
                "ycoord": "y_pos",
                "barcodes": "cell_bc",
                "barcode": "cell_bc",
            }
        )
        .set_index("cell_bc")
        .loc[:, ["x_pos", "y_pos"]]
    )

    bc = bc.loc[~bc.index.duplicated(keep="first")]

    bc = bc.loc[~bc.index.duplicated(keep="first")]

    return bc


def attach_barcode_file(adata, barcode_file):
    bc = parse_barcode_file(barcode_file)

    # new obs has only the indices of the exact barcode matches
    new_obs = adata.obs.merge(bc, left_index=True, right_index=True, how="inner")
    adata = adata[new_obs.index, :]
    adata.obs = new_obs
    adata.obsm["spatial"] = adata.obs[["x_pos", "y_pos"]].to_numpy()

    return adata


def attach_puck_variables(adata, puck_variables):
    if "spatial" not in adata.obsm.keys():
        raise SpacemakeError(
            f"this dataset has no spatial information "
            + "available. Please attach the spatial information using the "
            + "spacemake.preprocess.attach_barcode_file() function first"
        )

    adata.uns["puck_variables"] = puck_variables

    x_pos_max, y_pos_max = tuple(adata.obsm["spatial"].max(axis=0))
    x_pos_min, y_pos_min = tuple(adata.obsm["spatial"].min(axis=0))
    #print(f"PUCK VARS {puck_variables} X MIN {x_pos_min} X MAX {x_pos_max} Y MIN {y_pos_min} Y MAX {y_pos_max}")

    width_um = adata.uns["puck_variables"]["width_um"]
    coord_by_um = (x_pos_max - x_pos_min) / width_um

    # this can be NaN if only one coordinate (only one cell, will fail)
    if coord_by_um > 0:
        height_um = int((y_pos_max - y_pos_min) / coord_by_um)
    else:
        height_um = 1 # avoid division by zero and error in reports
        coord_by_um = 1

    adata.uns["puck_variables"]["height_um"] = height_um
    adata.uns["puck_variables"]["coord_by_um"] = coord_by_um

    return adata


def attach_puck(adata, puck):
    attach_puck_variables(adata, puck.variables)
    adata.uns["puck_name"] = puck.name

    return adata