import logging
import novosparc

from anndata import AnnData
from spacemake.util import message_aggregation
from numpy import ndarray

logger_name = 'spacemake.spatial.novosparc_reconstruction'
logger = logging.getLogger(logger_name)

def get_parser():
    import argparse

    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description='reconstruct 2-d tissue with novosparc',
    )

    parser.add_argument(
        '--single_cell_dataset', type=str, help='path to single_cell.h5ad file',
        required = True,
    )

    parser.add_argument(
        '--spatial_dataset', type=str, help='path to spatial.h5ad file',
        required=False,
    )

    parser.add_argument(
        '--output', type=str, help='path to output.h5ad file',
        required=True,
    )

    return parser

# copy of a yet-to-be-pushed novosparc function
def quantify_clusters_spatially(tissue: novosparc.cm.Tissue, cluster_key: str='clusters') -> ndarray:
    """Maps the annotated clusters obtained from the scRNA-seq analysis onto
    the tissue space.

    :param tissue:    the novosparc tissue object containing the gene expression data,
                      the clusters annotation and the spatial reconstruction. Assumes
                      that the cluster annotation exists in the underlying anndata object.
    :type tissue: novosparc.cm.Tissue
    :param cluster_key: the key of the clustering to be used for deduction
    :type cluster_key: str
    rtype: numpy.ndarray
    """
    import numpy as np

    clusters = tissue.dataset.obs[cluster_key].to_numpy().flatten()
    cluster_names = np.unique(clusters)
    ixs = np.array(
            [np.argmax(
                np.array(
                    [np.median(
                        np.array(
                            tissue.gw[:, location][np.argwhere(clusters == cluster).flatten()]))
                                         for cluster in cluster_names]))
                    for location in range(len(tissue.locations))])
    return [cluster_names[ix] for ix in ixs]
    
def novosparc_denovo(
        adata: AnnData,
        num_spatial_locations: int=5000,
        num_input_cells: int=30000,
        locations=None,
    ) -> novosparc.cm.Tissue:
    """reconstruct spatial information de-novo with novosparc.

    :param adata: A spacemake processed sample.
    :type adata: AnnData
    :param num_spatial_locations:
    :type num_spatial_locations: int
    :param num_input_cells:
    :type num_input_cells: int
    :param locations: Numpy array of (x,y) locations (optional)
    :type locations: numpy.ndarray
    :rtype: novosparc.cm.Tissue
    """
    import numpy as np
    import pandas as pd
    import scanpy as sc

    import novosparc
    import anndata

    from scipy.sparse import issparse, csc_matrix

    logger.info('Reconstructing the tissue de-novo with novosparc') 
    
    gene_names = adata.var.index.tolist()

    num_cells, num_genes = adata.shape

    if num_cells > num_input_cells:
        sc.pp.subsample(adata, n_obs = num_input_cells)

    if num_cells < num_spatial_locations:
        num_spatial_locations = num_cells

    sc.pp.highly_variable_genes(adata, n_top_genes = 100)
    is_var_gene = adata.var['highly_variable']
    # select only 100 genes
    var_genes = list(is_var_gene.index[is_var_gene])

    dge_rep = adata.to_df()[var_genes]

    # if you uncomment this, adata.X will be stored not as sparse matrix, rather
    # as a regular numpy array. uncommenting this doesnt throw an error
    if issparse(adata.X):
        dense_adata = anndata.AnnData(
            adata.X.toarray(),
            obs = adata.obs,
            var = adata.var)
        adata = dense_adata
        del dense_adata

    if locations is None:
        # create circle locations
        locations = novosparc.gm.construct_circle(num_locations = num_spatial_locations)

    tissue = novosparc.cm.Tissue(dataset=adata, locations=locations)

    num_neighbors_s = num_neighbors_t = 5

    logger.info('Novosparc setup')
    tissue.setup_smooth_costs(dge_rep=dge_rep, num_neighbors_s=num_neighbors_s, num_neighbors_t=num_neighbors_t)

    tissue.reconstruct(alpha_linear=0, epsilon=5e-3)

    return tissue

def novosparc_mapping(sc_adata: AnnData, st_adata: AnnData) -> novosparc.cm.Tissue:
    """novosparc_mapping.

    :param sc_adata:
    :type sc_adata: AnnData
    :param st_adata:
    :type st_adata: AnnData
    :rtype: novosparc.cm.Tissue
    """
    import scanpy as sc
    import anndata
    import novosparc

    from scanpy._utils import check_nonnegative_integers
    from scipy.sparse import csc_matrix

    logger.info('Mapping single-cell data onto spatial data with novosparc')

    if (check_nonnegative_integers(sc_adata.X)
        or check_nonnegative_integers(st_adata.X)
    ):
        # if any of the inputs is count-data, raise error
        raise SpacemakeError(f'External dge seems to contain raw counts. '+
            'Normalised values are expected for both sc_adata and st_adata.')

    # calculate top 500 variable genes for both
    sc.pp.highly_variable_genes(sc_adata, n_top_genes=500)
    sc.pp.highly_variable_genes(st_adata, n_top_genes=500)

    sc_adata_hv = sc_adata.var_names[sc_adata.var.highly_variable].to_list()
    st_adata_hv = st_adata.var_names[st_adata.var.highly_variable].to_list()

    markers = list(set(sc_adata_hv).intersection(st_adata_hv))

    logger.info(f'{len(markers)} number of common markers found. Using them' +
        ' for reconstruction')

    # save sc dge as a pandas dataframe
    dge_rep = sc_adata.to_df()[sc_adata_hv]

    if not 'spatial' in st_adata.obsm:
        raise SpacemakeError(f'The object provided to st_adata is not spatial')

    locations = st_adata.obsm['spatial']
    atlas_matrix = st_adata.to_df()[markers].values

    # make dense dataset
    dense_dataset = anndata.AnnData(
        sc_adata.X.toarray(),
        obs = sc_adata.obs,
        var = sc_adata.var)

    marker_ix = [dense_dataset.var.index.get_loc(marker) for marker in markers]

    tissue = novosparc.cm.Tissue(dataset=dense_dataset, locations=locations)
    num_neighbors_s = num_neighbors_t = 5

    tissue.setup_linear_cost(markers_to_use=marker_ix, atlas_matrix=atlas_matrix,
                             markers_metric='minkowski', markers_metric_p=2)
    tissue.setup_smooth_costs(dge_rep = dge_rep,
                              num_neighbors_s=num_neighbors_s,
                              num_neighbors_t=num_neighbors_t)

    tissue.reconstruct(alpha_linear=0.5, epsilon=5e-3)

    return tissue

def save_novosparc_res(tissue, adata_original) -> AnnData:
    import anndata
    import pandas as pd
    import numpy as np

    from scipy.sparse import csc_matrix

    adata_reconst = anndata.AnnData(
        csc_matrix(tissue.sdge.T),
        var = pd.DataFrame(index=tissue.dataset.var_names))

    logger.info('Scaling mapped data to original data...')

    adata_reconst.X = np.sum(adata_original.X) * adata_reconst.X / np.sum(adata_reconst.X)

    adata_reconst.obsm['spatial'] = tissue.locations

    logger.info('Transferring original cluster labels...')

    for res_key in adata_original.obs.columns[adata_original.obs.columns.str.startswith('leiden_')]:
        adata_reconst.obs[f'novosparc_{res_key}'] = quantify_clusters_spatially(tissue, res_key)

    return adata_reconst

@message_aggregation(logger_name)
def cmdline():
    import scanpy as sc

    parser = get_parser()

    args = parser.parse_args()

    sc_adata = sc.read(args.single_cell_dataset)

    if args.spatial_dataset is not None:
        # make mapping
        st_adata = sc.read(args.spatial_dataset)

        tissue_reconst = novosparc_mapping(sc_adata, st_adata)
        adata = save_novosparc_res(tissue_reconst, sc_adata)
    else:
        tissue_reconst = novosparc_denovo(sc_adata)  
        adata = save_novosparc_res(tissue_reconst, sc_adata)

    adata.write(args.output)

if __name__ == "__main__":
    cmdline()
