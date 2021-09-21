import numpy as np
import pandas as pd
import scanpy as sc

import novosparc
import anndata

from scipy.sparse import issparse, csc_matrix

data_dir = '/data/local/rajewsky/home/tsztank/projects/spatial/repos/sts-paper/'
data_dir += 'projects/slideseq_v2/processed_data/slideseq_v2_mouse_hippo_new/'
data_dir += 'illumina/complete_data/automated_analysis/slideseq/umi_cutoff_50'

dataset = sc.read(f'{data_dir}/results.h5ad')

dataset = dataset[dataset.obs.total_counts, :]

gene_names = dataset.var.index.tolist()

num_cells, num_genes = dataset.shape

is_var_gene = dataset.var['highly_variable']
# select only 100 genes
var_genes = list(is_var_gene.index[is_var_gene])

dge_rep = dataset.to_df()[var_genes]

# if you uncomment this, dataset.X will be stored not as sparse matrix, rather
# as a regular numpy array. uncommenting this doesnt throw an error
if issparse(dataset.X):
    dense_dataset = anndata.AnnData(
        dataset.X.toarray(),
        obs = dataset.obs,
        var = dataset.var)
    dataset = dense_dataset
    del dense_dataset

num_locations = 5000
locations_circle = novosparc.gm.construct_circle(num_locations = num_locations)

tissue = novosparc.cm.Tissue(dataset=dataset, locations=locations_circle)

num_neighbors_s = num_neighbors_t = 5

tissue.setup_smooth_costs(dge_rep=dge_rep, num_neighbors_s=num_neighbors_s, num_neighbors_t=num_neighbors_t)

tissue.reconstruct(alpha_linear=0, epsilon=5e-3)

dataset_reconst = anndata.AnnData(
    csc_matrix(tissue.sdge.T),
    var = pd.DataFrame(index=gene_names))

# copy of a yet-to-be-pushed novosparc function
def quantify_clusters_spatially(tissue, cluster_key='clusters'):
    """Maps the annotated clusters obtained from the scRNA-seq analysis onto
    the tissue space.
    Args:
        tissue: the novosparc tissue object containing the gene expression data,
                the clusters annotation and the spatial reconstruction. Assumes
                that the cluster annotation exists in the underlying anndata object.
    Returns:
        [numpy array]: An array of the cluster annotation per tissue position.
    """
    clusters = tissue.dataset.obs[cluster_key].to_numpy().flatten()
    return np.array([np.argmax(np.array([np.median(np.array(tissue.gw[:, location][np.argwhere(clusters == cluster).flatten()]))
                                         for cluster in np.unique(clusters)])) for location in range(len(tissue.locations))])
