import numpy as np
import pandas as pd
import scanpy as sc
import itertools
from scipy.spatial.distance import cdist
from scipy.sparse import csr_matrix

#############
# FUNCTIONS #
#############
def compute_neighbors(adata, min_dist=None, max_dist=None):
    '''Compute all direct neighbors of all spots in the adata. Currently
        tailored for 10X Visium.
    Args:
        adata: an AnnData object
        min_dist: int, minimum distance to consider neighbors
        max_dist: int, maximum distance to consider neighbors
    Returns:
        neighbors: a dictionary holding the spot IDs for every spot
    '''
    # Calculate all Euclidean distances on the adata
    dist_mtrx = cdist(adata.obsm['spatial'],
                      adata.obsm['spatial'])
    neighbors = dict()

    for i in range(adata.obs.shape[0]):
        neighbors[i] = np.where((dist_mtrx[i,:] > min_dist) & (dist_mtrx[i,:] < max_dist))[0]

    return neighbors

def compute_islands(adata, min_umi):
    '''Find contiguity islands

    Args:
        adata: an AnnData object
        cluster: a cell type label to compute the islands for
    Returns:
        islands: A list of lists of spots forming contiguity islands
    '''
    # this is hard coded for now for visium, to have 6 neighbors per spot
    # TODO: define an iterative approach where he key is to have around 6 
    # neighbors per spot on average
    neighbors = compute_neighbors(adata, min_dist = 0, max_dist=3)
    spots_cluster = np.where(np.array(adata.obs['total_counts']) < min_umi)[0]

    # create a new dict holding only neighbors of the cluster
    islands = []

    for spot in neighbors:
        if spot not in spots_cluster:
            continue
        islands.append({spot}.union({x for x in neighbors[spot] if x in spots_cluster}))

    # merge islands with common spots
    island_spots = set(itertools.chain.from_iterable(islands)) 

    for each in island_spots:
        components = [x for x in islands if each in x]
        for i in components:
            islands.remove(i)
        islands += [list(set(itertools.chain.from_iterable(components)))]

    return islands

def detect_tissue(adata, min_umi):
    ''' Detect tissue: first find beads with at least min_umi UMIs, then detect island in the rest
    
    Args
        adata: an AnnData object, with spatial coordinates
        min_umi: integer, the min umi to be assigned as tissue bead by default
    Returns:
        tissue_indices: a list of indices which should be kept for this AnnData object
    '''

    islands = compute_islands(adata, min_umi)

    # find the sizes of the islands. remove the biggest, as either the tissue has a big hole in it
    # or there are not so many big islands in which case removal is OK.
    # to be evaluated later...
    island_sizes = [len(island) for island in islands]

    tissue_islands = np.delete(islands, np.argmax(island_sizes))

    # get the indices of the islands
    tissue_indices = np.where(np.array(adata.obs['total_counts']) >= min_umi)[0]

    tissue_indices = np.append(tissue_indices, np.hstack(tissue_islands))

    return tissue_indices


########
# CODE #
########
dge_path = snakemake.input['dge']

# umi cutoff 
umi_cutoff = int(snakemake.wildcards['umi_cutoff'])

adata = sc.read_text(dge_path, delimiter='\t').T
adata.X = csr_matrix(adata.X)

print('data read')

has_barcode_file = 'barcode_file' in snakemake.input.keys()

# ATTACH BARCODE FILE #
if has_barcode_file:
    barcode_file = snakemake.input['barcode_file']
    bc = pd.read_csv(barcode_file, sep='[,|\t]', engine='python')

    # rename columns
    bc = bc.rename(columns={'xcoord':'x_pos', 'ycoord':'y_pos','barcodes':'cell_bc', 'barcode':'cell_bc'})\
        .set_index('cell_bc')\
        .loc[:, ['x_pos', 'y_pos']]

    bc = bc.loc[~bc.index.duplicated(keep='first')]

    bc = bc.loc[~bc.index.duplicated(keep='first')]
    # new obs has only the indices of the exact barcode matches
    new_obs = adata.obs.merge(bc, left_index=True, right_index=True, how='inner')
    adata = adata[new_obs.index, :]
    adata.obs = new_obs
    adata.obsm['spatial'] = adata.obs[['x_pos', 'y_pos']].to_numpy()


# filter out cells based on umi, and genes based on number of cells
sc.pp.filter_cells(adata, min_genes=1)
sc.pp.filter_genes(adata, min_cells=3)

print('data filtered')

# calculate mitochondrial gene percentage
adata.var['mt'] = adata.var_names.str.startswith('Mt-') |\
        adata.var_names.str.startswith('mt-') |\
        adata.var_names.str.startswith('MT-')

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# DETECT TISSUE #
# if there is no barcode file, filter adata based on UMI, otherwise detect tissue with UMI cutoff
if has_barcode_file and snakemake.params['downstream_variables']['detect_tissue']:
    tissue_indices = detect_tissue(adata, umi_cutoff)
    adata = adata[tissue_indices, :]
else:
    adata = adata[adata.obs.total_counts > umi_cutoff, :]

# save the raw counts
adata.raw = adata

# identify highly variable genes if we have any observations
nrow, ncol = adata.shape

# require at least 1000 genes expressed in the sample
if nrow > 1 and ncol >= 1000:
    print('starting analysis')
    try:
        sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=2000)
    except ValueError:
        sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=1000, span = 1)
    
    # calculate log(cpm)
    print('normalising and log-scaling')
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata, base=2)
    
    # PCA ANALYSIS
    print('calculating pca components')
    sc.tl.pca(adata, svd_solver='arpack')
    
    # get number of pcs-s identified. Sometimes can be smaller than 50, if
    # less than 50 cells pass the threshold
    n_pcs = adata.uns['pca']['variance'].size
    # limit the number of pcs to 50
    n_pcs = n_pcs if n_pcs < 40 else 40
    
    # Compute the neighborhood graph
    print('computing neighborhood graph')
    sc.pp.neighbors(adata, n_pcs=n_pcs)
    
    # compute UMAP
    # for a very low number of cells, scanpy will throw an error here
    try:
        print('dimensionality reduction')
        sc.tl.umap(adata)   
    except TypeError:
        pass
    
    # find out the clusters
    # restrict to max 20 clusters
    resolution = [0.2, 0.4, 0.6, 0.8, 1]

    print('clustering')
    
    for res in resolution:
        res_key = 'leiden_' + str(res)
        
        sc.tl.leiden(adata, resolution = res, key_added = res_key)
        
        # finding marker genes
        print(f'ranking genes for resolution {res}')
        sc.tl.rank_genes_groups(adata, res_key, method='t-test', key_added = 'rank_genes_groups_' + res_key, pts=True,
            use_raw = False)

adata.write(snakemake.output[0])
