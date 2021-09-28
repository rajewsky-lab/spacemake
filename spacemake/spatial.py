from spacemake.preprocess import calculate_adata_metrics

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
    import numpy as np
    from scipy.spatial.distance import cdist

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
    import numpy as np
    import itertools

    # this is hard coded for now for visium, to have 6 neighbors per spot
    # TODO: define an iterative approach where the key is to have around 6 
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
    import numpy as np

    islands = compute_islands(adata, min_umi)

    # find the sizes of the islands. remove the biggest, as either the tissue has a big hole in it
    # or there are not so many big islands in which case removal is OK.
    # to be evaluated later...
    island_sizes = [len(island) for island in islands]

    tissue_islands = np.delete(islands, np.argmax(island_sizes))

    # get the indices of the islands
    tissue_indices = np.where(np.array(adata.obs['total_counts']) >= min_umi)[0]

    tissue_indices = np.append(tissue_indices, np.hstack(tissue_islands))

    adata = adata[tissue_indices, :]

    return adata

def create_meshed_adata(adata,
        width_um,
        spot_diameter_um = 55,
        spot_distance_um = 100,
        bead_diameter_um = 10
    ):
    import pandas as pd
    import scanpy as sc
    import numpy as np
    import anndata

    from sklearn.metrics.pairwise import euclidean_distances
    from scipy.sparse import csr_matrix, csc_matrix, vstack

    coords = adata.obsm['spatial']

    top_left_corner = np.min(coords, axis=0)
    bottom_right_corner = np.max(coords, axis=0)
    width_px = abs(top_left_corner[0] - bottom_right_corner[0])
    height_px = abs(top_left_corner[1] - bottom_right_corner[1])

    # capture area width
    spot_radius_um = spot_diameter_um / 2

    um_by_px = width_um/width_px
    spot_diameter_px = spot_diameter_um / um_by_px
    spot_radius_px = spot_radius_um / um_by_px
    spot_distance_px = spot_distance_um / um_by_px

    height_um = height_px * um_by_px

    def create_mesh(
        width,
        height,
        diameter,
        distance,
        push_x = 0,
        push_y = 0):
        distance_y = np.sqrt(3) * distance
        
        x_coord = np.arange(push_x, width + diameter, distance)
        y_coord = np.arange(push_y, height + diameter, distance_y)
        
        X, Y = np.meshgrid(x_coord, y_coord)
        xy = np.vstack((X.flatten(), Y.flatten())).T
        
        return xy


    # create meshgrid with one radius push
    xy = create_mesh(width_um,
        height_um,
        spot_diameter_um,
        spot_distance_um,
        push_x = spot_radius_um,
        push_y = spot_radius_um)
    # create pushed meshgrid, pushed by 
    xy_pushed = create_mesh(width_um,
        height_um,
        spot_diameter_um,
        spot_distance_um,
        push_x = spot_radius_um + spot_distance_um / 2,
        push_y = spot_radius_um + np.sqrt(3) * spot_distance_um / 2)

    mesh = np.vstack((xy, xy_pushed))
    mesh_px = mesh/um_by_px
    # add the top_left_corner
    mesh_px = mesh_px + top_left_corner

    # example: the diameter of one visium spot is 55um. if we take any bead, which center
    # falls into that, assuming that a bead is 10um, we would in fact take all beads
    # within 65um diameter. for this reason, the max distance should be 45um/2 between
    # visium spot center and other beads
    max_distance_px = (spot_diameter_um - bead_diameter_um)/um_by_px/2

    distance_M = euclidean_distances(mesh_px, coords)

    # new_ilocs contains the indices of the columns of the sparse matrix to be created
    # original_ilocs contains the column location of the original adata.X (csr_matrix)
    new_ilocs, original_ilocs = np.nonzero(distance_M < max_distance_px)

    joined_C = adata.X[original_ilocs]

    # at which indices does the index in the newly created matrix change
    change_ix = np.where(new_ilocs[:-1] != new_ilocs[1:])[0] + 1

    #array of indices, split by which row they should go together

    ix_array = np.asarray(np.split(np.arange(new_ilocs.shape[0]), change_ix, axis=0), dtype='object')

    joined_C_sumed = vstack([csr_matrix(joined_C[ix_array[n].astype(int), :].sum(0)) for n in range(len(ix_array))])

    joined_coordinates = mesh_px[np.unique(new_ilocs)]

    adata_out = anndata.AnnData(csc_matrix(joined_C_sumed), 
        obs = pd.DataFrame({'x_pos': joined_coordinates[:, 0],
                            'y_pos': joined_coordinates[:, 1]}),
        var = adata.var)

    adata_out.obsm['spatial'] = joined_coordinates
    
    # rename index
    adata_out.obs.index.name = 'cell_bc'

    def summarise_adata_obs_column(adata, column, summary_fun=sum):
        vals_to_join = adata.obs[column].to_numpy()[original_ilocs]
        vals_joined = np.array(
            [summary_fun(vals_to_join[ix_array[n].astype(int)])
                for n in range(len(ix_array))])
        return vals_joined
    print(adata)

    # summarise and attach n_reads, calculate metrics (incl. pcr)
    adata_out = calculate_adata_metrics(adata_out,
        # provide the n_reads as a parameter
        n_reads = summarise_adata_obs_column(adata, 'n_reads'))

    adata_out.obs['n_joined'] = [len(x) for x in ix_array]

    from statistics import mean

    for column in ['exact_entropy', 'theoretical_entropy', 'exact_compression',\
        'theoretical_compression']:
        adata_out.obs[column] = summarise_adata_obs_column(adata, column, mean)

    return adata_out

def run_novosparc(dataset, num_spatial_locations=5000, num_input_cells=30000):
    import numpy as np
    import pandas as pd
    import scanpy as sc

    import novosparc
    import anndata

    from scipy.sparse import issparse, csc_matrix

    dataset = dataset[dataset.obs.total_counts, :]

    gene_names = dataset.var.index.tolist()

    num_cells, num_genes = dataset.shape

    if num_cells > num_input_cells:
        sc.pp.subsample(dataset, n_obs = num_input_cells)

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

    locations_circle = novosparc.gm.construct_circle(num_locations = num_spatial_locations)

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


