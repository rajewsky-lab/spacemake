def ensure_path(path):
    import os

    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path


def read_fq(fname):
    import gzip
    from more_itertools import grouper

    if fname.endswith(".gz"):
        src = gzip.open(fname, mode="rt")
    elif type(fname) is str:
        src = open(fname)
    else:
        src = fname  # assume its a stream or file-like object already

    for name, seq, _, qual in grouper(src, 4):
        yield name.rstrip()[1:], seq.rstrip(), qual.rstrip()

def calculate_adata_metrics(adata, dge_summary_path=None, n_reads=None):
    import scanpy as sc
    import pandas as pd
    # calculate mitochondrial gene percentage
    adata.var['mt'] = adata.var_names.str.startswith('Mt-') |\
            adata.var_names.str.startswith('mt-') |\
            adata.var_names.str.startswith('MT-')

    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    add_reads = False
    if dge_summary_path is not None:
        dge_summary = pd.read_csv(dge_summary_path, skiprows=7, sep ='\t', index_col = 'cell_bc', names = ['cell_bc', 'n_reads', 'n_umi', 'n_genes'])

        adata.obs = pd.merge(adata.obs, dge_summary[['n_reads']], left_index=True, right_index=True)

        add_reads = True

    if n_reads is not None:
        adata.obs['n_reads'] = n_reads
        add_reads = True

    if add_reads:
        adata.obs['reads_per_counts'] = adata.obs.n_reads / adata.obs.total_counts

    return adata

def dge_to_sparse_adata(dge_path, dge_summary_path):
    import anndata
    import numpy as np
    import gzip
    import pandas as pd
    from scipy.sparse import csr_matrix, vstack

    ix = 0
    mix = 0
    matrices = []
    gene_names = []

    with gzip.open(dge_path, 'rt') as dge:
        first_line = dge.readline().strip().split('\t')

        barcodes = first_line[1:]
        ncol = len(barcodes)
        M = np.zeros((1000, ncol))

        for line in dge:
            vals = line.strip().split('\t')
            # first element contains the gene name
            gene_names.append(vals[0])
            vals = vals[1:]
            vals = np.array(vals, dtype=np.int64)

            M[ix] = vals

            if ix % 1000 == 999:
                mix = mix + 1
                matrices.append(csr_matrix(M))
                ix = 0
                M = np.zeros((1000, ncol))
            else:
                ix = ix + 1

        # get the leftovers
        M = M[:ix]
        matrices.append(csr_matrix(M))

        # sparse expression matrix
        X = vstack(matrices, format='csr')
        
        adata = anndata.AnnData(X.T, obs = pd.DataFrame(index=barcodes), var = pd.DataFrame(index=gene_names))

        adata.obs.index.name = 'cell_bc'

        adata = calculate_adata_metrics(adata, dge_summary_path)

        return adata

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

    return tissue_indices


def attach_barcode_file(adata, barcode_file):
    import pandas as pd

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
    from scipy.sparse import csr_matrix
    from scipy.sparse import vstack

    from spacemake.util import attach_barcode_file

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

    # the diameter of one visium spot is 55um. if we take any bead, which center
    # falls into that, assuming that a bead is 10um, we would in fact take all beads
    # within 65um diameter. for this reason, the max distance should be 45um/2 between
    # visium spot center and other beads
    max_distance_px = 45/um_by_px/2

    distance_M = euclidean_distances(mesh_px, coords)

    # new_ilocs contains the indices of the columns of the sparse matrix to be created
    # original_ilocs contains the column location of the original adata.X (csr_matrix)
    new_ilocs, original_ilocs = np.nonzero(distance_M < max_distance_px)

    joined_C = adata.X[original_ilocs]

    # at which indices does the index in the newly created matrix change
    change_ix = np.where(new_ilocs[:-1] != new_ilocs[1:])[0] + 1

    #array of indices, split by which row they should go together

    ix_array = np.asarray(np.split(np.arange(new_ilocs.shape[0]), change_ix, axis=0), dtype='object')

    joined_C_sumed = vstack([csr_matrix(joined_C[ix_array[n], :].sum(0)) for n in range(len(ix_array))])

    joined_coordinates = mesh_px[np.unique(new_ilocs)]

    adata_out = anndata.AnnData(joined_C_sumed,
        obs = pd.DataFrame({'x_pos': joined_coordinates[:, 0],
                            'y_pos': joined_coordinates[:, 1]}),
        var = adata.var)

    adata_out.obsm['spatial'] = joined_coordinates
    
    # rename index
    adata_out.obs.index.name = 'cell_bc'

    # summarise and attach n_reads, calculate metrics (incl. pcr)
    n_reads = adata.obs.n_reads.to_numpy()[original_ilocs]
    joined_n_reads = np.array([sum(n_reads[ix_array[n]]) for n in range(len(ix_array))])

    adata_out = calculate_adata_metrics(adata_out, n_reads = joined_n_reads)
    adata_out.obs['n_joined'] = [len(x) for x in ix_array]

    return adata_out
