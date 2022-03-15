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

def calculate_shannon_entropy_scompression(adata):
    import math
    import itertools
    import numpy as np
    from collections import Counter

    def compute_shannon_entropy(barcode):
        prob, length = Counter(barcode), float(len(barcode))
        return -sum( count/length * math.log(count/length, 2) for count in prob.values())

    def compute_string_compression(barcode):
        compressed_barcode = ''.join(
                letter + str(len(list(group)))
                for letter, group in itertools.groupby(barcode))

        return len(compressed_barcode)

    bc = adata.obs.index.to_numpy()
    bc_len = len(bc[0])
    theoretical_barcodes = np.random.choice(['A','C', 'T', 'G'],
        size = (bc.shape[0], bc_len))

    adata.obs['exact_entropy'] = np.round(np.array(
        [compute_shannon_entropy(cell_bc) for cell_bc in bc]), 2)
    adata.obs['theoretical_entropy'] = np.round(np.array(
        [compute_shannon_entropy(cell_bc) for cell_bc in theoretical_barcodes]), 2)
    adata.obs['exact_compression'] = np.round(np.array(
        [compute_string_compression(cell_bc) for cell_bc in bc]), 2)
    adata.obs['theoretical_compression'] = np.round(np.array(
        [compute_string_compression(cell_bc) for cell_bc in theoretical_barcodes]), 2)

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

        # create an intermediate matrix, which will contain 1000 rows
        # and cell_number columns. we will parse the dge file and fill
        # this matrix line by line. Like this we create gene_number/1000 number
        # of CSR (compressed sparse row) matrices which we join at the end
        M = np.zeros((1000, ncol))

        # read DGE line by line
        # first row: contains CELL BARCODEs
        # each next row contains one gene name, and the counts of that gene
        for line in dge:
            vals = line.strip().split('\t')
            # first element contains the gene name
            # append gene name to gene_names
            gene_names.append(vals[0])
            vals = vals[1:]

            # store counts as np.array
            vals = np.array(vals, dtype=np.int64)

            # update the 1000xcell_number matrix 
            M[ix] = vals

            if ix % 1000 == 999:
                # if we reached the end of M, make it sparse and append to other
                # already read matrices
                mix = mix + 1
                matrices.append(csr_matrix(M))
                ix = 0

                # reset M
                M = np.zeros((1000, ncol))
            else:
                ix = ix + 1

        # get the leftovers: these are the overhang lines, when gene_number is
        # not divisible by 1000
        M = M[:ix]
        matrices.append(csr_matrix(M))

        # sparse expression matrix
        X = vstack(matrices, format='csr')
        
        # create anndata object, but we get the transpose of X, so matrix will
        # be in CSC format
        adata = anndata.AnnData(X.T,
                obs = pd.DataFrame(index=barcodes),
                var = pd.DataFrame(index=gene_names))

        # name the index
        adata.obs.index.name = 'cell_bc'

        # attach metrics such as: total_counts, pct_mt_counts, etc
        # also attach n_genes, and calculate pcr
        calculate_adata_metrics(adata, dge_summary_path)

        # calculate per shannon_entropy and string_compression per bead
        calculate_shannon_entropy_scompression(adata)

        return adata

def load_external_dge(dge_path):
    import anndata
    import scanpy as sc

    from scanpy._utils import check_nonnegative_integers
    from scipy.sparse import issparse, csc_matrix
    from spacemake.errors import SpacemakeError

    adata = sc.read(dge_path)

    if not check_nonnegative_integers(adata.X):
        raise SpacemakeError(f'External dge seems to contain values '+
            'which are already normalised. Raw-count matrix expected.')

    if not issparse(adata.X):
        adata.X = csc_matrix(adata.X)

    # name the index
    adata.obs.index.name = 'cell_bc'

    # attach metrics such as: total_counts, pct_mt_counts, etc
    # also attach n_genes, and calculate pcr
    calculate_adata_metrics(adata)

    return adata

def parse_barcode_file(barcode_file):
    import pandas as pd

    bc = pd.read_csv(barcode_file, sep='[,|\t]', engine='python')

    # rename columns
    bc = bc.rename(columns={'xcoord':'x_pos', 'ycoord':'y_pos','barcodes':'cell_bc', 'barcode':'cell_bc'})\
        .set_index('cell_bc')\
        .loc[:, ['x_pos', 'y_pos']]

    bc = bc.loc[~bc.index.duplicated(keep='first')]

    bc = bc.loc[~bc.index.duplicated(keep='first')]

    return bc

def attach_barcode_file(adata, barcode_file):
    import pandas as pd

    bc = parse_barcode_file(barcode_file)

    # new obs has only the indices of the exact barcode matches
    new_obs = adata.obs.merge(bc, left_index=True, right_index=True, how='inner')
    adata = adata[new_obs.index, :]
    adata.obs = new_obs
    adata.obsm['spatial'] = adata.obs[['x_pos', 'y_pos']].to_numpy()

    return adata

def attach_puck_variables(adata,puck_variables):
    if 'spatial' not in adata.obsm.keys():
        raise SpacemakeError(f'this dataset has no spatial information '+
            'available. Please attach the spatial information using the ' +
            'spacemake.preprocess.attach_barcode_file() function first')

    adata.uns['puck_variables'] = puck_variables

    x_pos_max, y_pos_max = tuple(adata.obsm['spatial'].max(axis=0))

    width_um = adata.uns['puck_variables']['width_um']
    coord_by_um = x_pos_max / width_um
    height_um = int(y_pos_max / coord_by_um)
    
    adata.uns['puck_variables']['height_um'] = height_um
    adata.uns['puck_variables']['coord_by_um'] = coord_by_um

    return adata
