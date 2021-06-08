import numpy as np
import pandas as pd
import scanpy as sc

print(snakemake.input)

dge_path = snakemake.input['dge']

# umi cutoff 
umi_cutoff = int(snakemake.wildcards['umi_cutoff'])

adata = sc.read_text(dge_path, delimiter='\t').T

# filter out cells based on umi, and genes based on number of cells
sc.pp.filter_cells(adata, min_counts=umi_cutoff)
sc.pp.filter_cells(adata, min_genes=1)
sc.pp.filter_cells(adata, max_genes=4000)
sc.pp.filter_genes(adata, min_cells=3)

# calculate mitochondrial gene percentage
adata.var['mt'] = adata.var_names.str.startswith('Mt-') | adata.var_names.str.startswith('mt-')

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# save the raw counts
adata.raw = adata

# identify highly variable genes if we have any observations
nrow, ncol = adata.shape

# require at least 1000 genes expressed in the sample
if nrow > 1 and ncol >= 1000:
    try:
        sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=2000)
    except ValueError:
        sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=1000, span = 1)
    
    # calculate log(cpm)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata, base=2)
    
    # PCA ANALYSIS
    sc.tl.pca(adata, svd_solver='arpack')
    
    # get number of pcs-s identified. Sometimes can be smaller than 50, if
    # less than 50 cells pass the threshold
    n_pcs = adata.uns['pca']['variance'].size
    # limit the number of pcs to 50
    n_pcs = n_pcs if n_pcs < 40 else 40
    
    # Compute the neighborhood graph
    sc.pp.neighbors(adata, n_pcs=n_pcs)
    
    # compute UMAP
    # for a very low number of cells, scanpy will throw an error here
    try:
        sc.tl.umap(adata)   
    except TypeError:
        pass
    
    # find out the clusters
    # restrict to max 20 clusters
    resolution = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
    
    for res in resolution:
        res_key = 'leiden_' + str(res)
        
        sc.tl.leiden(adata, resolution = res, key_added = res_key)
        
        # finding marker genes
        sc.tl.rank_genes_groups(adata, res_key, method='t-test', key_added = 'rank_genes_groups_' + res_key, pts=True)

#######################
# ATTACH BARCODE FILE #
#######################
if 'barcode_file' in snakemake.input:
    bc = pd.read_csv(snakemake.input['barcode_file'], sep='[,|\t]', engine='python')

    # rename columns
    bc = bc.rename(columns={'xcoord':'x_pos', 'ycoord':'y_pos','barcodes':'cell_bc', 'barcode':'cell_bc'})\
        .set_index('cell_bc')\
        .loc[:, ['x_pos', 'y_pos']]

    bc = bc.loc[~bc.index.duplicated(keep='first')]

    adata.obs = adata.obs.merge(bc, left_index=True, right_index=True, how='left')
    adata.obsm['spatial'] = adata.obs[['x_pos', 'y_pos']].to_numpy()

adata.write(snakemake.output[0])
