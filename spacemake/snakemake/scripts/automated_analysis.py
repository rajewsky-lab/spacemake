import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

from spacemake.spatial import detect_tissue

# expect fitlered .h5ad dge, with spatial coords attached, tissue detected
adata = sc.read_h5ad(snakemake.input[0])
umi_cutoff = int(snakemake.wildcards['umi_cutoff'])

# filter_umi or detect tissue
# if data spatial and detect_tissue=True
if 'spatial' in adata.obsm.keys() and snakemake.params['run_mode_variables']['detect_tissue']:
    adata = detect_tissue(adata, umi_cutoff)
    print('tissue detection')
else:
    print(f'filtering by umi cutoff: {umi_cutoff}')
    adata = adata[adata.obs.total_counts > umi_cutoff, :]

# make the var indices (gene names) and obs indices (cell barcode) unique
adata.obs_names_make_unique()
adata.var_names_make_unique()

# save the raw counts
adata.raw = adata

# identify highly variable genes if we have any observations
nrow, ncol = adata.shape

# require at least 1000 genes expressed in the sample and at least 100 cells
if nrow > 100 and ncol >= 1000:
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
    resolution = [0.4, 0.6, 0.8, 1.0, 1.2]

    print('clustering')


    if snakemake.params['is_spatial']:
        sq.gr.spatial_neighbors(adata, coord_type="generic")

    for res in resolution:
        try:
            res_key = 'leiden_' + str(res)
            
            sc.tl.leiden(adata, resolution = res, key_added = res_key)
            
            # finding marker genes
            print(f'ranking genes for resolution {res}')
            sc.tl.rank_genes_groups(adata, res_key, method='t-test', key_added = 'rank_genes_groups_' + res_key, pts=True,
                use_raw = False)
            if snakemake.params['is_spatial']:
                # calculate nhood enrichment from squidpy
                sq.gr.nhood_enrichment(adata, cluster_key=res_key)
        except ZeroDivisionError as e:
            pass



adata.write(snakemake.output[0])
