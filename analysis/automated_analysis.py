import numpy as np
import pandas as pd
import scanpy as sc

dge_path = snakemake.input['dge']

# umi cutoff 
umi_cutoff = int(snakemake.wildcards['umi_cutoff'])

adata = sc.read_text(dge_path, delimiter='\t').T

# filter out cells based on umi, and genes based on number of cells
sc.pp.filter_cells(adata, min_counts=umi_cutoff)
sc.pp.filter_cells(adata, min_genes=1)
sc.pp.filter_cells(adata, max_genes=4000)
sc.pp.filter_genes(adata, min_cells=10)

# calculate mitochondrial gene percentage
adata.var['mt'] = adata.var_names.str.startswith('Mt-') | adata.var_names.str.startswith('mt-')

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# save the log-normalised counts in the .raw attribute
adata.raw = adata

# identify highly variable genes if we have any observations
nrow, ncol = adata.shape

if nrow > 1 and ncol > 1:
    sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=5000)

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
    # for a very low number of cells, scipy will throw an error here
    try:
        sc.tl.umap(adata)   
    except TypeError:
        pass
    
    # find out the clusters
    # restrict to max 20 clusters
    resolution = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
    
    top_10_marker_dfs = []
    
    for res in resolution:
        res_key = 'leiden_' + str(res)
        
        sc.tl.leiden(adata, resolution = res, key_added = res_key)
        
        # finding marker genes
        sc.tl.rank_genes_groups(adata, res_key, method='t-test', key_added = 'rank_genes_groups_' + res_key, pts=True)
        
        #save top10 cluster markers
        df = pd.DataFrame(adata.uns['rank_genes_groups_'+res_key]['names'])\
                .head(10).melt(var_name = 'cluster', value_name = 'gene')

        for key in ['logfoldchanges', 'pvals', 'pvals_adj', 'pts', 'pts_rest']:
            df_key = pd.DataFrame(adata.uns['rank_genes_groups_'+res_key][key])\
                .head(10).melt(var_name = 'cluster', value_name = key)
            df[key] = df_key[key]
             
        df['resolution'] = res
        
        top_10_marker_dfs.append(df)

pd.concat(top_10_marker_dfs).to_csv(snakemake.output['cluster_markers'], index=False, sep = '\t')

adata.write(snakemake.output['res_file'])


# save obs_df with umap coordinates
obs_df = sc.get.obs_df(adata, obsm_keys=[('X_umap', 0), ('X_umap', 1)])\
        .join(adata.obs)\
        .rename(columns={'X_umap-0':'umap_0', 'X_umap-1':'umap_1'})

obs_df.index.set_names('cell_bc', inplace=True)

obs_df.to_csv(snakemake.output['obs_df'], sep = '\t')

# save expression values as a long_df
def create_long_df(expr_matrix, id_vars = ['cell_bc']):
    long_df = expr_matrix.melt(id_vars = id_vars, var_name = 'gene', value_name = 'expr')

    long_df = long_df[long_df.expr > 0]
    
    return long_df

# create expr_matrix from raw

expr_matrix = pd.DataFrame(adata.raw.X)
expr_matrix.columns = adata.raw.var.index
expr_matrix['cell_bc'] = adata.obs.index

create_long_df(expr_matrix).to_csv(snakemake.output['long_expr_df'], index=False, sep = '\t')
