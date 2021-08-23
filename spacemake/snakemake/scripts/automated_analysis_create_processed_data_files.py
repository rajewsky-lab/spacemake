import numpy as np
import pandas as pd
import scanpy as sc

# save expression values as a long_df
def create_long_df(expr_matrix, id_vars = ['cell_bc']):
    long_df = expr_matrix.melt(id_vars = id_vars, var_name = 'gene', value_name = 'expr')
    long_df = long_df[long_df.expr > 0]
    return long_df

##############
# LOAD ADATA #
##############

adata = sc.read(snakemake.input[0])

uns_keys = ['hvg', 'leiden', 'log1p', 'neighbors', 'pca', 'umap']

# all the keys have to be in adata.uns
adata_complete = any([key in adata.uns.keys() for key in uns_keys])
#################
# TOP10 markers #
#################
if not adata_complete:
    pd.DataFrame().to_csv(snakemake.output['cluster_markers'])
else:
    resolution = [0.2, 0.4, 0.6, 0.8, 1]
    
    top_10_marker_dfs = []
    
    for res in resolution:
        res_key = 'leiden_' + str(res)
        rank_key = 'rank_genes_grops_' + res_key
        
        if res_key not in adata.obs.columns or rank_key not in adata.uns.keys():
            sc.tl.leiden(adata, resolution = res, key_added = res_key)
            sc.tl.rank_genes_groups(adata, res_key, method='t-test',
                    key_added = rank_key, pts=True,
                    use_raw = False)
        
        df = pd.DataFrame(adata.uns[rank_key]['names'])\
                .head(10).melt(var_name = 'cluster', value_name = 'gene')
        
        for key in ['logfoldchanges', 'pvals', 'pvals_adj']:
            df_key = pd.DataFrame(adata.uns[rank_key][key])\
                .head(10).melt(var_name = 'cluster', value_name = key)
            df[key] = df_key[key]
            # set the index to gene-cluster pair
            
        df.set_index(['gene', 'cluster'], inplace=True)
            
        for key in ['pts', 'pts_rest']:
            # get the percentage expressed in cluster and rest
            df2 = adata.uns[rank_key][key]
            df2['gene'] = df2.index
            df2 = df2.melt(var_name='cluster', id_vars='gene')\
                .set_index(['gene', 'cluster'])
            
            df[key] = df2.loc[df.index].value
             
        df['resolution'] = res
        df.reset_index(inplace=True)
         
        top_10_marker_dfs.append(df)
        
    pd.concat(top_10_marker_dfs).to_csv(snakemake.output['cluster_markers'], index=False, sep = '\t')

###############
# SAVE OBS DF #
###############
obs_df = adata.obs

if adata_complete:
    obs_df = sc.get.obs_df(adata, obsm_keys=[('X_umap', 0), ('X_umap', 1)])\
        .join(obs_df)\
        .rename(columns={'X_umap-0':'umap_0', 'X_umap-1':'umap_1'})

obs_df.index.set_names('cell_bc', inplace=True)

obs_df.to_csv(snakemake.output['obs_df'], sep = '\t')

###############
# SAVE VAR DF #
###############
adata.var.index.set_names('gene_name', inplace=True)

adata.var.to_csv(snakemake.output['var_df'], sep = '\t')
