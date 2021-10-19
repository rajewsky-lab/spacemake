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
    pd.DataFrame().to_csv(snakemake.output['nhood_enrichment'])
else:
    res_keys = adata.obs.columns[adata.obs.columns.str.startswith('leiden_')]
    
    top_10_marker_dfs = []
    nhood_enrichment_dfs = []
    
    for res_key in res_keys:
        rank_key = 'rank_genes_groups_' + res_key
        
        df = pd.DataFrame(adata.uns[rank_key]['names'])\
                .melt(var_name = 'cluster', value_name = 'gene')
        
        for key in ['logfoldchanges', 'pvals', 'pvals_adj']:
            df_key = pd.DataFrame(adata.uns[rank_key][key])\
                .melt(var_name = 'cluster', value_name = key)
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
             
        df['resolution'] = res_key.split('_')[1]
        df.reset_index(inplace=True)
         
        top_10_marker_dfs.append(df)

        if snakemake.params['is_spatial']:
        # get nhood data
            df = pd.DataFrame(adata.uns[f'{res_key}_nhood_enrichment']['zscore'])
            df = pd.melt(df.reset_index(), id_vars='index')\
                .rename(columns={'index': 'cluster_a',
                                 'variable': 'cluster_b',
                                 'value': 'zscore'})
            df['resolution'] = res_key.split('_')[1]

            nhood_enrichment_dfs.append(df)
        
    pd.concat(top_10_marker_dfs).to_csv(snakemake.output['cluster_markers'], index=False)

    if snakemake.params['is_spatial']:
        pd.concat(nhood_enrichment_dfs).to_csv(snakemake.output['nhood_enrichment'], index=False)
    else:
        # output empty csv file
        pd.DataFrame().to_csv(snakemake.output['nhood_enrichment'])


###############
# SAVE OBS DF #
###############
obs_df = adata.obs

if adata_complete:
    obs_df = sc.get.obs_df(adata, obsm_keys=[('X_umap', 0), ('X_umap', 1)])\
        .join(obs_df)\
        .rename(columns={'X_umap-0':'umap_0', 'X_umap-1':'umap_1'})

obs_df.index.set_names('cell_bc', inplace=True)

obs_df.to_csv(snakemake.output['obs_df'])

###############
# SAVE VAR DF #
###############
adata.var.index.set_names('gene_name', inplace=True)

adata.var.to_csv(snakemake.output['var_df'])
