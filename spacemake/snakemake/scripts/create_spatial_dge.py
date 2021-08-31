import numpy as np
import pandas as pd
import scanpy as sc
from spacemake.util import detect_tissue

dge_path = snakemake.input['dge']

# umi cutoff 
umi_cutoff = int(snakemake.wildcards['umi_cutoff'])

adata = sc.read_h5ad(dge_path)

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

adata.write(snakemake.output[0])
