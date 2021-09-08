import numpy as np
import pandas as pd
import scanpy as sc
from spacemake.util import detect_tissue, attach_barcode_file

dge_path = snakemake.input['dge']

# umi cutoff 
umi_cutoff = int(snakemake.wildcards['umi_cutoff'])

adata = sc.read_h5ad(dge_path)

print('data read')

has_barcode_file = 'barcode_file' in snakemake.input.keys()

# ATTACH BARCODE FILE #
if has_barcode_file:
    adata = attach_barcode_file(adata, snakemake.input['barcode_file'])

# filter out cells based on umi, and genes based on number of cells
sc.pp.filter_cells(adata, min_genes=1)
sc.pp.filter_genes(adata, min_cells=3)

print('data filtered')

# DETECT TISSUE #
# if there is no barcode file, filter adata based on UMI, otherwise detect tissue with UMI cutoff
if has_barcode_file and snakemake.params['downstream_variables']['detect_tissue']:
    tissue_indices = detect_tissue(adata, umi_cutoff)
    adata = adata[tissue_indices, :]
else:
    adata = adata[adata.obs.total_counts > umi_cutoff, :]

adata.write(snakemake.output[0])
