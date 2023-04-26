import pytest
import os
import scanpy as sc
import spacemake.tk as tk

def test_tk():
    fname = os.path.join(os.path.dirname(__file__), "../test_data/mini_adata.h5ad")
    adata = sc.read_h5ad(fname)
    adata = tk.add_common_metrics(adata)
    adata = tk.pre_filter(adata, umi_cutoff=1, normalize=True)
