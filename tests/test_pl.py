import pytest
import numpy as np
import spacemake.pl as pl
import spacemake.tk as tk
import scipy.sparse
import scanpy as sc

def make_test_adata(n_cells=1000, n_genes=2000, n_data=10000, mt_genes=10, n_coding=900):
    
    counts = np.random.randint(1, 300, size=n_data)
    reads = counts + np.random.randint(0, 150, size=n_data)

    i = np.random.randint(0, n_cells, size=n_data)
    j = np.random.randint(0, n_genes, size=n_data)

    X = scipy.sparse.coo_matrix((counts, (i, j))).tocsr()
    R = scipy.sparse.coo_matrix((reads, (i, j))).tocsr()
    import anndata
    adata = anndata.AnnData(X, dtype=np.float32)
    adata.layers["exonic_counts"] = adata.X.copy()
    adata.layers["exonic_reads"] = R
    
    obs_names = [random_str() for i in range(n_cells)]
    var_names = [random_str(L=5, letters="ABCDEFGHIJKLMNOPQRSTUVWXYZ") for i in range(n_genes)]
    
    # make some genes mt
    for i in np.random.randint(0, len(var_names), size=mt_genes):
        var_names[i] = "mt-" + var_names[i]

    # simulate ambient RNA
    obs_names[-1] = 'NA'

    adata.obs_names = obs_names
    adata.obs_names_make_unique()

    # select some genes as protein coding 
    adata.var_names = var_names
    adata.var_names_make_unique()

    coding_genes = set(adata.var_names[i] for i in np.random.randint(0, len(var_names), size=n_coding))
    print(f"coding genes {coding_genes}")
    adata.uns['coding_genes'] = sorted(coding_genes)


    return adata

def random_str(L=10, letters="ACGT"):
    letters = np.array(list(letters))
    i = np.random.randint(0, len(letters), size=L)

    word = "".join(letters[i])
    return word

def test_pl():
    adata = make_test_adata()
    print(adata)
    print(adata.obs_names)
    print(adata.var_names)
    print(adata.X.sum())
    adata = tk.add_common_metrics(adata)

    fig, x_cut = pl.loglog_knee(adata)
    fig2 = pl.dge_stats(adata)
    ratio = pl.reads_per_UMI(adata)
    fig3 = pl.cell_metrics(adata)

def test_miRNA_pl(tmp_path_factory, seed=47110815):
    np.random.seed(seed)
    tmp = tmp_path_factory.mktemp("mirna_pl").as_posix()
    adata =  tk.add_common_metrics(make_test_adata())

    fcod = f"{tmp}/coding_genes.txt"
    with open(fcod, 'w') as f:
        for gene in adata.uns['coding_genes']:
            f.write(gene + '\n') 

    fma = f"{tmp}/adata.h5ad"
    adata.uns = {}
    adata.write(fma)

    midata = tk.add_common_metrics(make_test_adata(n_cells=1200, n_genes=100))
    obs_names = midata.obs_names.to_list()
    obs_names[:len(adata.obs_names) - 1] = adata.obs_names.to_list()[:-1]
    midata.obs_names = obs_names
    midata.obs_names_make_unique()

    midata.var_names = [f"Hsa-Mir-{n+1}-P1-v1_5p" for n in range(len(midata.var_names))]
    midata.var['reference'] = "miRNA"

    fmi = f"{tmp}/midata.h5ad"
    midata.uns = {}
    midata.write(fmi)

    print(midata.obs.index)
    adata = tk.stitch_mRNA_and_miRNA(fma, fmi, mrna_umi_cutoff=100, mirna_umi_cutoff=10, simplify_mirnames=True, protein_coding_genes_path=fcod),

