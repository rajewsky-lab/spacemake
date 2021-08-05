def ensure_path(path):
    import os

    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path


def read_fq(fname):
    import gzip
    from more_itertools import grouper

    if fname.endswith(".gz"):
        src = gzip.open(fname, mode="rt")
    elif type(fname) is str:
        src = open(fname)
    else:
        src = fname  # assume its a stream or file-like object already

    for name, seq, _, qual in grouper(src, 4):
        yield name.rstrip()[1:], seq.rstrip(), qual.rstrip()

def dge_to_sparse(dge_path, out_path):
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
        M = np.zeros((1000, ncol))

        for line in dge:
            vals = line.strip().split('\t')
            # first element contains the gene name
            gene_names.append(vals[0])
            vals = vals[1:]
            vals = np.array(vals, dtype=np.int64)

            M[ix] = vals

            if ix % 1000 == 999:
                print(mix)
                mix = mix + 1
                matrices.append(csr_matrix(M))
                ix = 0
                M = np.zeros((1000, ncol))
            else:
                ix = ix + 1

        # get the leftovers
        M = M[:ix]
        matrices.append(csr_matrix(M))

        # sparse expression matrix
        X = vstack(matrices, format='csr')
        print(len(gene_names))
        
        adata = anndata.AnnData(X, var = pd.DataFrame(index=barcodes), obs = pd.DataFrame(index=gene_names))

        adata.write(out_path)
