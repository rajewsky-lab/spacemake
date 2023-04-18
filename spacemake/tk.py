"""
Some useful functions for single cell & spatial transcriptomic analysis with scanpy
"""

import scanpy as sc
import numpy as np
import pandas as pd

def bootstrap(adata, cell_mask=None, genes=[], n_bootstrap=100):
    if cell_mask is not None:
        barcodes = adata.obs_names[cell_mask]
    else:
        barcodes = adata.obs_names
    if not len(genes):
        genes = adata.var_names
    
    counts = []
    raw = adata.raw[barcodes, genes].X #.copy()
    for i in range(n_bootstrap):
        bootstrap = np.random.randint(0, high=len(barcodes), size=len(barcodes))
        counts.append(np.array(raw[bootstrap, :].sum(axis=0), dtype=int)[0])
    
    counts = np.array(counts)
    return counts

def pseudo_bulk(adata, cell_mask=None, genes=[], n_bootstrap=100, q=[2.5,50,97.5]):
    """
    samples cells from the adata object n_bootstrap times with replacement. 
    Then determines mean, sum and low, mid, and upper quantiles (default=[2.5, 50, 97.5])
    across the bootstraps. With the default values hi - low will give 95% confidence interval
    and the main value is the median.
    
    The results are returned as a DataFrame.
    """
    if not len(genes):
        genes = adata.var_names.to_list()

    counts = bootstrap(adata, cell_mask, genes, n_bootstrap)
    quantiles = np.percentile(counts, axis=0, q=q)

    # print(counts.shape, raw.shape, quantiles.shape)
    d = {
        'name': genes,
        # 'sum': np.array(raw.sum(axis=0), dtype=int)[0],
        'mean': np.mean(counts, axis=0),
        'lo': quantiles[0],
        'med': quantiles[1],
        'hi': quantiles[2]
    }
    # for k, v in d.items():
    #     print(k, len(v))

    return pd.DataFrame(d).set_index('name').sort_values('mean', ascending=False)


def pseudo_bulk_clusters(adata, genes=[], groupby='leiden', **kw):
    """
    runs pseudo_bulk on cells grouped by a clustering method (default='leiden'),
    but can also be used for user-defined labels. Just change groupby=... to
    sth in .obs[]
    Results are again returned as a DataFrame
    """

    data = {}
    print(f"sampling all")
    df = pseudo_bulk(adata, genes=genes, **kw)
    data['all'] = df['med']
    data['all_lo'] = df['lo']
    data['all_hi'] = df['hi']

    for clst in set(adata.obs[groupby]):
        print(f"sampling {clst}")
        df = pseudo_bulk(adata, cell_mask=(adata.obs[groupby] == clst), genes=genes, **kw)
        data[f"cluster_{clst}"] = df['med']
        data[f"cluster_{clst}_lo"] = df['lo']
        data[f"cluster_{clst}_hi"] = df['hi']

    df = pd.DataFrame(data)
    return df


def sample_from_clusters(adata, clusters=['15'], groupby='leiden', n_cells=1):
    barcodes = []
    for clst in clusters:
        m = adata.obs[groupby] == clst
        clst_cells = adata.obs_names[m]
        cbs = np.random.choice(clst_cells, replace=False, size=n_cells)
        #print(cbs)
        barcodes.extend(cbs)
    
    data = {'name': mirnames}
    for cb in barcodes:
        counts = np.array(raw[cb, mirnames].X.todense(), dtype=int)[0]
        data[cb] = counts
    
    return pd.DataFrame(data).set_index('name').sort_values(cb, ascending=False)
        

def merge_adata_objects_on_cells(adata1, adata2, genes1=None, genes2=None):
    obs1 = set(adata1.obs_names)
    obs2 = set(adata2.obs_names)
    obs_both = obs1 & obs2
    obs_both -= set(['NA'])

    n1 = len(obs1)
    n2 = len(obs2)
    n_both = len(obs_both)

    f = n_both/min(n1, n2)
    print(f"Fraction of shared cell barcodes between adata1 and adata2 is {f:.3f}")

    # create aligned copies of the data
    obs_names = sorted(obs_both)
    _1 = adata1[obs_names,:].copy()
    _2 = adata2[obs_names,:].copy()

    print("about to merge")
    print(_1.X.shape)
    print(_2.X.shape)

    # make a new AnnData object, concatenated along the var (genes) axis
    import scipy.sparse

    adata = sc.AnnData(X=scipy.sparse.hstack([
        scipy.sparse.csr_matrix(_1.X), 
        scipy.sparse.csr_matrix(_2.X)
        ]), dtype=np.float32)

    adata.obs_names = obs_names
    adata.var_names = _1.var_names.tolist() + _2.var_names.tolist()

    adata.uns['merge_barcode_overlap'] = f
    return adata


def aggregate_gene_counts(adata, func=None):

    def simplify(mirname):
        """
        Remove the precursor name and modification id from the miRNA name.
        Example: 
        
            Mir-124-P1-v1_3p -> Mir-124_3p

        """
        import re
        M = re.search(r"^(\w+)-(\w+)\-(\d+)(.*?)_(5|3)p", mirname)
        if not M:
            return mirname
        else:
            species, a, b, c, d = M.groups()
            return f"{species}-{a}-{b}_{d}p"

    if func is None:
        func = simplify

    df = pd.DataFrame({'full': adata.var_names, 'short': [simplify(m) for m in adata.var_names]})
    new_counts = {}
    for name, group in df.groupby('short'):
        #print(name)
        new_counts[name] = adata[:, group['full']].X.sum(axis=1)

    simpler_names = df['short'].drop_duplicates().sort_values()
    import scipy.sparse
    X = scipy.sparse.csr_matrix(np.array([new_counts[s] for s in simpler_names]).T[0])

    agg_data = sc.AnnData(X=X, dtype=int)
    agg_data.obs_names = adata.obs_names
    agg_data.var_names = simpler_names

    agg_data.uns = adata.uns.copy()
    agg_data.uns['name_map'] = df
    agg_data.obs = adata.obs.copy()
    agg_data.obsm = adata.obsm.copy()
    agg_data.obsp = adata.obsp.copy()

    adata.uns['name_map'] = df

    print(f"aggregated counts. We started with {adata.X.shape} and now we have {agg_data.X.shape}")
    return agg_data


def stitch_mRNA_and_miRNA(mdata_path, midata_path, mrna_umi_cutoff=1000, mirna_umi_cutoff=25, simplify_mirnames=True, protein_coding_genes_path='', normalize=True, mt_gene_pattern="^mt-"):
    # load
    import re
    mdata = sc.read_h5ad(mdata_path)
    if 'NA' in mdata.obs_names:
        ambient_mrna = mdata['NA'].to_df().T['NA']
    else:
        ambient_mrna = 0
    mdata.var['ambient'] = ambient_mrna
    # print(f"initial mdata {mdata}")
    
    midata = sc.read_h5ad(midata_path)
    # print(f"initial midata {midata}")
    if 'NA' in midata.obs_names:
        ambient_mirna = midata['NA'].to_df().T['NA']
    else:
        ambient_mirna = 0
    midata.var['ambient'] = ambient_mirna

    # pre-filter miRNA -> apply UMI cutoff and select only miRNA genes
    mirna_genes = midata.var['reference'] == 'miRNA'
    print(f"found {mirna_genes.sum()} miRNA genes")
    midata.obs['miRNA_counts'] = midata[:,mirna_genes].X.sum(axis=1)
    m = (midata.obs['miRNA_counts'] >= mirna_umi_cutoff) #& (adata.obs['n_counts'] < 10000000)
    midata = midata[m, mirna_genes].copy()
    print(f"pre-filtered midata {midata.X.shape}")
    print(f"pre-filtered midata var {midata.var['reference']}")

    # pre-filter mRNA -> apply UMI cutoff and select only genes from the genome index
    # print(f"mdata {mdata}")
    genome_genes = mdata.var['reference'] == 'genome'
    mdata.obs['genome_counts'] = mdata[:, genome_genes].X.sum(axis=1)
    # print(f"genome_counts quantiles: {np.percentile(mdata.obs['genome_counts'], [1, 5, 25, 50, 75, 95, 99])}")
    m = mdata.obs['genome_counts'] >= mrna_umi_cutoff
    mdata = mdata[m, genome_genes].copy()
    print(f"pre-filtered mdata {mdata.X.shape}")
    print(f"pre-filtered mdata var {mdata.var['reference']}")

    # simplify miRNA names if requested
    print(f"copying midata")
    midata_original = midata.copy()
    if simplify_mirnames:
        print(f"simplifying")
        midata = aggregate_gene_counts(midata)
        print(f"simplified midata {midata.X.shape}")
    
    # normalize separately, if requested
    if normalize:
        print(f"copy before norm")
        mdata_norm = mdata.copy()
        midata_norm = midata.copy()
        # print(f">> mdata_norm before normalize_total {mdata_norm}")
        print(f"normalizing")
        sc.pp.normalize_total(mdata_norm, 1e4)
        sc.pp.normalize_total(midata_norm, 1e4)

        sc.pp.log1p(mdata_norm)
        sc.pp.log1p(midata_norm)

        # print(f">> mdata_norm before scaling {mdata_norm}")
        sc.pp.scale(mdata_norm, max_value=10)
        sc.pp.scale(midata_norm, max_value=10)
    else:
        mdata_norm = mdata
        midata_norm = midata

    print("merging")
    # merge both AnnData objects over common cells
    adata = merge_adata_objects_on_cells(mdata_norm, midata_norm)
    print(f"merged, normalized data {adata}")
    raw = merge_adata_objects_on_cells(mdata, midata)
    print(f"raw data: {raw}")
    adata.raw = raw

    adata.var['reference'] = pd.Series(
        mdata.var['reference'].tolist() + (['miRNA'] * midata.X.shape[1])
    ).astype('category')

    adata.obs['n_counts'] = raw.X.sum(axis=1)

    print(f"we have the following references {adata.var['reference'].drop_duplicates()}.")
    miRNA_genes = adata.var_names[adata.var['reference'] == 'miRNA'].to_list()
    genome_genes = adata.var_names[adata.var['reference'] == 'genome'].to_list()
    mt_genes = [g for g in adata.var_names if re.search(mt_gene_pattern, g.lower())]
    adata.uns['mt_genes'] = sorted(mt_genes)
    adata.uns['miRNA_genes'] = sorted(miRNA_genes)
    adata.uns['genome_genes'] = sorted(genome_genes)
    print(f"detected {len(mt_genes)} mitochondrial gene names")
    print(f"detected {len(miRNA_genes)} miRNA gene names")
    print(f"detected {len(genome_genes)} genome gene names")

    protein_coding = set([g.strip() for g in open(protein_coding_genes_path)])
    print(f"detected {len(protein_coding)} protein coding gene names")
    adata.uns['protein_coding_genes'] = sorted(protein_coding)

    adata.obs['n_miRNA_counts'] = raw[:, miRNA_genes].X.sum(axis=1)
    adata.obs['n_genome_counts'] = raw[:, genome_genes].X.sum(axis=1)
    adata.obs['n_mt_counts'] = raw[:, mt_genes].X.sum(axis=1)
    adata.var['protein_coding'] = adata.var_names.isin(protein_coding)
    adata.obs['n_coding_counts'] = raw[:, adata.var['protein_coding']].X.sum(axis=1)

    adata.obs['pct_coding'] = (100.0 * adata.obs['n_coding_counts']) / adata.obs['n_counts']
    adata.obs['pct_genome'] = (100.0 * adata.obs['n_genome_counts']) / adata.obs['n_counts']
    adata.obs['pct_miRNA'] = (100.0 * adata.obs['n_miRNA_counts']) / adata.obs['n_counts']
    adata.obs['pct_mt'] = (100.0 * adata.obs['n_mt_counts']) / adata.obs['n_counts']

    raw.obs = adata.obs.copy()
    # sanity check
    assert adata[:, miRNA_genes].var['protein_coding'].sum() == 0

    adata.uns['ambient_miRNA'] = midata_original.var['ambient']
    adata.uns['ambient_mRNA'] = mdata.var['ambient']
    return adata, raw, mdata[adata.obs_names], midata_original[adata.obs_names]


def assign_species(adata, species_patterns, pct_thresh=80):
    """
    adds up the counts for species-specific genes (e.g. by prefix hg_...) and populates the following adata.obs columns:
    
        - adata.obs[f'n_{species}_counts'] sum of counts of genes that match the pattern
        - adata.obs[f'n_all_species_counts'] above, but summed over all species that are in the list
        - adata.obs[f'pct_{species}'] percent of species-assignable counts that come from that species
        - adata.obs['species'] a categorical that is set to the name of a species if pct_{species} is >= pct_thresh
        - adata.obs[f'pct_species'] percent of species-assignable counts that come from the assigned (above threshold) species

    If the threshold is not reached for any species, the cell is assigned 'nan' as value of the 'species' field and 0 for 'pct_species'
    """
    import re
    adata.obs['species'] = 'nan'
    adata.obs['n_all_species_counts'] = 0
    adata.obs['pct_species'] = np.nan

    for species, var_name_pattern in species_patterns:
        species_genes = sorted([gene for gene in adata.var_names if re.search(var_name_pattern, gene)])
        adata.uns[f"{species}_genes"] = species_genes
        adata.obs[f"n_{species}_counts"] = np.array(adata.raw[:, species_genes].X.sum(axis=1))[:, 0]
        adata.obs['n_all_species_counts'] += adata.obs[f"n_{species}_counts"]

    masks = []
    for species, var_name_pattern in species_patterns:
        adata.obs[f"pct_{species}"] = 100.0 * adata.obs[f"n_{species}_counts"] / adata.obs["n_all_species_counts"]
    
        # assign the species to cells where the fraction of counts coming from the species' genes exceeds
        # the threshold
        m = adata.obs[f"pct_{species}"] >= pct_thresh
        adata.obs.loc[m, 'pct_species'] = adata.obs.loc[m, f'pct_{species}']
        adata.obs.loc[m, "species"] = species
        masks.append(m)

    return masks


def pre_filter(adata, umi_cutoff=0, umi_key="n_counts", normalize=True, mt_gene_pattern="^mt-", protein_coding_genes_path='/data/rajewsky/projects/sc_smRNA_marvin/sm/species_data/mouse/genome/protein_coding_genes.txt'):
    import re
    genome_genes = adata.var['reference'] == 'genome'
    adata.uns['genome_genes'] = genome_genes
    adata.obs['genome_counts'] = adata[:,genome_genes].X.sum(axis=1)
#    m = (mdata.obs['genome_counts'] > mrna_umi_cutoff) #& (mdata.obs['genome_counts'] < 10000)
    mt_genes = [g for g in adata.var_names if re.search(mt_gene_pattern, g.lower())]
    adata.uns['mt_genes'] = sorted(mt_genes)
    print(f"detected {len(mt_genes)} mitochondrial gene names")

    protein_coding = set([g.strip() for g in open(protein_coding_genes_path)])
    print(f"detected {len(protein_coding)} protein coding gene names")
    adata.uns['protein_coding_genes'] = sorted(protein_coding)

    adata.obs['n_genome_counts'] = adata[:, genome_genes].X.sum(axis=1)
    adata.obs['n_mt_counts'] = adata[:, mt_genes].X.sum(axis=1)
    adata.var['protein_coding'] = adata.var_names.isin(protein_coding)
    adata.obs['n_coding_counts'] = adata[:, adata.var['protein_coding']].X.sum(axis=1)

    adata.obs['pct_coding'] = (100.0 * adata.obs['n_coding_counts']) / adata.obs['n_counts']
    adata.obs['pct_mt'] = (100.0 * adata.obs['n_mt_counts']) / adata.obs['n_counts']
    
    # remove the low-count barcodes lumped together already at the quant.py stage
    m = adata.obs_names != 'NA'
    adata = adata[m, :]
    
    if umi_cutoff:
        m = (adata.obs[umi_key] > umi_cutoff)
        adata = adata[m, :]

    raw = adata.copy()
    if normalize:
        sc.pp.normalize_total(adata, 1e4)
        # print(adata.X)
        sc.pp.log1p(adata)

    adata.raw = raw
    return adata, raw



def compare_to_bulk(adata, bulk, genes=[]):

    import scipy
    bulk = pd.DataFrame({
        'name' : bulk.var_names, 
        'bulk': np.array(bulk.X.sum(axis=0))[0]
    }).set_index('name')

    if not len(genes):
        genes = adata.var_names
    gdata = adata.raw[:, genes]

    pseudo = pd.DataFrame({
        'name' : gdata.var_names,
        'pseudo' : np.array(gdata.X.sum(axis=0))[0]
    }).set_index('name')

    df = pd.concat([bulk, pseudo], axis=1).fillna(0)

    print(scipy.stats.pearsonr(np.log10(df['bulk'] + 1), np.log10(df['pseudo'] + 1)))
    
    df["log2FC"] = np.log2((df["pseudo"] + 1)/ (df["bulk"] + 1) )
    df = df.sort_values("log2FC", ascending=False)
    df.sort_values('log2FC', ascending=False).head(30)  
    
    return df


def add_common_metrics(adata, mt_gene_pattern="^mt-", protein_coding_genes_path='/data/rajewsky/projects/sc_smRNA_marvin/sm/species_data/mouse/genome/protein_coding_genes.txt'):
    import scanpy as sc
    import re
    # lose all zero rows/columns
    sc.pp.filter_cells(adata, min_counts=1)
    sc.pp.filter_genes(adata, min_counts=1)

    # re-generate the metrics/marginals of the adata matrix
    adata.obs["n_counts"] = adata.X.sum(axis=1)
    adata.obs["n_genes"] = (adata.X > 0).sum(axis=1)
    adata.var["n_cells"] = (adata.X > 0).sum(axis=0).T

    # and all extra channels/layers
    for l in sorted(adata.layers.keys()):
        adata.obs[f"n_{l}"] = adata.layers[l].sum(axis=1)[:, 0]
    
    # this assumes that we have a layer["reads"], which should be the case 
    # if the h5ad comes from quant.py output
    adata.obs["reads_per_counts"] = adata.obs["n_reads"] / adata.obs["n_counts"]

    genome_genes = adata.var['reference'] == 'genome'
    adata.uns['genome_genes'] = genome_genes
    adata.obs['genome_counts'] = adata[:,genome_genes].X.sum(axis=1)
#    m = (mdata.obs['genome_counts'] > mrna_umi_cutoff) #& (mdata.obs['genome_counts'] < 10000)
    mt_genes = [g for g in adata.var_names if re.search(mt_gene_pattern, g.lower())]
    adata.uns['mt_genes'] = sorted(mt_genes)
    print(f"detected {len(mt_genes)} mitochondrial gene names")

    protein_coding = set([g.strip() for g in open(protein_coding_genes_path)])
    print(f"detected {len(protein_coding)} protein coding gene names")
    adata.uns['protein_coding_genes'] = sorted(protein_coding)

    adata.obs['n_genome_counts'] = adata[:, genome_genes].X.sum(axis=1)
    adata.obs['n_mt_counts'] = adata[:, mt_genes].X.sum(axis=1)
    adata.var['protein_coding'] = adata.var_names.isin(protein_coding)
    adata.obs['n_coding_counts'] = adata[:, adata.var['protein_coding']].X.sum(axis=1)

    adata.obs['pct_coding'] = (100.0 * adata.obs['n_coding_counts']) / adata.obs['n_counts']
    adata.obs['pct_mt'] = (100.0 * adata.obs['n_mt_counts']) / adata.obs['n_counts']

    return adata

