__version__ = '0.1.0'
__author__ = ['Nikos Karaiskos', 'Tamas Ryszard Sztanka-Toth']
__licence__ = 'GPL'
__email__ = ['nikolaos.karaiskos@mdc-berlin.de', 'tamasryszard.sztanka-toth@mdc-berlin.de']


import pandas as pd
import yaml
import scanpy as sc
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from fpdf import FPDF
import pickle

figures_root= snakemake.params['fig_root']
results_file = snakemake.input['res_file']

figures_suffix = '_' + snakemake.wildcards['united_sample'] + '_' + snakemake.wildcards['puck'] + '.png'

sc.settings.figdir = figures_root
sc.settings.file_format_figs = 'png'

################
# LOAD RESULTS #
################
adata = sc.read(results_file)
nrow, ncol = adata.shape

if 'rank_genes_groups' in adata.uns_keys():
    top_markers = adata.uns['rank_genes_groups']['names'][:2]

##############
# SAVE PLOTS #
##############
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
                     jitter=0.4, multi_panel=True, save='_filtered' + figures_suffix,
                     log=True)

# only print the first 5 pairs of components
if 'neighbors' in adata.uns_keys():
    pcs = [str(x) for x in range(1, adata.uns['neighbors']['params']['n_pcs'] + 1)]
else:
    pcs = []

pc_components = [','.join([a,b]) for a, b in zip(pcs[::2], pcs[1::2])] 

if pcs == [] or pc_components == []:
    plt.plot()
    plt.savefig(figures_root + '/pca_first_components' + figures_suffix)
    plt.close()
else:
    sc.pl.pca(adata, color='leiden', components = pc_components[:6], save='_first_components' + figures_suffix, ncols=3, title=['']*6)

if 'umap' in adata.uns_keys() and 'X_umap' in adata.obsm_keys():
    # save umap
    sc.pl.umap(adata, color = 'leiden', save='_clusters' + figures_suffix, size = 10, title='')

    # save small umap per marker gene per cluster
    sc.pl.umap(adata, color = top_markers[0], save='_top1_markers' + figures_suffix, color_map=mpl.cm.Reds)
    sc.pl.umap(adata, color = top_markers[1], save='_top2_markers' + figures_suffix, color_map=mpl.cm.Reds)
else:
    plt.plot()
    # save two empty figures as top marker plots
    plt.savefig(figures_root + '/umap_clusters' + figures_suffix)
    plt.savefig(figures_root + '/umap_top1_markers' + figures_suffix)
    plt.savefig(figures_root + '/umap_top2_markers' + figures_suffix)

##############
# SAVE LISTS #
##############
# save umap coordinates, top cluster
if ('X_umap', 0) in adata.obsm_keys():
    obsm_keys = [('X_umap', 0), ('X_umap', 1)]
else:
    obsm_keys = []

sc.get.obs_df(adata, keys=adata.obs_keys(), obsm_keys=obsm_keys)\
    .rename(columns={'X_umap-0':'umap_0', 'X_umap-1':'umap_1', 'leiden':'cluster', 'n_counts':'n_umis'})\
    .to_csv(snakemake.output['results_metadata'], index_label='cell_bc')

#######################
# GENERATE PDF REPORT #
#######################
with open(snakemake.input['parameters_file']) as f:
    parameters = yaml.load(f, Loader=yaml.FullLoader)

read_statistics = {}
with open(snakemake.input['star_log'], 'r') as fi:
    idx = 0
    for line in fi.readlines():
        entry = line.strip('\n').split('\t')
        if idx == 5:
            read_statistics['input_reads'] = float(entry[1])
        if idx == 8:
            read_statistics['uniquely_mapped'] = float(entry[1])
        idx += 1


pdf = FPDF()
pdf.add_page()
pdf.set_xy(0, 0)

# title
pdf.set_font('Arial', 'B', 16)
pdf.cell(90, 20, " ", 0, 1, 'C')
pdf.cell(10)
pdf.cell(100, 8, "Automated analysis of a spatial illumina sequencing run", 0, 2, 'L')
pdf.set_font('Arial', '', 11)
pdf.cell(100, 8, "tamasryszard.sztanka-toth@mdc-berlin.de", 0, 2, 'L')

# sample info
pdf.set_font('Arial', 'B', 14)
pdf.cell(100, 10, "Sample info", 0, 2, 'L')
pdf.set_font('Arial', '', 11)
pdf.cell(100, 5, "project_id: %s" % (parameters['project_id']), 0, 2, 'L')
pdf.cell(100, 5, "sample_id: %s" % (parameters['sample_id']), 0, 2, 'L')
pdf.cell(100, 5, "puck_id: %s" % (parameters['puck_id']), 0, 2, 'l')
# sometimes we have more experiment descriptions. this is when we have merged samples.
# in this case we need to print each experiment desc, one below the other
# first just print one cell titled experiment
pdf.cell(25, 5, "experiment: ", 0, 0, 'L')
# then right of it print one cell per experiment description. for merged they are split by ','
for experiment in parameters['experiment'].split(','):
    pdf.cell(75, 5, experiment, 0, 2, 'L')
pdf.cell(75, 0, ' ', 0, 1, 'L')
pdf.cell(10)
pdf.cell(100, 5, "sequencing_date: %s" % (', '.join(parameters['sequencing_date'].split(','))), 0, 2, 'L')
pdf.cell(100, 5, "investigator(s): %s" % (', '.join(parameters['investigator'].split(','))), 0, 2, 'L')
# print mapping stats
input_reads = int(read_statistics['input_reads'])
mapped_reads = int(read_statistics['uniquely_mapped'])
pdf.cell(100, 5, "input reads: %s" % (f'{input_reads:,}'), 0, 2, 'l')
pdf.cell(100, 5, "uniquely mapped reads: %s (%s)" % (f'{mapped_reads:,}', round(mapped_reads/input_reads, 3)), 0, 2, 'l')

# Downstream analysis
pdf.add_page()
pdf.cell(10)
pdf.set_font('Arial', 'B', 14)
pdf.cell(100, 10, 'Downstream results', 0, 2, 'L')
pdf.set_font('Arial', '', 11)

n_beads = adata.obs.shape[0]

try:
    n_clusters = len(adata.obs.leiden.dtypes.categories)
except AttributeError:
    n_clusters = 0

pdf.cell(100, 5, "%s Beads passing filtering, in %s clusters." % (n_beads, n_clusters), 0, 2, 'L')
pdf.cell(100, 5,'Filtered with percent_mito < 0.05, n_genes < 2500 per bead, and min ' + snakemake.wildcards['umi_cutoff'] + ' UMIs per bead', 0, 2, 'L')

pdf.set_font('Arial', 'B', 11)
pdf.cell(100, 10, 'Violin plots of genes, umis, percent_mito after filtering', 0, 2, 'L')
pdf.image(figures_root + '/violin_filtered' + figures_suffix, w = 150, h=50)
pdf.cell(100, 10, 'PCA plots for PC1-PC8 based on highly variable genes', 0, 2, 'L')
pdf.image(figures_root + '/pca_first_components' + figures_suffix, w = 160, h=100)

# UMAP
pdf.add_page()
pdf.cell(10)
pdf.set_font('Arial', 'B', 14)
pdf.cell(100, 10, "Clustering analysis", 0, 2, 'L')

# UMAP plot
pdf.set_font('Arial', 'B', 11)
pdf.cell(100, 10, 'UMAP', 0, 2, 'L')
pdf.image(figures_root + '/umap_clusters' + figures_suffix, w = 100, h = 80)

# clustering parameters 
pdf.set_font('Arial', 'B', 11)
pdf.cell(100, 10, 'Clustering parameters used', 0, 2, 'L')
pdf.set_font('Arial', '', 11)
try:
    pdf.cell(100, 5, "# of highly variable genes: %s" % (adata.var.highly_variable.value_counts().loc[True]), 0, 2, 'L')
except AttributeError:
    pass

if 'leiden' in adata.uns_keys():
    pdf.cell(100, 5, "clustering resolution: %s" % (adata.uns['leiden']['params']['resolution']), 0, 2, 'L')
if 'neighbors' in adata.uns_keys():
    pdf.cell(100, 5, "number of PCs used for clustering: %s" % (adata.uns['neighbors']['params']['n_pcs']), 0, 2, 'L')

# top marker genes table
pdf.set_font('Arial', 'B', 11)
pdf.cell(100, 10, 'Top marker genes per cluster', 0, 2, 'L')
pdf.set_font('Arial', '', 11)
pdf.cell(100, 10, 'Top 15 marker genes for each cluster, when compared to other clusters.', 0, 2, 'L')
pdf.cell(100, 10, 'The higher the gene is in the table, the more differential it is', 0, 1, 'L')

if 'rank_genes_groups' in adata.uns_keys():
    marker_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(15)
    table_split = np.int(np.log10(marker_genes.columns.shape[0]))

    for split in range(table_split + 1):
        pdf.set_font('Arial', 'B', 10)
        pdf.cell(10)
        for c in marker_genes.columns[split*10:(split+1)*10]:
            pdf.cell(18, 4, c, 1, 0, 'C')
        
        pdf.cell(0, 4, " ", 0, 1, 'L') 
        pdf.set_font('Arial', '', 10)
        for index, row in marker_genes.iloc[:, split*10:(split+1)*10].iterrows():
            pdf.cell(10)
            for item in row:
                pdf.cell(18, 4, item, 1, 0, 'C')
            pdf.cell(0, 4, " ", 0, 1, 'L') 
        pdf.cell(0, 4, " ", 0, 1, 'L')

# Marker gene plots
pdf.add_page()
pdf.cell(10)
pdf.set_font('Arial', 'B', 14)
pdf.cell(90, 20, " ", 0, 2, 'C')
pdf.cell(100, 10, "Top marker gene plots - top marker", 0, 2, 'L')
pdf.image(figures_root + '/umap_top1_markers' + figures_suffix, w = 180, h = 100)
pdf.cell(100, 20, " ", 0, 2, 'c')

pdf.add_page()
pdf.cell(10)
pdf.set_font('Arial', 'B', 14)
pdf.cell(90, 20, " ", 0, 2, 'C')
pdf.cell(100, 10, "Top marker gene plots - second top marker", 0, 2, 'L')
pdf.image(figures_root + '/umap_top2_markers' + figures_suffix, w = 180, h = 100)

pdf.output(snakemake.output['report'], 'F')
