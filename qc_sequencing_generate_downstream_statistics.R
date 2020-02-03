###
# file is called from snakemake. snakemake is a global variable holding all the required input
# and output file names
###


library(data.table)
library(yaml)

parameters = read_yaml(snakemake@input[['parameters']])

threshold = parameters$threshold

print('Loading the DGE...')
dge = fread(paste0('zcat < ', snakemake@input[['dge']]))
dge_reads = fread(paste0('zcat < ', snakemake@input[['dgeReads']]))

# find beads which have the minimum number of UMIs
sum_umis = colSums(dge[, 3:dim(dge)[2]])
beads = which(sum_umis >= threshold)

# compute reads per bead
reads_per_bead = colSums(dge_reads[, names(beads), with=FALSE])

# compute genes per bead
genes_per_bead = colSums(dge[, names(beads), with=FALSE] > 0)

# compute UMIs per bead
umis_per_bead = colSums(dge[, names(beads), with=FALSE])

# write to disk to read it and plot in python
# (might want to plot directly in R at some point)
df = data.frame('reads'=reads_per_bead,
                'genes'=genes_per_bead,
		'umis'=umis_per_bead)
fwrite(df, snakemake@output[[1]], row.names=TRUE)
