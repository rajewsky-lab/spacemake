library(data.table)

# pass the arguments from command line
args <- commandArgs(trailingOnly = TRUE)
folder = args[1]
threshold = as.integer(args[2])

dge = fread(paste0('zcat < ', folder, '/dge.txt.gz'))
dge_reads = fread(paste0('zcat < ', folder, '/dgeReads.txt.gz'))

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
fwrite(df, paste0(folder, '/output_qc_sheet/downstream_statistics.csv'),
       row.names=TRUE)
