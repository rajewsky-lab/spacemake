import pysam
import pandas as pd
from spacemake.util import FASTQ_src
from byo.track import load_track
from collections import defaultdict
import re

gene_starts = defaultdict(lambda: int(2 ** 32))
gene_ends = defaultdict(int)

flank = 200

# # extract gene loci from GTF
# with open("test_annotation.gtf", "wt") as f:
for line in open("gencode.v38.chr22.gtf", "rt"):
    if line.startswith("#"):
        continue

    if not "protein_coding" in line:
        continue

    # print(line)
    parts = line.rstrip().split("\t")
    chrom, source, rec, start, end = parts[:5]

    if chrom != "chr22":
        continue

    if rec != "gene":
        continue

    start = int(start)
    end = int(end)

    m = re.search(r'gene_name "(\S+)"', parts[-1])
    if not m:
        continue

    gene = m.groups()[0]

    gene_starts[gene] = min(gene_starts[gene], start - flank)
    gene_ends[gene] = max(gene_ends[gene], end + flank)

    # kb = (gene_ends[gene] - gene_starts[gene]) / 1000
    # print(f" {gene} {gene_starts[gene]}-{gene_ends[gene]} -> L={kb:.1f}kb")


for gene in sorted(gene_starts.keys()):
    print(
        f"{gene}\t{gene_starts[gene]}\t{gene_ends[gene]}\t{gene_ends[gene] - gene_starts[gene]}"
    )
