import pysam
import pandas as pd
from spacemake.util import FASTQ_src, read_fq
from byo.track import load_track

genome = load_track("/data/rajewsky/genomes/hg38/hg38.fa")

# load output of gene_loci_to_gtf.py: a table with gene-name, start and end coordinates
df_genes = pd.read_csv(
    "chr22_gene_bounds.csv", sep="\t", names=["gene", "start", "end", "L"]
).set_index("gene")

# keep track of the few dozen reads that Nikos has selected
selected_reads = {}
for fa_id, seq, qual in read_fq("reads_chr22_R2.fastq.gz"):
    if "IGLC3" in fa_id:
        print(f"selecting IGLC3 read {fa_id} -> {fa_id.split('_')[0]}")
    selected_reads[fa_id.split("_")[0]] = fa_id

# then go through the SAM file (no header) to get the mapping position of these reads (not ideal)
df = pd.read_csv(
    "/data/rajewsky/home/nkarais/murphy/fc_sts/collect_reads_chr22/final.polyA_adapter_trimmed_chr22.sam",
    sep="\t",
    header=None,
)

starts = []
for row in df.itertuples():
    # print(row)
    qname = row[1]
    if "A00643:496:HFJ5MDRX2:1:2101:12888:1172" in qname:
        print(f"YAY! detected IGLC3 read {qname}")

    if qname in selected_reads:
        print(f"selecting read {qname}")
        starts.append(row[4])

# find genes that overlap the selected reads mapping position
# this intersection code is very crude, but effective
intervals = set()
starts = set(starts)
for row in df_genes.itertuples():
    next_starts = set(starts)
    for x in starts:
        if row.start < x and row.end > x:
            intervals.add((row.start, row.end))
            print(f"selecting gene entry '{row}'")
            next_starts.discard(x)

    starts = next_starts

print(
    f"we have the following start coordinates left. selecting buffer regions around these"
)
print(starts)


def do_merge(s, e, intervals):
    keep = []
    for j, (s2, e2) in enumerate(intervals):
        new = (s2, e2)
        if s2 <= e and e2 >= e:
            print(f"overlap on the right, s={s} e={e} s2={s2} e2={e2}")
            new = (min(s2, s), max(e, e2))
        elif e2 >= s and s2 <= s:
            print(f"overlap on the left, s={s} e={e} s2={s2} e2={e2}")
            new = (min(s2, s), max(e, e2))
        elif s2 >= s and e2 <= e:
            print(f"contained in other interval. discard, s={s} e={e} s2={s2} e2={e2}")
            continue

        keep.append(new)

    return keep


# merge intervals that have some overlap
intervals = sorted(list(intervals), key=lambda x: x[1] - x[0], reverse=True)
print(intervals)

while True:
    changed = False
    for i, (s, e) in enumerate(intervals):
        others = intervals[i + 1 :]
        keep = do_merge(s, e, others)
        if keep != others:
            print("we had a change!")
            print("before")
            for s, e in intervals:
                print(f"{s} - {e}")

            intervals = intervals[: i + 1] + keep
            intervals = sorted(list(intervals), key=lambda x: x[1] - x[0], reverse=True)
            print("after")
            for s, e in intervals:
                print(f"{s} - {e}")

            changed = True
            break  # the for loop

    if not changed:
        break

intervals = sorted(list(set(intervals)), key=lambda x: x[1] - x[0], reverse=True)
print("remaining intervals")
for s, e in intervals:
    print(f"{s} - {e}")


# Okay, now we know the gene loci which are needed to map the test reads!
intervals = list(sorted(intervals))
print(f"relevant intervals found: {len(intervals)}")

# extract the genomic sequence for the loci we need and save as a mini-"genome"
with open("test_genome.fa", "wt") as f:
    for start, end in intervals:
        seq = genome.get_oriented("chr22", start, end, "+")
        f.write(f">test_chr22.{start}-{end}\n{seq}\n")

# cut down the GTF annotation to only those parts that pertain to the genic regions we care about
with open("test_annotation.gtf", "wt") as f:
    for line in open("gencode.v38.chr22.gtf", "rt"):
        if line.startswith("#"):
            continue

        parts = line.rstrip().split("\t")
        chrom, source, rec, start, end = parts[:5]
        if chrom != "chr22":
            continue

        start = int(start)
        end = int(end)

        for s, e in intervals:
            if (
                (start <= s and end > s)  # overlap the start of interval
                or (start > s and end < e)  # internal to interval
                or (start < e and end > e)  # overlap the end of interval
                or (start < s and end > e)  # overlap the entire interval
            ):
                # the name of the pseudo-chromosome this is on (see excision of genomic sequence above)
                chrom = f"test_{chrom}.{s}-{e}"

                start = max(
                    0, start - s
                )  # translate start and end coordinates from whole chr22 to the gene region
                end = min(e - s, end - s)

                parts[0:5] = (chrom, "test_data", rec, str(start), str(end))
                f.write("\t".join(parts) + "\n")

# Done! We now have:
# * test_genome.fa.gz with the gene sequences
# * test_annotation.gtf.gz with the gene models (exon/intron, CDS/UTR etc.)
# * reads_chr22_R1.fastq.gz with test read barcodes mapping to a few CBs and UMIs
# * reads_chr22_R2.fastq.gz with test read cDNAs mapping to the genic regions we care about
