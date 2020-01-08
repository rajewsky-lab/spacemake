#!/bin/bash
#
# Programs used:
# - picard-tools
# - the Drop-Seq toolkit
#
# The pipeline assumes that:
# - two (paired end) fastq files are present in the working directory
# - an indexed reference genome exists somewhere (with its dictionary)
# - an annotation .gtf file exists somewhere
#
# Arguments passed:
# - species type (hek3t3, flies)
# - prefix of the sample
# - number of STAMPs expected
#
# Outputs:
# -

dropSeqDir=/data/rajewsky/shared_bins/Drop-seq_tools-2.3.0

# Setting up species type
if [ $1 == "hek3t3" ]; then
  genome_dir=/data/rajewsky/home/nkarais/hg38_mm10_GRCh38.84_GRCm38.84_STAR_2.4.0/
  annotation_file=/data/rajewsky/home/nkarais/hg38_mm10_GRCh38.84_GRCm38.84_STAR_2.4.0/Homo_sapiens_84_Mus_musculus_84.gtf
  ref_sequence=/data/rajewsky/home/nkarais/hg38_mm10_GRCh38.84_GRCm38.84_STAR_2.4.0/MegaGenome.fa
elif [ $1 == "flies" ]; then
  genome_dir=/data/rajewsky/home/nkarais/dmel_dvir_star/
  annotation_file=${genome_dir}/dmel_dvir_pDsRed_pVm26Aa.gtf
  ref_sequence=/data/rajewsky/home/nkarais/dmel_dvir_star/MegaGenome.fa
elif [ $1 == 'human' ]; then
  genome_dir=/data/rajewsky/home/nkarais/hg38_GRCh38.84_STAR/
  annotation_file=/data/rajewsky/home/nkarais/hg38_GRCh38.84_STAR/Homo_sapiens.GRCh38.84.chr.gtf
  ref_sequence=/data/rajewsky/home/nkarais/hg38_GRCh38.84_STAR/GRCh38.primary_assembly.genome.fa
elif [ $1 == 'mouse' ]; then
  genome_dir=/data/rajewsky/home/nkarais/mm10_GRCm38.p5_STAR/
  annotation_file=/data/rajewsky/home/nkarais/mm10_GRCm38.p5_STAR/gencode.vM12.annotation.gtf
  ref_sequence=/data/rajewsky/home/nkarais/mm10_GRCm38.p5_STAR/GRCm38.primary_assembly.genome.fa
elif [ $1 == 'ercc' ]; then
  genome_dir=/data/rajewsky/home/nkarais/ercc/
  annotation_file=/data/rajewsky/home/nkarais/ercc/ERCC92.gtf
  ref_sequence=/data/rajewsky/home/nkarais/ercc/ERCC92.fa
elif [ $1 == 'ercc_gfp' ]; then
  genome_dir=/data/rajewsky/home/nkarais/ercc_gfp/
  annotation_file=/data/rajewsky/home/nkarais/ercc_gfp/ercc_gfp.gtf
  ref_sequence=/data/rajewsky/home/nkarais/ercc_gfp/ercc_gfp.fa
fi


# The prefix of sample
name=$2
stamps_expected=$3

if [ ! -d ${name} ]; then
	mkdir ${name}
fi

# Convert fastq to unaligned sam
fastq_to_sam="java -jar /data/rajewsky/shared_bins/picard-tools-2.18.17/picard.jar FastqToSam \
F1=${name}_1.fastq.gz \
F2=${name}_2.fastq.gz \
SM=${name} O=/dev/stdout \
SO=queryname \
TMP_DIR=`pwd`/tmp"

tag_cells="${dropSeqDir}/TagBamWithReadSequenceExtended \
SUMMARY=${name}/${name}_summary_unaligned_tag_cell.txt \
BASE_RANGE=1-12 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=false \
TAG_NAME=XC \
NUM_BASES_BELOW_QUALITY=1 \
INPUT=/dev/stdin \
OUTPUT=/dev/stdout \
COMPRESSION_LEVEL=0"

tag_molecular_barcodes="${dropSeqDir}/TagBamWithReadSequenceExtended \
SUMMARY=${name}/${name}_summary_unaligned_tag_molcell.txt \
BASE_RANGE=13-20 \
BASE_QUALITY=10 \
BARCODED_READ=1 \
DISCARD_READ=true \
TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 \
INPUT=/dev/stdin \
OUTPUT=/dev/stdout \
COMPRESSION_LEVEL=0"

filter_bam="${dropSeqDir}/FilterBam \
TAG_REJECT=XQ \
INPUT=/dev/stdin \
OUTPUT=/dev/stdout \
COMPRESSION_LEVEL=0"

trim_SMART="${dropSeqDir}/TrimStartingSequence \
INPUT=/dev/stdin \
OUTPUT=/dev/stdout \
OUTPUT_SUMMARY=${name}/${name}_summary_adapter_trim.txt \
SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
MISMATCHES=0 \
NUM_BASES=5 \
COMPRESSION_LEVEL=0"

trim_polyA="${dropSeqDir}/PolyATrimmer \
INPUT=/dev/stdin \
OUTPUT=${name}/${name}_polyA.bam \
OUTPUT_SUMMARY=${name}/${name}_summary_polyA_trim.txt \
MISMATCHES=0 \
NUM_BASES=6"

bam_to_fastq="java -jar /data/rajewsky/shared_bins/picard-tools-2.18.17/picard.jar SamToFastq \
INPUT=${name}/${name}_polyA.bam \
FASTQ=/dev/stdout \
TMP_DIR=`pwd`/tmp"

star_align="STAR \
--genomeDir ${genome_dir}/STAR_index \
--readFilesIn /dev/stdin \
--outFileNamePrefix ${name}/star_ \
--runThreadN 20"

sort_sam_to_bam_old="samtools sort -n -@ 20 -m 8G -T ./tmp \
${name}/star_Aligned.out.sam \
-o ${name}/star_Aligned_sorted.bam"

sort_sam_to_bam="java -jar /data/rajewsky/shared_bins/picard-tools-2.18.17/picard.jar SortSam \
I=${name}/star_Aligned.out.sam \
O=${name}/star_Aligned_sorted.bam \
SO=queryname"

merge_alignments="java -jar /data/rajewsky/shared_bins/picard-tools-2.18.17/picard.jar MergeBamAlignment \
REFERENCE_SEQUENCE=${ref_sequence} \
UNMAPPED_BAM=${name}/${name}_polyA.bam \
ALIGNED_BAM=${name}/star_Aligned_sorted.bam \
OUTPUT=/dev/stdout \
INCLUDE_SECONDARY_ALIGNMENTS=false \
PAIRED_RUN=false \
TMP_DIR=`pwd`/tmp \
COMPRESSION_LEVEL=0"

# Previous version that only quantified exonic reads
add_gene_annotations_old="${dropSeqDir}/TagReadWithGeneExon \
I=/dev/stdin \
O=${name}/star_gene_exon_tagged.bam \
ANNOTATIONS_FILE=${annotation_file} \
TAG=GE"

# Current version that quantifies functions to include intronic reads as well
add_gene_annotations="${dropSeqDir}/TagReadWithGeneFunction \
I=/dev/stdin \
O=${name}/star_gene_exon_tagged.bam \
ANNOTATIONS_FILE=${annotation_file}"

detect_bead_substitution_errors="${dropSeqDir}/DetectBeadSubstitutionErrors \
I=${name}/star_gene_exon_tagged.bam \
O=${name}/star_gene_exon_tagged_bead_substitution_corrected.bam \
OUTPUT_REPORT=${name}/${name}_substitution_errors_report.txt \
OUTPUT_SUMMARY=${name}/${name}_substitution_errors_summary.txt \
NUM_THREADS=10"

detect_bead_synthesis_errors="${dropSeqDir}/DetectBeadSynthesisErrors \
I=${name}/star_gene_exon_tagged_bead_substitution_corrected.bam \
O=${name}/star_gene_exon_tagged_corrected.bam \
REPORT=${name}/${name}_indel_report.txt \
OUTPUT_STATS=${name}/${name}_synthesis_stats.txt \
SUMMARY=${name}/${name}_synthesis_stats_summary.txt \
PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC \
NUM_THREADS=10"

bam_tag_histogram="${dropSeqDir}/BamTagHistogram \
I=${name}/star_gene_exon_tagged_corrected.bam \
O=${name}/out_readcounts.txt.gz \
TAG=XC"

dge_exon_only="${dropSeqDir}/DigitalExpression \
I=${name}/star_gene_exon_tagged_corrected.bam \
O=${name}/dge.txt.gz \
SUMMARY=${name}/dge_summary.txt \
CELL_BC_FILE=${name}/topBarcodes.txt"

dge_intron_only="${dropSeqDir}/DigitalExpression \
I=${name}/star_gene_exon_tagged_corrected.bam \
O=${name}/dge_intron.txt.gz \
SUMMARY=${name}/dge_intron_summary.txt \
CELL_BC_FILE=${name}/topBarcodes.txt \
LOCUS_FUNCTION_LIST=null \
LOCUS_FUNCTION_LIST=INTRONIC"

dge_exon_intron="${dropSeqDir}/DigitalExpression \
I=${name}/star_gene_exon_tagged_corrected.bam \
O=${name}/dge_all.txt.gz \
SUMMARY=${name}/dge_all_summary.txt \
CELL_BC_FILE=${name}/topBarcodes.txt \
LOCUS_FUNCTION_LIST=INTRONIC"

dgeReads_exon_only="${dropSeqDir}/DigitalExpression \
I=${name}/star_gene_exon_tagged_corrected.bam \
O=${name}/dgeReads.txt.gz \
SUMMARY=${name}/dgeReads_summary.txt \
CELL_BC_FILE=${name}/topBarcodes.txt \
OUTPUT_READS_INSTEAD=true"

dgeReads_intron_only="${dropSeqDir}/DigitalExpression \
I=${name}/star_gene_exon_tagged_corrected.bam \
O=${name}/dgeReads_intron.txt.gz \
SUMMARY=${name}/dgeReads_intron_summary.txt \
CELL_BC_FILE=${name}/topBarcodes.txt \
OUTPUT_READS_INSTEAD=true \
LOCUS_FUNCTION_LIST=null \
LOCUS_FUNCTION_LIST=INTRONIC"

dgeReads_exon_intron="${dropSeqDir}/DigitalExpression \
I=${name}/star_gene_exon_tagged_corrected.bam \
O=${name}/dgeReads_all.txt.gz \
SUMMARY=${name}/dgeReads_all_summary.txt \
CELL_BC_FILE=${name}/topBarcodes.txt \
OUTPUT_READS_INSTEAD=true \
LOCUS_FUNCTION_LIST=INTRONIC"

$fastq_to_sam | \
  $tag_cells | \
  $tag_molecular_barcodes | \
  $filter_bam | \
  $trim_SMART | \
  $trim_polyA

$bam_to_fastq | \
  $star_align

$sort_sam_to_bam

$merge_alignments | \
  $add_gene_annotations

$detect_bead_substitution_errors
$detect_bead_synthesis_errors

$bam_tag_histogram

# Rscript which
# - plots the cuulative fraction of reads over the number of cell barcodes and
# - computes the "elbow" of the graph and writes it in a .txt file
Rscript ~/src/dropseq_pipelines/cell_number.R $name $(($stamps_expected * 10))
cell_number=$(cat ${name}/cell_number.txt)

# Create a txt file with the top $cellnumber barcodes so that there's no
# need to iterate over BAM again
if [ $cell_number == 2 ]; then
	echo "no knee found! calculating gene expression for the top 5,000 barcodes"
	zcat ${name}/out_readcounts.txt.gz | cut -f2 | head -5000 > ${name}/topBarcodes.txt
elif [ $cell_number > 100000 ]; then
  echo "too many bead barcodes found! calculating gene expression for the top 50,000 barcodes"
  zcat ${name}/out_readcounts.txt.gz | cut -f2 | head -50000 > ${name}/topBarcodes.txt
else
	zcat ${name}/out_readcounts.txt.gz | cut -f2 | head -${cell_number} > ${name}/topBarcodes.txt
fi

$dge_exon_only
$dge_intron_only
$dge_exon_intron
$dgeReads_exon_only
$dgeReads_intron_only
$dgeReads_exon_intron

cd ${name}

# determine percentages of intronic/intergenic etc
~/src_murphy/determine_type_of_mapped_reads.sh star_gene_exon_tagged_corrected.bam

# index the bamfile
samtools index star_gene_exon_tagged_corrected.bam

rm ${name}_polyA.bam
rm star_Aligned.out.sam
rm star_Aligned_sorted.bam
rm star_gene_exon_tagged.bam
rm star_gene_exon_tagged_bead_substitution_corrected.bam
