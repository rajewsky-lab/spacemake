#bowtie2 -x /data/rajewsky/indices/hsa_rrna_ncbi_bowtie_2.3.2/hsa_rrna_ncbi -U /data/rajewsky/projects/slide_seq/projects/sts_043/raw_data/illumina/reads/raw/sts_043_3_R2.fastq.gz -p 20 --very-fast-local --un-gz out/sts_043_3_RiboDepleted.fastq.gz > /dev/null 2> out/sts_043_3_ribo_log.txt

#bowtie2 -x /data/rajewsky/indices/hsa_rrna_ncbi_bowtie_2.3.2/hsa_rrna_ncbi -U /data/rajewsky/projects/slide_seq/projects/sts_043/raw_data/illumina/reads/raw/sts_043_5_R2.fastq.gz -p 20 --very-fast-local --un-gz out/sts_043_5_RiboDepleted.fastq.gz > /dev/null 2> out/sts_043_5_ribo_log.txt

#bowtie2 -x /data/rajewsky/indices/hsa_rrna_ncbi_bowtie_2.3.2/hsa_rrna_ncbi -U /data/rajewsky/projects/slide_seq/projects/sts_048/raw_data/illumina/reads/raw/sts_048_2_R2.fastq.gz -p 20 --very-fast-local --un-gz out/sts_048_2_RiboDepleted.fastq.gz > /dev/null 2> out/sts_048_2_ribo_log.txt

#bowtie2 -x /data/rajewsky/indices/mm10_rRNA_bowtie2_2.3.3.1/mouse_rRNA -U /data/rajewsky/projects/slide_seq/projects/sts_048/raw_data/illumina/reads/raw/sts_048_4_R2.fastq.gz -p 20 --very-fast-local --un-gz out/sts_048_4_RiboDepleted.fastq.gz > /dev/null 2> out/sts_048_4_ribo_log.txt

bowtie2 -x /data/rajewsky/indices/hsa_rrna_ncbi_bowtie_2.3.2/hsa_rrna_ncbi -U /data/rajewsky/sequencing/human/dropseq_runs/ds_046_HeLa_20170707/ds_046_2.fastq.gz -p 20 --very-fast-local --un-gz out/ds_046_HeLa_RiboDepleted.fastq.gz > /dev/null 2> out/ds_046_HeLa_ribo_log.txt

bowtie2 -x /data/rajewsky/indices/hsa_rrna_ncbi_bowtie_2.3.2/hsa_rrna_ncbi -U /data/rajewsky/sequencing/human/dropseq_runs/og_009_1234_20180904/ARybak-Wolf/NR_AR_047/SVS1C_1_S3_R2_001.fastq.gz -p 20 --very-fast-local --un-gz out/NR_AR_047_RiboDepleted.fastq.gz > /dev/null 2> out/NR_AR_047_ribo_log.txt

bowtie2 -x /data/rajewsky/indices/hsa_rrna_ncbi_bowtie_2.3.2/hsa_rrna_ncbi -U /data/rajewsky/sequencing/human/dropseq_runs/ds_041_043_muscle_satellite_cells/ds_041_rep1_2.fastq.gz -p 20 --very-fast-local --un-gz out/ds_041_muscle_cells_rep1_RiboDepleted.fastq.gz > /dev/null 2> out/ds_041_muscle_cells_rep1_ribo_log.txt

bowtie2 -x /data/rajewsky/indices/mm10_rRNA_bowtie2_2.3.3.1/mouse_rRNA -U /data/rajewsky/sequencing/mouse/kidney_glomeruli/ds_040_042_045_20170601/ds_042_2.fastq.gz -p 20 --very-fast-local --un-gz out/ds_042_mouse_kidney_RiboDepleted.fastq.gz > /dev/null 2> out/ds_042_mouse_kidney_ribo_log.txt
