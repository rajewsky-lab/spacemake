# snakemake root file to create pipeline for the spatial sequencing illumina data
#
# author: tsztank
# email: tamasryszard.sztanka-toth@mdc-berlin.de
#
# ###

####
# import necessary python packages
####
import os
from sample_sheet import SampleSheet

def get_samples(runs):
    out_dict = {}
    for key in runs.keys():
        out_dict[key] = [sample['Sample_ID'] for sample in SampleSheet(runs[key]['samplesheet']).samples]

    return(out_dict)

####
# this file should contain all sample information, sample name etc.
####
configfile: 'config.yaml'

###############
# Global vars #
###############
# set root output dir
root_dir = 'sequencing_runs'

# runs as a dictionary
runs = config['illumina_runs']

# get the samples
samples = get_samples(runs) 

##############
# Demux vars #
##############
# demultiplexing root. for each run we demultiplex only once
demux_out = root_dir + '/{run}/demux_data'

# Undetermined files pattern
# they are the output of bcl2fastq, and serve as an indicator to see if the demultiplexing has finished
demux_indicator = root_dir + '/{run}/demux_data/indicator.log'

####################################
# FASTQ file linking and reversing #
####################################
reads_suffix = '.fastq.gz'

raw_reads_prefix = root_dir + '/{run}/reads/raw/{sample}_R'
raw_reads_pattern = raw_reads_prefix + '{mate}' + reads_suffix
raw_reads_mate_1 = raw_reads_prefix + '1' + reads_suffix
raw_reads_mate_2 = raw_reads_prefix + '2' + reads_suffix

reverse_reads_prefix = root_dir + '/{run}/reads/reversed/{sample}_reversed_R'
reverse_reads_pattern = reverse_reads_prefix + '{mate}' + reads_suffix
reverse_reads_mate_1 = reverse_reads_prefix + '1' + reads_suffix
reverse_reads_mate_2 = reverse_reads_prefix + '2' + reads_suffix

###############
# Fastqc vars #
###############
fastqc_root = root_dir + '/{run}/fastqc/'
fastqc_pattern = fastqc_root + '{sample}_reversed_R{mate}_fastqc.{ext}'
fastqc_command = '/data/rajewsky/shared_bins/FastQC-0.11.2/fastqc'

# create fastqc out list
fastqc_files_out = [expand(fastqc_pattern, run=key, sample=value, mate=[1,2], ext=['html', 'zip']) for key, value in samples.items()]

# flatten the list
fastqc_files_out = [item for sublist in fastqc_files_out for item in sublist]

#########################
# Dropseq pipeline vars #
#########################
# set the tool script directories
picard_tools = '/data/rajewsky/shared_bins/picard-tools-2.21.6/picard.jar'
dropseq_tools = '/data/rajewsky/shared_bins/Drop-seq_tools-2.3.0'

# set per sample vars
dropseq_root = root_dir + '/{run}/data/{sample}'
data_root = dropseq_root
dropseq_reports_dir = dropseq_root + '/reports'
dropseq_tmp_dir = dropseq_root + '/tmp'
smart_adapter = config['adapters']['smart']

# file containing R1 and R2 merged
dropseq_merged_reads = dropseq_root + '/unaligned.bam'

# tag reads with umis and cells
dropseq_cell_tagged = dropseq_root + '/unaligned_tagged_umi_cell.bam'
dropseq_umi_tagged = dropseq_root + '/unaligned_tagged_umi.bam'

# filter out XC tag
dropseq_tagged_filtered = dropseq_root + '/unaligned_tagged_filtered.bam'

# trim smart adapter from the reads
dropseq_tagged_filtered_trimmed = dropseq_root + '/unaligned_tagged_filtered_trimmed.bam'

# trim polyA overheang if exists
dropseq_tagged_filtered_trimmed_polyA = dropseq_root + '/unaligned_tagged_filtered_trimmed_polyA.bam'

# create fastq file from the previous .bam to input into STAR
dropseq_star_input = dropseq_root + '/unaligned_reads_star_input.fastq'

# mapped reads
dropseq_mapped_reads = dropseq_root + '/star_Aligned.out.sam'
star_log_file = dropseq_root + '/star_Log.final.out'

# sort reads and create bam
dropseq_mapped_sorted_reads = dropseq_root + '/star_Aligned.sorted.bam'

# merge bam files
dropseq_merged = dropseq_root + '/merged.bam'

# tag gene with exon
dropseq_gene_exon_tagged = dropseq_root + '/star_gene_exon_tagged.bam'

# detect bead substitution errors
dropseq_bead_substitution_cleaned = dropseq_root + '/clean_substitution.bam'

# detect bead synthesis errors
dropseq_mapped_clean_reads = dropseq_root + '/clean.bam'
synthesis_stats_summary = dropseq_reports_dir + '/detect_bead_synthesis_error.stats.txt'
substitution_error_report = dropseq_reports_dir + '/detect_bead_substitution_error.report.txt'

# index bam file
dropseq_mapped_clean_reads_ix = dropseq_mapped_clean_reads + '.bai'

# create readcounts file
dropseq_out_readcounts = dropseq_root + '/out_readcounts.txt.gz'

# create a file with the top barcodes
dropseq_top_barcodes = dropseq_root + '/topBarcodes.txt'

# dges
dge_root = dropseq_root + '/dge'
dge_exon_only =     dge_root + '/dge.txt.gz'
dge_intron_only =   dge_root + '/dge_intron.txt.gz'
dge_exon_intron =     dge_root + '/dge_all.txt.gz'

dgeReads_exon_only =   dge_root + '/dgeReads.txt.gz'
dgeReads_intron_only =   dge_root + '/dgeReads_intron.txt.gz'
dgeReads_exon_intron =     dge_root + '/dgeReads_all.txt.gz'

#######################
# post dropseq and QC #
#######################
reads_type_out = dropseq_root + '/uniquely_mapped_reads_type.txt'
cell_number = data_root + '/cell_number.txt'
cell_cummulative_plot = data_root + '/cell_cummulative.png'
qc_sheet_parameters_file = data_root + '/qc_sheet_parameters.yaml'
qc_sheet = data_root + '/qc_sheet.pdf'

################################
# Final output file generation #
################################

def get_final_output_files(pattern):
    # create fastqc out list
    out_files = [expand(pattern, run=key, sample=value, mate=[1,2]) for key, value in samples.items()]

    # flatten the list
    out_files = [item for sublist in out_files for item in sublist]

    return out_files

#############
# Main rule a
#############
rule all:
    input:
       get_final_output_files(dge_exon_only),
       get_final_output_files(reads_type_out),
       get_final_output_files(cell_number),
       get_final_output_files(cell_cummulative_plot),
       get_final_output_files(qc_sheet_parameters_file),
       get_final_output_files(dropseq_mapped_clean_reads_ix)

rule demultiplex_data:
    params:
        samplesheet=lambda wildcards: runs[wildcards.run]['samplesheet'],
        flowcell_id=lambda wildcards: runs[wildcards.run]['flowcell_id'],
        output_dir= lambda wildcards: expand(demux_out, run=wildcards.run)
    output:
        demux_indicator
    shell:
        """
        bcl2fastq \
            --no-lane-splitting --fastq-compression-level=9 \
            --mask-short-adapter-reads 15 \
            --barcode-mismatch 1 \
            --output-dir {params.output_dir} \
            --sample-sheet {params.samplesheet} \
            --runfolder-dir /data/remote/basecalls/{params.flowcell_id}

        echo "demux finished: $(date)" > {output}
        """

rule link_raw_reads:
    output:
        raw_reads_pattern
    input:
        demux_indicator
    # isntead of hard links the link is now relative 
    shell:
        """
        mkdir -p {root_dir}/{wildcards.run}/reads/raw

        find {root_dir}/{wildcards.run}/demux_data -type f -wholename '*/{wildcards.sample}/*R{wildcards.mate}*.fastq.gz' -exec ln -s ../../../../{{}} {output} \; 
        """

rule reverse_first_mate:
    input:
        raw_reads_mate_1
    output:
        reverse_reads_mate_1
    params:
        tmp_file_pattern = lambda wildcards: root_dir + '/' + wildcards.run + '/reads/reversed/' + wildcards.sample + '_small'
    script:
        'reverse_fastq_file.py'

rule reverse_second_mate:
    input:
        raw_reads_mate_2
    output:
        reverse_reads_mate_2
    shell:
        """
        mkdir -p {root_dir}/{wildcards.run}/reads/reversed

        ln -s ../../../../{input} {output}
        """

rule run_fastqc:
    input:
        reverse_reads_pattern
    output:
        fastqc_pattern
    params:
        output_dir = fastqc_root 
    threads: 8
    shell:
        """
        mkdir -p {params.output_dir}

        {fastqc_command} -t {threads} -o {params.output_dir} {input}
        """
# #######################
# include dropseq rules #
# #######################
include: 'dropseq.smk'


rule determine_precentages:
    input:
        dropseq_mapped_clean_reads
    output:
        reads_type_out
    shell:
        ## Script taken from sequencing_analysis.sh
        """
        samtools view {input} | \
          awk '!/GE:Z:/ && $5 == "255" && match ($0, "XF:Z:") split(substr($0, RSTART+5), a, "\t") {{print a[1]}}' | \
          awk 'BEGIN {{ split("INTRONIC INTERGENIC CODING UTR", keyword)
                      for (i in keyword) count[keyword[i]]=0
                    }}
              /INTRONIC/  {{ count["INTRONIC"]++ }}
              /INTERGENIC/  {{ count["INTERGENIC"]++ }}
              /CODING/ {{count["CODING"]++ }}
              /UTR/ {{ count["UTR"]++ }}
              END   {{
                      for (i in keyword) print keyword[i], count[keyword[i]]
                    }}' > {output}
        """

rule index_bam_file:
    input:
        dropseq_mapped_clean_reads
    output:
        dropseq_mapped_clean_reads_ix 
    shell:
       "samtools index {input}"
rule estimate_cell_number:
    input:
        dropseq_out_readcounts = dropseq_root + '/out_readcounts.txt.gz'
    output:
        cell_number=cell_number,
        cummulative_plot = cell_cummulative_plot
    params:
        knee_limit = 250000
    script:
        "estimate_cell_number.R"

rule create_qc_parameters:
    input:
        samplesheet=lambda wildcards: runs[wildcards.run]['samplesheet'],
    output:
        qc_sheet_parameters_file
    script:
        "qc_sequencing_create_parameters_from_sample_sheet.py"

rule create_qc_sheet:
    input:
        star_log = star_log_file,
        reads_type_out=reads_type_out,
        synthesys_stats_summary=synthesis_stats_summary,
        parameters_file=qc_sheet_parameters_file
    output:
        qc_sheet 
    script:
        "qc_sequencing_create_sheet.py"
