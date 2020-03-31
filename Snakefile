#########
# about #
#########
__version__ = '0.1.1'
__author__ = ['Nikos Karaiskos', 'Tamas Ryszard Sztanka-Toth']
__licence__ = 'GPL'
__email__ = ['nikolaos.karaiskos@mdc-berlin.de', 'tamasryszard.sztanka-toth@mdc-berlin.de']

###########
# imports #
###########
import os
import pandas as pd
import numpy as np
import math

################
# Shell prefix #
################
shell.prefix('set +o pipefail; JAVA_TOOL_OPTIONS="-Xmx8g -Xss2560k" ;')

#############
# FUNCTIONS #
#############
include: 'snakemake_helper_functions.py'

####
# this file should contain all sample information, sample name etc.
####
# configfile should be loaded from command line

###############
# Global vars #
###############
repo_dir = os.path.dirname(workflow.snakefile)

# set root dir where the processed_data goes
project_dir = config['root_dir'] + '/projects/{project}'

microscopy_root = '/data/rajewsky/slideseq_microscopy'
microscopy_raw = microscopy_root + '/raw'
microscopy_qc = microscopy_root + '/qc'

illumina_projects = config['illumina_projects']

# get the samples
project_df = pd.concat([read_sample_sheet(ip['sample_sheet'], ip['flowcell_id']) for ip in illumina_projects], ignore_index=True)

# add additional samples from config.yaml, which have already been demultiplexed. add none instead of NaN
project_df = project_df.append(config['additional_illumina_projects'], ignore_index=True).replace(np.nan, 'none', regex=True)

samples = create_lookup_table(project_df)
samples_list = project_df.T.to_dict().values()

project_puck_df = project_df.merge(get_sample_info(microscopy_raw), how='inner', on ='puck_id')

demux_dir2project = {s['demux_dir']: s['project_id'] for s in samples_list}

# create lookup table for flowcell-to-samplesheet
# flowcell_id2samplesheet = {value['flowcell_id'] : value['sample_sheet'] for key, value in samples.items()}
#################
# DIRECTORY STR #
#################
raw_data_root = project_dir + '/raw_data'
raw_data_illumina = raw_data_root + '/illumina'
raw_data_illumina_reads = raw_data_illumina + '/reads/raw'
raw_data_illumina_reads_reversed = raw_data_illumina + '/reads/reversed'

raw_data_optical = raw_data_root + '/optical'
raw_data_optical_images = raw_data_optical + '/{sample}'

processed_data_root = project_dir + '/processed_data/{sample}'
processed_data_illumina = processed_data_root + '/illumina'

processed_data_optical = processed_data_root + '/optical'

# metadata file created during the linking rule
projects_puck_info = config['root_dir'] + '/.config/projects_puck_info.csv'


##############
# Demux vars #
##############
# Undetermined files pattern
# they are the output of bcl2fastq, and serve as an indicator to see if the demultiplexing has finished
demux_dir_pattern = config['root_dir'] + '/demultiplex_data/{demux_dir}'
demux_indicator = demux_dir_pattern + '/indicator.log'

####################################
# FASTQ file linking and reversing #
####################################
reads_suffix = '.fastq.gz'

raw_reads_prefix = raw_data_illumina_reads + '/{sample}_R'
raw_reads_pattern = raw_reads_prefix + '{mate}' + reads_suffix
raw_reads_mate_1 = raw_reads_prefix + '1' + reads_suffix
raw_reads_mate_2 = raw_reads_prefix + '2' + reads_suffix

reverse_reads_prefix = raw_data_illumina_reads_reversed + '/{sample}_reversed_R'
reverse_reads_pattern = reverse_reads_prefix + '{mate}' + reads_suffix
reverse_reads_mate_1 = reverse_reads_prefix + '1' + reads_suffix
reverse_reads_mate_2 = reverse_reads_prefix + '2' + reads_suffix

###############
# Fastqc vars #
###############
fastqc_root = raw_data_illumina + '/fastqc'
fastqc_pattern = fastqc_root + '/{sample}_reversed_R{mate}_fastqc.{ext}'
fastqc_command = '/data/rajewsky/shared_bins/FastQC-0.11.2/fastqc'
fastqc_ext = ['zip', 'html']

########################
# UNIQUE PIPELINE VARS #
########################
# set the tool script directories
picard_tools = '/data/rajewsky/shared_bins/picard-tools-2.21.6/picard.jar'
dropseq_tools = '/data/rajewsky/shared_bins/Drop-seq_tools-2.3.0'

# set per sample vars
dropseq_root = processed_data_illumina + '/complete_data'

data_root = dropseq_root
dropseq_reports_dir = dropseq_root + '/reports'
dropseq_tmp_dir = dropseq_root + '/tmp'
smart_adapter = config['adapters']['smart']

# subsample vars
downsample_root = processed_data_illumina + '/downsampled_data'

# file containing R1 and R2 merged
dropseq_merge_in_mate_1 = reverse_reads_mate_1
dropseq_merge_in_mate_2 = reverse_reads_mate_2
dropseq_merged_reads = dropseq_root + '/unaligned.bam'

#######################
# post dropseq and QC #
#######################
qc_sheet_parameters_file = data_root + '/qc_sheet/qc_sheet_parameters.yaml'
qc_sheet = data_root + '/qc_sheet/qc_sheet_{sample}_{puck}.pdf'

# reads type
reads_type_out = dropseq_root + '/uniquely_mapped_reads_type.txt'

############################
# Link optical to illumina #
############################
optical_linked = []

for i, row in project_puck_df.iterrows():
    optical_linked = optical_linked + expand(raw_data_optical_images,
                                       project = row['project_id'],
                                       sample = row['sample_id'])
    
    optical_linked = optical_linked + expand(processed_data_optical,
                                       project = row['project_id'],
                                       sample = row['sample_id'])

# #######################
# include dropseq rules #
# #######################
include: 'dropseq.smk'

################################
# Final output file generation #
################################

def get_final_output_files(pattern, projects = 'all', **kwargs):
    if projects == 'all':
        samples = samples_list
    else:
        samples = [s for s in samples_list if s['project_id'] in projects]

    out_files = [expand(pattern,
            project=s['project_id'], 
            sample=s['sample_id'],
            puck=s['puck_id'], **kwargs) for s in samples]

    out_files = [item for sublist in out_files for item in sublist]
    
    return out_files

#############
# Main rule #
#############
rule all:
    input:
        #get_final_output_files(dge_out, dge_type = dge_types),
        get_final_output_files(dropseq_final_bam_ix),
        get_final_output_files(qc_sheet),
        get_final_output_files(fastqc_pattern, ext = fastqc_ext, mate = [1,2])

###################
# LINK TO OPTICAL #
###################
rule link_optical:
    input:
        optical_linked

########################
# CREATE METADATA FILE #
########################
rule create_projects_puck_info:
    input:
        projects_puck_info

################
# DOWNSAMPLING #
################
include: 'downsample.smk'

rule downsample:
    input:
        get_final_output_files(downsample_saturation_analysis, projects = config['downsample_projects'])

#################
# MERGE SAMPLES #
#################
include: 'merge_samples.smk'

samples_to_merge = []

# expect a list of lists in the config file. samples in each list will be merged
if 'samples_to_merge' in config:
    for merge_group in config['samples_to_merge']
        samples_to_merge = expand(merged_qc_sheet, merged_name = '.'.join(merge_group))

rule merge_samples:
    input:
        expand(merged_qc_sheet, merged_name = samples_to_merge)

#########
# RULES #
#########
ruleorder: link_raw_reads > link_demultiplexed_reads 

rule demultiplex_data:
    params:
        demux_barcode_mismatch=lambda wildcards: samples[demux_dir2project[wildcards.demux_dir]]['demux_barcode_mismatch'],
        sample_sheet=lambda wildcards: samples[demux_dir2project[wildcards.demux_dir]]['sample_sheet'],
        flowcell_id=lambda wildcards: samples[demux_dir2project[wildcards.demux_dir]]['flowcell_id'],
        output_dir= lambda wildcards: expand(demux_dir_pattern, demux_dir=wildcards.demux_dir)
    input:
        unpack(get_basecalls_dir)
    output:
        demux_indicator
    threads: 16
    shell:
        """
        bcl2fastq \
            --no-lane-splitting --fastq-compression-level=9 \
            --mask-short-adapter-reads 15 \
            --barcode-mismatch {params.demux_barcode_mismatch} \
            --output-dir {params.output_dir} \
            --sample-sheet {params.sample_sheet} \
            --runfolder-dir  {input} \
            -r {threads} -p {threads} -w {threads}
            
            echo "demux finished: $(date)" > {output}
        """

rule link_demultiplexed_reads:
    input:
        ancient(unpack(get_demux_indicator))
    output:
        raw_reads_pattern
    params:
        demux_dir = lambda wildcards: expand(demux_dir_pattern, demux_dir=get_demux_dir(wildcards)),
        reads_folder = raw_data_illumina_reads
    shell:
        """
        mkdir -p {params.reads_folder}

        find {params.demux_dir} -type f -wholename '*/{wildcards.sample}/*R{wildcards.mate}*.fastq.gz' -exec ln -sr {{}} {output} \; 
        """

def get_reads(wildcards):
    return [samples[wildcards.project]['samples'][wildcards.sample]['R'+wildcards.mate]]

rule link_raw_reads:
    input:
        unpack(get_reads)
    output:
        raw_reads_pattern
    shell:
        """
        ln -s {input} {output}
        """

rule reverse_first_mate:
    input:
        raw_reads_mate_1
    output:
        reverse_reads_mate_1
    script:
        'reverse_fastq_file.py'

rule reverse_second_mate:
    input:
        raw_reads_mate_2
    output:
        reverse_reads_mate_2
    params:
        reads_folder = raw_data_illumina_reads_reversed
    shell:
        """
        mkdir -p {params.reads_folder}

        ln -sr {input} {output}
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

rule determine_precentages:
    input:
        dropseq_final_bam
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
        dropseq_final_bam
    output:
        dropseq_final_bam_ix 
    shell:
       "samtools index {input}"

rule create_qc_parameters:
    params:
        sample_id = lambda wildcards: wildcards.sample,
        project_id = lambda wildcards: wildcards.project,
        puck_id = lambda wildcards: samples[wildcards.project]['samples'][wildcards.sample]['puck'],
        experiment = lambda wildcards: samples[wildcards.project]['samples'][wildcards.sample]['experiment'],
        sequencing_date = lambda wildcards: samples[wildcards.project]['sequencing_date'],
        input_beads = '60k-100k',
        threshold= '100'
    output:
        qc_sheet_parameters_file
    script:
        "qc_sequencing_create_parameters_from_sample_sheet.py"

rule create_qc_sheet:
    input:
        star_log = star_log_file,
        reads_type_out=reads_type_out,
        parameters_file=qc_sheet_parameters_file,
        read_counts = dropseq_out_readcounts,
        dge_all_summary = dge_root + '/dge_all_summary.txt'
    output:
        qc_sheet
    script:
        "qc_sequencing_create_sheet.py"

rule link_raw_data_images:
    input:
        unpack(get_raw_data_optical_images_input)
    output:
        directory(raw_data_optical_images)
    shell:
        """
        ln -s {input} {output}
        """

rule linked_processed_data_optical:
    input:
        unpack(get_processed_data_optical)
    output:
        directory(processed_data_optical)
    shell:
        """
        ln -s {input} {output}
        """

rule create_projects_metadata:
    output:
        projects_puck_info
    run:
        project_puck_df.to_csv(output[0], index=False) 
