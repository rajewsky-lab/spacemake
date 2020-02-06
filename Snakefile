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
import math

#############
# FUNCTIONS #
#############
def hamming_distance(string1, string2):
    return sum(c1 != c2 for c1, c2 in zip(string1, string2))

def compute_max_barcode_mismatch(indices):
    """computes the maximum number of mismatches allowed for demultiplexing based
    on the indices present in the sample sheet."""
    num_samples = len(indices)
    
    if num_samples == 1:
        return 4
    else:
        max_mismatch = 3
        for i in range(num_samples-1):
            for j in range(i+1, num_samples):
                hd = hamming_distance(indices[i], indices[j])
                max_mismatch = min(max_mismatch, math.ceil(hd/2)-1)
    return max_mismatch

def read_sample_sheet(sample_sheet_path, flowcell_id):
    with open(sample_sheet_path) as sample_sheet:
        ix = 0
        for line in sample_sheet:
            if '[Data]' in line:
                break
            else:
                ix = ix + 1

    df = pd.read_csv(sample_sheet_path, skiprows = ix+1)
    df['species'] = df['Description'].str.split('_').str[1]

    df.rename(columns={"Sample_ID":"sample_id", "Sample_Name":"puck_id", "Sample_Project":"project_id"}, inplace=True)
    
    df['flowcell_id'] = flowcell_id
    df['demux_barcode_mismatch'] = compute_max_barcode_mismatch(df['index'])
    df['sample_sheet'] = sample_sheet_path

    return df[['sample_id', 'puck_id', 'project_id', 'sample_sheet', 'flowcell_id',
               'species', 'demux_barcode_mismatch']]    

def create_lookup_table(df):
    samples_lookup = {}

    projects = df.project_id.unique()

    for p in projects:
        sample_ids = df[df.project_id.eq(p)].sample_id.to_list()
        species = df[df.project_id.eq(p)].species.to_list()
        pucks = df[df.project_id.eq(p)].puck_id.to_list()
        sample_sheet = df[df.project_id.eq(p)].sample_sheet.to_list()[0]
        flowcell_id = df[df.project_id.eq(p)].flowcell_id.to_list()[0]
        demux_barcode_mismatch = df[df.project_id.eq(p)].demux_barcode_mismatch.to_list()[0]

        samples = {}
        for i in range(len(sample_ids)):
            samples[sample_ids[i]] = {
                'species': species[i],
                'puck': pucks[i]
            }

        samples_lookup[p] = {
            'sample_sheet': sample_sheet,
            'flowcell_id': flowcell_id,
            'samples': samples,
            'demux_barcode_mismatch': demux_barcode_mismatch
        }

    return samples_lookup

def get_species_info(wildcards):
    # This function will return 3 things required by STAR:
    #    - annotation (.gtf file)
    #    - genome (.fa file)
    #    - index (a directory where the STAR index is)
    species = samples[wildcards.project]['samples'][wildcards.sample]['species']

    return {
        'annotation': config['knowledge']['annotations'][species],
        'genome': config['knowledge']['genomes'][species],
        'index': config['knowledge']['indices'][species]['star']
    }

def get_dge_extra_params(wildcards):
    dge_type = wildcards.dge_type

    if dge_type == '_exon':
        return ''
    elif dge_type == '_intron':
        return "LOCUS_FUNCTION_LIST=null LOCUS_FUNCTION_LIST=INTRONIC"
    elif dge_type == '_all':
        return "LOCUS_FUNCTION_LIST=INTRONIC"
    if dge_type == 'Reads_exon':
        return "OUTPUT_READS_INSTEAD=true"
    elif dge_type == 'Reads_intron':
        return "OUTPUT_READS_INSTEAD=true LOCUS_FUNCTION_LIST=null LOCUS_FUNCTION_LIST=INTRONIC"
    elif dge_type == 'Reads_all':
        return "OUTPUT_READS_INSTEAD=true LOCUS_FUNCTION_LIST=INTRONIC"
    
####
# this file should contain all sample information, sample name etc.
####
# configfile should be loaded from command line

###############
# Global vars #
###############
# set root output dir
project_dir = '{project}'

illumina_projects = config['illumina_projects']

# get the samples
project_df = pd.concat([read_sample_sheet(ip['sample_sheet'], ip['flowcell_id']) for ip in illumina_projects], ignore_index=True)

samples = create_lookup_table(project_df)
samples_list = project_df.T.to_dict().values()

# create lookup table for flowcell-to-samplesheet
# flowcell_id2samplesheet = {value['flowcell_id'] : value['sample_sheet'] for key, value in samples.items()}

##############
# Demux vars #
##############
# Undetermined files pattern
# they are the output of bcl2fastq, and serve as an indicator to see if the demultiplexing has finished
demux_dir = project_dir + '/demultiplex_data'
demux_indicator = demux_dir + '/indicator.log'

####################################
# FASTQ file linking and reversing #
####################################
reads_suffix = '.fastq.gz'

raw_reads_prefix = project_dir + '/reads/raw/{sample}_R'
raw_reads_pattern = raw_reads_prefix + '{mate}' + reads_suffix
raw_reads_mate_1 = raw_reads_prefix + '1' + reads_suffix
raw_reads_mate_2 = raw_reads_prefix + '2' + reads_suffix

reverse_reads_prefix = project_dir + '/reads/reversed/{sample}_reversed_R'
reverse_reads_pattern = reverse_reads_prefix + '{mate}' + reads_suffix
reverse_reads_mate_1 = reverse_reads_prefix + '1' + reads_suffix
reverse_reads_mate_2 = reverse_reads_prefix + '2' + reads_suffix

###############
# Fastqc vars #
###############
fastqc_root = project_dir + '/reads/fastqc/'
fastqc_pattern = fastqc_root + '{sample}_reversed_R{mate}_fastqc.{ext}'
fastqc_command = '/data/rajewsky/shared_bins/FastQC-0.11.2/fastqc'
fastqc_ext = ['zip', 'html']

########################
# UNIQUE PIPELINE VARS #
########################
# set the tool script directories
picard_tools = '/data/rajewsky/shared_bins/picard-tools-2.21.6/picard.jar'
dropseq_tools = '/data/rajewsky/shared_bins/Drop-seq_tools-2.3.0'

# set per sample vars
dropseq_root = project_dir + '/data/{sample}'

data_root = dropseq_root
dropseq_reports_dir = dropseq_root + '/reports'
dropseq_tmp_dir = dropseq_root + '/tmp'
smart_adapter = config['adapters']['smart']

# file containing R1 and R2 merged
dropseq_merge_in_mate_1 = reverse_reads_mate_1
dropseq_merge_in_mate_2 = reverse_reads_mate_2
dropseq_merged_reads = dropseq_root + '/unaligned.bam'

#######################
# post dropseq and QC #
#######################
qc_sheet_parameters_file = data_root + '/qc_sheet/qc_sheet_parameters.yaml'
qc_sheet = data_root + '/qc_sheet/qc_sheet_{sample}_{puck}.pdf'

# #######################
# include dropseq rules #
# #######################
include: 'dropseq.smk'

################################
# Final output file generation #
################################

def get_final_output_files(pattern, **kwargs):
    out_files = [expand(pattern,
            project=s['project_id'], 
            sample=s['sample_id'],
            puck=s['puck_id'], **kwargs) for s in samples_list]

    out_files = [item for sublist in out_files for item in sublist]
    
    return out_files

#############
# Main rule a
#############
rule all:
    input:
        get_final_output_files(fastqc_pattern, ext = fastqc_ext, mate = [1,2]),
        get_final_output_files(dge_out, dge_type = dge_types),
        get_final_output_files(dropseq_final_bam_ix),
        get_final_output_files(qc_sheet)

rule demultiplex_data:
    params:
        sample_sheet=lambda wildcards: samples[wildcards.project]['sample_sheet'],
        flowcell_id=lambda wildcards: samples[wildcards.project]['flowcell_id'],
        demux_barcode_mismatch=lambda wildcards: samples[wildcards.project]['demux_barcode_mismatch'],
        output_dir=lambda wildcards: expand(demux_dir, project=wildcards.project)
    output:
        demux_indicator
    shell:
        """
        bcl2fastq \
            --no-lane-splitting --fastq-compression-level=9 \
            --mask-short-adapter-reads 15 \
            --barcode-mismatch {params.demux_barcode_mismatch} \
            --output-dir {params.output_dir} \
            --sample-sheet {params.sample_sheet} \
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
        mkdir -p {wildcards.project}/reads/raw

        find {wildcards.project}/demultiplex_data -type f -wholename '*/{wildcards.sample}/*R{wildcards.mate}*.fastq.gz' -exec ln -sr {{}} {output} \; 
        """

rule reverse_first_mate:
    input:
        raw_reads_mate_1
    output:
        reverse_reads_mate_1
    params:
        tmp_file_pattern = lambda wildcards: wildcards.project + '/reads/reversed/' + wildcards.sample + '_small'
    script:
        'reverse_fastq_file.py'

rule reverse_second_mate:
    input:
        raw_reads_mate_2
    output:
        reverse_reads_mate_2
    shell:
        """
        mkdir -p {wildcards.project}/reads/reversed

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
    input:
        samplesheet=lambda wildcards: samples[wildcards.project]['sample_sheet'],
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
