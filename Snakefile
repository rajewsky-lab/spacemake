#########
# about #
#########
__version__ = '0.1.1'
__author__ = ['Nikos Karaiskos', 'Tamas Ryszard Sztanka-Toth']
__license__ = 'GPL'
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
shell.prefix('set +o pipefail; JAVA_TOOL_OPTIONS="-Xmx8g -Xss2560k" ; umask g+w; ')

#############
# FUNCTIONS #
#############
include: 'snakemake/snakemake_helper_functions.py'

####
# this file should contain all sample information, sample name etc.
####
# configfile should be loaded from command line

###############
# Global vars #
###############
temp_dir = config['temp_dir']
repo_dir = os.path.dirname(workflow.snakefile)
# create puck_data root directory from pattern
config['puck_data']['root'] = config['puck_data']['root'].format(root_dir = config['root_dir'])

# set root dir where the processed_data goes
project_dir = config['root_dir'] + '/projects/{project}'

projects = config.get('projects', None)

# moved barcode_flavor assignment here so that additional samples/projects are equally processed
project_df = create_project_df()

project_df = df_assign_bc_flavor(project_df)

#################
# DIRECTORY STR #
#################
raw_data_root = project_dir + '/raw_data'
raw_data_illumina = raw_data_root + '/illumina'
raw_data_illumina_reads = raw_data_illumina + '/reads/raw'
raw_data_illumina_reads_reversed = raw_data_illumina + '/reads/reversed'
processed_data_root = project_dir + '/processed_data/{sample}'
processed_data_illumina = processed_data_root + '/illumina'

reports_root = config['root_dir'] + '/reports'
project_df_file = reports_root + '/project_df.csv'
sample_overview_file = reports_root + '/sample_overview.html'
sample_read_metrics_db = reports_root + '/sample_read_metrics_db.tsv'

united_illumina_root = config['root_dir'] + '/projects/{united_project}/processed_data/{united_sample}/illumina'
united_complete_data_root = united_illumina_root + '/complete_data'

##############
# Demux vars #
##############
# Undetermined files pattern
# they are the output of bcl2fastq, and serve as an indicator to see if the demultiplexing has finished
demux_dir_pattern = config['root_dir'] + '/raw_data/demultiplex_data/{demux_dir}'
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
reverse_reads_mate_1 = reverse_reads_prefix + '1' + reads_suffix

###############
# Fastqc vars #
###############
bin_dir = config['bin_dir']
fastqc_root = raw_data_illumina + '/fastqc'
fastqc_pattern = fastqc_root + '/{sample}_R{mate}_fastqc.{ext}'
fastqc_command = f'{bin_dir}/FastQC-0.11.2/fastqc'
fastqc_ext = ['zip', 'html']

########################
# UNIQUE PIPELINE VARS #
########################
# set the tool script directories
picard_tools = f'{bin_dir}/picard-tools-2.21.6/picard.jar'
dropseq_tools = f'{bin_dir}/Drop-seq_tools-2.3.0'

# set per sample vars
dropseq_root = processed_data_illumina + '/complete_data'

data_root = dropseq_root
dropseq_reports_dir = dropseq_root + '/reports'
dropseq_tmp_dir = dropseq_root + '/tmp'
smart_adapter = config['adapters']['smart']

dropseq_merged_reads = dropseq_root + '/unaligned.bam'

###
# splitting the reads
###

united_split_reads_root = united_complete_data_root + '/split_reads/'
united_unmapped_bam = united_split_reads_root + 'unmapped.bam'

united_split_reads_sam_names = ['plus_plus', 'plus_minus', 'minus_minus', 'minus_plus', 'plus_AMB', 'minus_AMB']
united_split_reads_sam_pattern = united_split_reads_root + '{file_name}.sam'
united_split_reads_bam_pattern = united_split_reads_root + '{file_name}.bam'

united_split_reads_sam_files = [united_split_reads_root + x + '.sam' for x in united_split_reads_sam_names]

united_split_reads_strand_type = united_split_reads_root + 'strand_type_num.txt'
united_split_reads_read_type = united_split_reads_root + 'read_type_num.txt'

#######################
# post dropseq and QC #
#######################
# qc generation for ALL samples, merged and non-merged
# umi cutoffs. used by qc-s and automated reports
umi_cutoffs = [10, 50, 100]

# parameters file for not merged samples
qc_sheet_parameters_file = data_root + '/qc_sheet_parameters.yaml'

united_qc_sheet = united_complete_data_root +'/qc_sheet_{united_sample}_{puck}.html'
united_star_log = united_complete_data_root + '/star_Log.final.out'
united_reads_type_out = united_split_reads_read_type
united_qc_sheet_parameters_file = united_complete_data_root + '/qc_sheet_parameters.yaml'
united_barcode_readcounts = united_complete_data_root + '/out_readcounts.txt.gz'
united_strand_info = united_split_reads_strand_type

# united final.bam
united_final_bam = united_complete_data_root + '/final.bam'
united_top_barcodes = united_complete_data_root + '/topBarcodes.txt'
united_top_barcodes_clean = united_complete_data_root + '/topBarcodesClean.txt'

# united dge
dge_root = united_complete_data_root + '/dge'
dge_out_prefix = dge_root + '/dge{dge_type}{dge_cleaned}'
dge_out = dge_out_prefix + '.txt.gz'
dge_out_summary = dge_out_prefix + '_summary.txt'
dge_types = ['_exon', '_intron', '_all', 'Reads_exon', 'Reads_intron', 'Reads_all']

dge_all_summary = united_complete_data_root +  '/dge/dge_all_summary.txt'
dge_all_cleaned_summary = united_complete_data_root +  '/dge/dge_all_cleaned_summary.txt'
dge_all_summary_fasta= united_complete_data_root +  '/dge/dge_all_summary.fa'
dge_all = united_complete_data_root +  '/dge/dge_all.txt.gz'
dge_all_cleaned = united_complete_data_root +  '/dge/dge_all_cleaned.txt.gz'

# primer analysis
united_primer_tagged_final_bam = united_complete_data_root +  '/primer_tagged_final.bam'
united_primer_tagged_summary = united_complete_data_root + '/primer_tagged_summary.txt'

# kmer stats per position
kmer_stats_file = data_root + '/kmer_stats/{kmer_len}mer_counts.csv'

# map paired-end to check errors
paired_end_prefix = data_root + '/mapped_paired_end/'
paired_end_sam = paired_end_prefix + '{sample}_paired_end.sam'
paired_end_bam = paired_end_prefix + '{sample}_paired_end.bam'
paired_end_flagstat = paired_end_prefix + '{sample}_paired_end_flagstat.txt'
paired_end_log = paired_end_prefix + '{sample}_paired_end.log'
paired_end_mapping_stats = paired_end_prefix + '{sample}_paired_end_mapping_stats.txt'

# automated analysis
automated_analysis_root = united_complete_data_root + '/automated_analysis/umi_cutoff_{umi_cutoff}'
automated_report = automated_analysis_root + '/{united_sample}_{puck}_illumina_automated_report.html'

automated_result_files = {
    'res_file': '/results.h5ad',
    'cluster_markers': '/top10_cluster_markers.csv',
    'obs_df': '/obs_df.csv',
    'long_expr_df': '/long_expr_df.csv'
    }

# prepend automated_result_root
automated_result_files = {key: automated_analysis_root + value for key, value in automated_result_files.items()}

# blast out
blast_db_primers = repo_dir + '/sequences/primers.fa'
blast_db_primers_files = [blast_db_primers + '.' + x for x in ['nhr', 'nin', 'nog', 'nsd', 'nsi', 'nsq']]
blast_header_out = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore"
united_barcode_blast_out = united_complete_data_root + '/cell_barcode_primer_blast_out.txt'

# downsample vars
downsample_root = united_illumina_root + '/downsampled_data'

# in silico repo depletion
ribo_depletion_log = data_root + '/ribo_depletion_log.txt'
united_ribo_depletion_log = united_complete_data_root + '/ribo_depletion_log.txt'
united_parsed_ribo_depletion_log = united_complete_data_root + '/parsed_ribo_depletion_log.txt'

# global wildcard constraints
wildcard_constraints:
    sample='(?!merged_).+',
    project='(?!merged_).+',
    dge_cleaned='|_cleaned',
    dge_type = '|'.join(dge_types),
    pacbio_ext = 'fq|fastq'

# #########################
#  dropseq rules and vars #
# #########################
dropseq_tagged = dropseq_root + '/unaligned_bc_tagged.bam'
dropseq_unassigned = dropseq_root + '/unaligned_bc_unassigned.bam'

# filter out XC tag
dropseq_tagged_filtered = dropseq_root + '/unaligned_tagged_filtered.bam'

# trim smart adapter from the reads
dropseq_tagged_trimmed = dropseq_root + '/unaligned_tagged_trimmed.bam'

# trim polyA overheang if exists
dropseq_tagged_trimmed_polyA = dropseq_root + '/unaligned_tagged_trimmed_polyA.bam'

# mapped reads
dropseq_mapped_reads_unsorted_headerless = dropseq_root + '/star_Aligned.unsorted.headerless.out.bam'
dropseq_mapped_reads_unsorted = dropseq_root + '/star_Aligned.unsorted.out.bam'
dropseq_mapped_reads = dropseq_root + '/star_Aligned.sorted.out.bam'
star_log_file = dropseq_root + '/star_Log.final.out'

# final dropseq bfinal dropseq bam
dropseq_final_bam = dropseq_root + '/final.bam'

# index bam file
dropseq_final_bam_ix = dropseq_final_bam + '.bai'

# include dropseq
include: 'snakemake/dropseq.smk'

################################
# Final output file generation #
################################
def get_final_output_files(pattern, projects = None, samples = None, **kwargs):
    samples_to_run = project_df.T.to_dict().values()

    if projects is not None:
        samples_to_run = [s for s in samples_to_run if s['project_id'] in projects]

    if samples is not None:
        samples_to_run = [s for s in samples_to_run if s['sample_id'] in samples]

    out_files = [expand(pattern,
            project=s['project_id'], 
            sample=s['sample_id'],
            puck=s['puck_id'], **kwargs) for s in samples_to_run]

    out_files = [item for sublist in out_files for item in sublist]
    
    return out_files

def get_united_output_files(pattern, projects = None, samples = None,
                            skip_projects = None, skip_samples = None, **kwargs):
    out_files = []
    df = project_df

    if projects is None and samples is None:
        # nothing is set, assume we need to use all samples
        pass
    elif samples is None and projects is not None:
        df = df[df.project_id.isin(projects)]
    elif samples is not None and projects is None:
        df = df[df.sample_id.isin(samples)]
    else:
        # both are set
        df = df[df.sample_id.isin(samples) | df.project_id.isin(projects)]


    if skip_samples is not None:
        df = df[~df.sample_id.isin(skip_samples)]

    if skip_projects is not None:
        df = df[~df.project_id.isin(skip_projects)]

    for index, row in df.iterrows():
        umi_cutoff = get_downstream_analysis_variables(project_id = row['project_id'],
            sample_id = row['sample_id'])['umi_cutoff']

        out_files = out_files + expand(pattern,
            united_project = row['project_id'],
            united_sample = row['sample_id'],
            puck=row['puck_id'], 
            umi_cutoff = umi_cutoff,
            **kwargs)

    return out_files

human_mouse_samples = project_df[project_df.species.isin(['human', 'mouse'])].sample_id.to_list()

##################
# include pacbio #
##################
processed_data_pacbio = processed_data_root + '/pacbio'
pacbio_fq = raw_data_root + '/pacbio/{sample}.{pacbio_ext}'
pacbio_report = processed_data_pacbio + '/{sample}.report.pdf'
pacbio_stats_file = processed_data_pacbio + '/{sample}.summary.tsv'
pacbio_run_summary = processed_data_pacbio + '/{sample}.examples.txt'
pacbio_rRNA_out = processed_data_pacbio + '/{sample}.rRNA.txt'
pacbio_overview = '/data/rajewsky/projects/slide_seq/.config/pacbio_overview.pdf'
pacbio_overview_csv = '/data/rajewsky/projects/slide_seq/.config/pacbio_overview.csv'
pacbio_bead_overview = '/data/rajewsky/projects/slide_seq/.config/pacbio_bead_overview.pdf'

include: 'snakemake/pacbio.smk' 

#############
# Main rule #
#############
rule all:
    input:
        #get_final_output_files(fastqc_pattern, ext = fastqc_ext, mate = [1,2]),
        # this will also create the clean dge
        get_united_output_files(automated_report),
        get_united_output_files(united_qc_sheet)

########################
# CREATE METADATA FILE #
########################
rule create_project_df_file:
    output:
        project_df_file
    run:
        project_df.to_csv(output[0], index=False)
        os.system('chmod 664 %s' % (output[0]))

rule create_sample_overview:
    input:
        project_df_file
    output:
        sample_overview_file
    script:
        'create_sample_overview.Rmd'

rule create_sample_db:
    input:
        project_df_file
    output:
        sample_read_metrics_db
    script:
        'snakemake/scripts/create_sample_db.R'

################
# DOWNSAMPLING #
################
include: 'snakemake/downsample.smk'

rule downsample:
    input:
        get_united_output_files(downsample_saturation_analysis, projects = config['downsample']['projects'], samples = config['downsample']['samples'])

#################
# MERGE SAMPLES #
#################
include: 'snakemake/merge_samples.smk'

#########
# RULES #
#########
ruleorder: link_raw_reads > link_demultiplexed_reads 

rule demultiplex_data:
    params:
        demux_barcode_mismatch=lambda wildcards: int(get_metadata('demux_barcode_mismatch', demux_dir = wildcards.demux_dir)),
        sample_sheet=lambda wildcards: get_metadata('sample_sheet', demux_dir = wildcards.demux_dir),
        flowcell_id=lambda wildcards:  get_metadata('flowcell_id', demux_dir = wildcards.demux_dir),
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
            -p {threads}            

            echo "demux finished: $(date)" > {output}
        """

rule link_demultiplexed_reads:
    input:
        ancient(unpack(get_demux_indicator))
    output:
        raw_reads_pattern
    params:
        demux_dir = lambda wildcards: expand(demux_dir_pattern,
            demux_dir = get_metadata('demux_dir', sample_id = wildcards.sample,
                                     project_id = wildcards.project)),
        reads_folder = raw_data_illumina_reads
    shell:
        """
        mkdir -p {params.reads_folder}

        find {params.demux_dir} -type f -wholename '*/{wildcards.sample}/*R{wildcards.mate}*.fastq.gz' -exec ln -sr {{}} {output} \; 
        """

def get_reads(wildcards):
    return([get_metadata('R'+ wildcards.mate, sample_id = wildcards.sample, project_id = wildcards.project)])

rule link_raw_reads:
    input:
        unpack(get_reads)
    output:
        raw_reads_pattern
    shell:
        """
        ln -s {input} {output}
        """

rule zcat_pipe:
    input: "{name}.fastq.gz"
    output: pipe("{name}.fastq")
    shell: "zcat {input} >> {output}"

dropseq_tagged_pipe = dropseq_tagged.replace('.bam', '.uncompressed.bam')

rule reverse_first_mate:
    input:
        # these implicitly depend on the raw reads via zcat_pipes
        R1_unpacked = raw_reads_mate_1.replace('fastq.gz', 'fastq'),
        R2_unpacked = raw_reads_mate_2.replace('fastq.gz', 'fastq')
    params:
        bc = lambda wildcards: get_bc_preprocess_settings(wildcards)
    output:
        assigned = pipe(dropseq_tagged_pipe),
        unassigned = dropseq_unassigned,
        bc_stats = reverse_reads_mate_1.replace(reads_suffix, ".bc_stats.tsv")
    log:
        reverse_reads_mate_1.replace(reads_suffix, ".preprocessing.log")
    threads: 2
    shell:
        "python {repo_dir}/snakemake/scripts/preprocess_read1.py "
        "--sample={wildcards.sample} "
        "--read1={input.R1_unpacked} "
        "--read2={input.R2_unpacked} "
        "--parallel={threads} "
        "--save-stats={output.bc_stats} "
        "--log-file={log} "
        "--bc1-ref={params.bc.bc1_ref} "
        "--bc2-ref={params.bc.bc2_ref} "
        "--bc1-cache={params.bc.bc1_cache} "
        "--bc2-cache={params.bc.bc2_cache} "
        "--threshold={params.bc.score_threshold} "
        "--cell='{params.bc.cell}' "
        "--cell-raw='{params.bc.cell_raw}' "
        "--out-format=bam "
        "--out-unassigned={output.unassigned} "
        "--out-assigned={output.assigned} "
        "--UMI='{params.bc.UMI}' "
        "--bam-tags='{params.bc.bam_tags}' "

rule compress_dropseq_tagged:
    input: dropseq_tagged_pipe
    output: dropseq_tagged
    shell: "sambamba view -h -l9 -f bam {input} > {output}"

rule run_fastqc:
    input:
        # we need to use raw reads here, as later during "reversing" we do the umi
        # extraction, and barcode identification (for the combinatorial barcoding)
        # in order for R1 to have the same pattern
        raw_reads_pattern
    output:
        fastqc_pattern
    params:
        output_dir = fastqc_root 
    threads: 4
    shell:
        """
        mkdir -p {params.output_dir}

        {fastqc_command} -t {threads} -o {params.output_dir} {input}
        """

rule index_bam_file:
    input:
        dropseq_final_bam
    output:
        dropseq_final_bam_ix 
    shell:
       "sambamba index {input}"

rule get_barcode_readcounts:
    # this rule takes the final.bam file (output of the dropseq pipeline) and creates a barcode read count file
    input:
        united_final_bam
    output:
        united_barcode_readcounts
    params:
        cell_barcode_tag = lambda wildcards: get_bam_tag_names(
            project_id = wildcards.united_project,
            sample_id = wildcards.united_sample)['{cell}'],
    shell:
        """
        {dropseq_tools}/BamTagHistogram \
        I= {input} \
        O= {output}\
        TAG={params.cell_barcode_tag}
        """

rule create_top_barcodes:
    input:
        united_barcode_readcounts
    output:
        united_top_barcodes
    params:
        n_beads=lambda wildcards: get_downstream_analysis_variables(sample_id = wildcards.united_sample, project_id = wildcards.united_project)['expected_n_beads'],
    shell:
        "set +o pipefail; zcat {input} | cut -f2 | head -{params.n_beads} > {output}"

rule clean_top_barcodes:
    input:
        united_top_barcodes
    output:
        united_top_barcodes_clean
    script:
        'snakemake/scripts/clean_top_barcodes.py'

rule create_dge:
    # creates the dge. depending on if the dge has _cleaned in the end it will require the
    # topBarcodesClean.txt file or just the regular topBarcodes.txt
    input:
        unpack(get_top_barcodes),
        reads=united_final_bam
    output:
        dge=dge_out,
        dge_summary=dge_out_summary
    params:
        dge_root = dge_root,
        dge_extra_params = lambda wildcards: get_dge_extra_params(wildcards),
        cell_barcode_tag = lambda wildcards: get_bam_tag_names(
            project_id = wildcards.united_project,
            sample_id = wildcards.united_sample)['{cell}'],
        umi_tag = lambda wildcards: get_bam_tag_names(
            project_id = wildcards.united_project,
            sample_id = wildcards.united_sample)['{UMI}']
    shell:
        """
        mkdir -p {params.dge_root}

        {dropseq_tools}/DigitalExpression \
        -m 16g \
        I= {input.reads}\
        O= {output.dge} \
        SUMMARY= {output.dge_summary} \
        CELL_BC_FILE={input.top_barcodes} \
        CELL_BARCODE_TAG={params.cell_barcode_tag} \
        MOLECULAR_BARCODE_TAG={params.umi_tag} \
        TMP_DIR={temp_dir} \
        {params.dge_extra_params}
        """

rule create_qc_parameters:
    params:
        sample_params=lambda wildcards: get_qc_sheet_parameters(wildcards.project, wildcards.sample)
    output:
        qc_sheet_parameters_file
    script:
        "analysis/qc_sequencing_create_parameters_from_sample_sheet.py"

rule parse_ribo_log:
    input: united_ribo_depletion_log
    output: united_parsed_ribo_depletion_log
    script: 'snakemake/scripts/parse_ribo_log.py'

rule create_qc_sheet:
    input:
        unpack(get_dge_type),
        star_log = united_star_log,
        reads_type_out=united_reads_type_out,
        parameters_file=united_qc_sheet_parameters_file,
        read_counts = united_barcode_readcounts,
        strand_info = united_strand_info,
        ribo_log=united_parsed_ribo_depletion_log
    output:
        united_qc_sheet
    script:
        "analysis/qc_sequencing_create_sheet.Rmd"

rule run_automated_analysis:
    input:
        unpack(get_dge_type)
    output:
       **automated_result_files
    script:
        'analysis/automated_analysis.py'
        
rule create_automated_report:
    input:
        unpack(get_puck_file),
        star_log=united_star_log,
        parameters_file=united_qc_sheet_parameters_file,
        **automated_result_files
    output:
        automated_report
    script:
        'analysis/automated_analysis_create_report.Rmd'

rule split_final_bam:
    input:
        united_final_bam
    output:
        temp(united_split_reads_sam_files),
        united_split_reads_read_type,
        united_split_reads_strand_type
    params:
        prefix=united_split_reads_root
    shell:
        "sambamba view -F 'mapping_quality==255' -h {input} | python {repo_dir}/snakemake/scripts/split_reads_by_strand_info.py --prefix {params.prefix} /dev/stdin"

rule split_reads_sam_to_bam:
    input:
        united_split_reads_sam_pattern
    output:
        united_split_reads_bam_pattern
    threads: 2
    shell:
        "sambamba view -S -h -f bam -t {threads} -o {output} {input}"

rule get_unmapped_reads:
    input:
        united_final_bam
    output:
        united_unmapped_bam
    threads: 2
    shell:
        "sambamba view -F 'unmapped' -h -f bam -t {threads} -o {output} {input}"


rule map_to_rRNA:
    input:
        raw_reads_mate_2
    output:
        ribo_depletion_log
    params:
        index= lambda wildcards: get_rRNA_index(wildcards)['rRNA_index'] 
    run:
        if wildcards.sample in human_mouse_samples:
            shell("bowtie2 -x {params.index} -U {input} -p 20 --very-fast-local > /dev/null 2> {output}")
        else:
            shell("echo 'no_rRNA_index' > {output}")

rule calculate_kmer_counts:
    input:
        raw_reads_mate_1
    output:
        kmer_stats_file
    params:
        kmer_len = lambda wildcards: wildcards.kmer_len
    script:
        'snakemake/scripts/kmer_stats_from_fastq.py'
