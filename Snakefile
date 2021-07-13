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
config['puck_data']['root'] = config['microscopy_out']

# set root dir where the processed_data goes
project_dir = config['root_dir'] + '/projects/{project}'


# moved barcode_flavor assignment here so that additional samples/projects are equally processed
project_df = create_project_df(config)

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

illumina_root = config['root_dir'] + '/projects/{project}/processed_data/{sample}/illumina'
complete_data_root = illumina_root + '/complete_data'

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
fastqc_root = raw_data_illumina + '/fastqc'
fastqc_pattern = fastqc_root + '/{sample}_R{mate}_fastqc.{ext}'
fastqc_ext = ['zip', 'html']

########################
# UNIQUE PIPELINE VARS #
########################
# set the tool script directories
picard_tools = config['external_bin']['picard_tools']
dropseq_tools = config['external_bin']['dropseq_tools']

reports_dir = complete_data_root + '/reports'
tmp_dir = complete_data_root + '/tmp'
smart_adapter = config['adapters']['smart']

merged_reads = complete_data_root + '/unaligned.bam'

###
# splitting the reads
###

split_reads_root = complete_data_root + '/split_reads{polyA_trimmed}/'

split_reads_sam_names = ['plus_plus', 'plus_minus', 'minus_minus', 'minus_plus', 'plus_AMB', 'minus_AMB']
split_reads_sam_pattern = split_reads_root + '{file_name}.sam'
split_reads_bam_pattern = split_reads_root + '{file_name}.bam'

split_reads_sam_files = [split_reads_root + x + '.sam' for x in split_reads_sam_names]

split_reads_strand_type = split_reads_root + 'strand_type_num.txt'
split_reads_read_type = split_reads_root + 'read_type_num.txt'

#######################
# post dropseq and QC #
#######################

qc_sheet = complete_data_root +'/qc_sheet_{sample}_{puck}.html'
reads_type_out = split_reads_read_type
barcode_readcounts = complete_data_root + '/out_readcounts{polyA_trimmed}.txt.gz'
strand_info = split_reads_strand_type

# united final.bam
top_barcodes = complete_data_root + '/topBarcodes{polyA_trimmed}.{n_beads}_beads.txt'
top_barcodes_clean = complete_data_root + '/topBarcodesClean{polyA_trimmed}.{n_beads}_beads.txt'

# united dgu
dge_root = complete_data_root + '/dge'
dge_out_prefix = dge_root + '/dge{dge_type}{dge_cleaned}'
dge_out_suffix = '{polyA_trimmed}{mm_included}.{n_beads}_beads'
dge_out = dge_out_prefix + dge_out_suffix + '.txt.gz'
dge_out_summary = dge_out_prefix + dge_out_suffix + '.summary.txt'
dge_types = ['.exon', '.intron', '.all', '.Reads_exon', '.Reads_intron', '.Reads_all']

dge_all_summary = complete_data_root +  '/dge/dge.all'+ dge_out_suffix + '.summary.txt'
# cleaned dge and summary
dge_all_cleaned = complete_data_root +  '/dge/dge.all.cleaned.txt.gz'
dge_all_cleaned_summary = complete_data_root +  '/dge/dge.all.cleaned' +dge_out_suffix + '.summary.txt'

# dge and summary
dge_all = complete_data_root +  '/dge/dge_all' + dge_out_suffix + '.txt.gz'

# kmer stats per position
kmer_stats_file = complete_data_root + '/kmer_stats/{kmer_len}mer_counts.csv'

# map paired-end to check errors
paired_end_prefix = complete_data_root + '/mapped_paired_end/'
paired_end_sam = paired_end_prefix + '{sample}_paired_end.sam'
paired_end_bam = paired_end_prefix + '{sample}_paired_end.bam'
paired_end_flagstat = paired_end_prefix + '{sample}_paired_end_flagstat.txt'
paired_end_log = paired_end_prefix + '{sample}_paired_end.log'
paired_end_mapping_stats = paired_end_prefix + '{sample}_paired_end_mapping_stats.txt'

# automated analysis
automated_analysis_root = complete_data_root + '/automated_analysis/{run_mode}/umi_cutoff_{umi_cutoff}'
automated_report = automated_analysis_root + '/{sample}_{puck}_illumina_automated_report.html'

automated_analysis_result_file = automated_analysis_root + '/results.h5ad'

automated_analysis_processed_data_files = {
    'cluster_markers': '/top10_cluster_markers.csv',
    'obs_df': '/obs_df.tsv',
    'var_df': '/var_df.tsv',
    'long_expr_df': '/long_expr_df.tsv'
    }

# prepend automated_result_root
automated_analysis_processed_data_files = {key: automated_analysis_root + value for key, value in automated_analysis_processed_data_files.items()}

# blast out
blast_db_primers = repo_dir + '/sequences/primers.fa'
blast_db_primers_files = [blast_db_primers + '.' + x for x in ['nhr', 'nin', 'nog', 'nsd', 'nsi', 'nsq']]
blast_header_out = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore"
barcode_blast_out = complete_data_root + '/cell_barcode_primer_blast_out.txt'

# downsample vars
downsample_root = illumina_root + '/downsampled_data'

# in silico repo depletion
ribo_depletion_log = complete_data_root + '/ribo_depletion_log.txt'
ribo_depletion_log = complete_data_root + '/ribo_depletion_log.txt'
parsed_ribo_depletion_log = complete_data_root + '/parsed_ribo_depletion_log.txt'

# global wildcard constraints
wildcard_constraints:
    sample='(?!merged_).+',
    project='(?!merged_).+',
    dge_cleaned='|.cleaned',
    dge_type = '|'.join(dge_types),
    pacbio_ext = 'fq|fastq',
    polyA_trimmed = '|.polyA_trimmed',
    mm_included = '|.mm_included',
    n_beads = '[0-9]+'

# #########################
#  dropseq rules and vars #
# #########################
tagged_bam = complete_data_root + '/unaligned_bc_tagged.bam'
unassigned = complete_data_root + '/unaligned_bc_unassigned.bam'

# trim smart adapter from the reads
tagged_trimmed_bam = complete_data_root + '/unaligned_bc_tagged_trimmed.bam'

# trim polyA overheang if exists
tagged_polyA_trimmed_bam = complete_data_root + '/unaligned_bc_tagged.polyA_trimmed.bam'

tagged_bam_pattern = complete_data_root + '/unaligned_bc_tagged{polyA_trimmed}.bam'

# mapped reads
#mapped_reads_sorted_headerless = complete_data_root + '/star{polyA_trimmed}.Aligned.headerless.out.bam'
mapped_reads_qname_sorted = complete_data_root + '/star{polyA_trimmed}.Aligned.out.bam'
star_log_file = complete_data_root + '/star{polyA_trimmed}.Log.final.out'
star_prefix  = complete_data_root + '/star{polyA_trimmed}.'

# final bam file
final_bam = complete_data_root + '/final{polyA_trimmed}.bam'
final_bam_mm_included = complete_data_root + '/final{dge_type}{dge_cleaned}{polyA_trimmed}.mm_included.bam'
final_bam_pattern = complete_data_root + '/final{polyA_trimmed}{mm_included}.bam'

# include dropseq
include: 'snakemake/dropseq.smk'

################################
# Final output file generation #
################################
def get_output_files(pattern, projects = [], samples = [],
                            skip_projects = [], skip_samples = [], **kwargs):
    out_files = []
    df = project_df

    if samples or projects:
        # one of the lists is not empty
        df = df.query('sample_id in @samples or project_id in @projects')

    if skip_samples:
        df = df.query('sample_id not in @skip_samples')

    if skip_projects is not None:
        df = df.query('project_id not in @skip_projects')

    for index, row in df.iterrows():
        for run_mode in row['run_mode']:
            run_mode_variables = get_run_mode_variables(run_mode)
            out_files = out_files + expand(pattern,
                project = index[0],
                sample = index[1],
                puck=row['puck_id'], 
                run_mode=run_mode,
                umi_cutoff=run_mode_variables['umi_cutoff'],
                **kwargs)

    return out_files

human_mouse_samples = project_df[project_df.species.isin(['human', 'mouse'])].index.get_level_values('sample_id')

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
        #get_final_output_files(fastqc_pattern, skip_merged = True, ext = fastqc_ext, mate = [1,2]),
        # this will also create the clean dge
        #get_output_files(automated_report),
        get_output_files(qc_sheet)

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
        'snakemake/scripts/create_sample_overview.Rmd'

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
#include: 'snakemake/downsample.smk'

#rule downsample:
#    input:
 #       get_output_files(downsample_saturation_analysis, projects = config['downsample']['projects'], samples = config['downsample']['samples'])

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
        output_dir= lambda wildcards: expand(demux_dir_pattern, demux_dir=wildcards.demux_dir)
    input:
        lambda wildcards: get_metadata('basecalls_dir', demux_dir = wildcards.demux_dir)
    output:
        demux_indicator
    threads: 8
    shell:
        """
        bcl2fastq \
            --no-lane-splitting --fastq-compression-level=9 \
            --mask-short-adapter-reads 15 \
            --barcode-mismatch {params.demux_barcode_mismatch} \
            --output-dir {params.output_dir} \
            --sample-sheet {params.sample_sheet} \
            --runfolder-dir {input} \
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
    ###
    # R1 and R2 for demultiplexed reads will return none
    ### 
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
    output: temp("{name}.fastq")
    threads: 2
    shell: "unpigz --keep --processes {threads} --stdout $(readlink {input}) >> {output}"

tagged_pipe = tagged_bam.replace('.bam', '.uncompressed.bam')

rule reverse_first_mate:
    input:
        # these implicitly depend on the raw reads via zcat_pipes
        R1 = raw_reads_mate_1,
        R2 = raw_reads_mate_2
    params:
        bc = lambda wildcards: get_bc_preprocess_settings(wildcards)
    output:
        assigned = temp(tagged_bam),
        unassigned = unassigned,
        bc_stats = reverse_reads_mate_1.replace(reads_suffix, ".bc_stats.tsv")
    log:
        reverse_reads_mate_1.replace(reads_suffix, ".preprocessing.log")
    threads: 16
    shell:
        "python {repo_dir}/snakemake/scripts/preprocess_read1.py "
        "--sample={wildcards.sample} "
        "--read1={input.R1} "
        "--read2={input.R2} "
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
        "--out-assigned=/dev/stdout "
        "--UMI='{params.bc.UMI}' "
        "--bam-tags='{params.bc.bam_tags}' "
        "--min-opseq-score={params.bc.min_opseq_score} "
        "| samtools view -bh /dev/stdin > {output.assigned} "

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

        fastqc -t {threads} -o {params.output_dir} {input}
        """

rule get_barcode_readcounts:
    # this rule takes the final.bam file (output of the dropseq pipeline) and creates a barcode read count file
    input:
        final_bam
    output:
        barcode_readcounts
    params:
        cell_barcode_tag = lambda wildcards: get_bam_tag_names(
            project_id = wildcards.project,
            sample_id = wildcards.sample)['{cell}'],
    shell:
        """
        {dropseq_tools}/BamTagHistogram \
        I= {input} \
        O= {output}\
        TAG={params.cell_barcode_tag}
        """

rule create_top_barcodes:
    input:
        barcode_readcounts
    output:
        top_barcodes
    shell:
        "set +o pipefail; zcat {input} | cut -f2 | head -{wildcards.n_beads} > {output}"

rule clean_top_barcodes:
    input:
        top_barcodes
    output:
        top_barcodes_clean
    script:
        'snakemake/scripts/clean_top_barcodes.py'

rule create_dge:
    # creates the dge. depending on if the dge has _cleaned in the end it will require the
    # topBarcodesClean.txt file or just the regular topBarcodes.txt
    input:
        unpack(get_top_barcodes),
        unpack(get_mapped_final_bam)
    output:
        dge=dge_out,
        dge_summary=dge_out_summary
    params:
        dge_root = dge_root,
        dge_extra_params = lambda wildcards: get_dge_extra_params(wildcards),
        cell_barcode_tag = lambda wildcards: get_bam_tag_names(
            project_id = wildcards.project,
            sample_id = wildcards.sample)['{cell}'],
        umi_tag = lambda wildcards: get_bam_tag_names(
            project_id = wildcards.project,
            sample_id = wildcards.sample)['{UMI}']
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

rule parse_ribo_log:
    input: ribo_depletion_log
    output: parsed_ribo_depletion_log
    script: 'snakemake/scripts/parse_ribo_log.py'

rule create_qc_sheet:
    input:
        unpack(get_dge_type),
        unpack(get_qc_sheet_input_files),
        ribo_log=parsed_ribo_depletion_log
    params:
        run_modes = lambda wildcards: get_run_modes_from_sample(wildcards.project,
            wildcards.sample)
    output:
        qc_sheet
    script:
        "analysis/qc_sequencing_create_sheet.Rmd"

rule run_automated_analysis:
    input:
        unpack(get_puck_file),
        unpack(get_dge_type)
    output:
        automated_analysis_result_file
    params:
        downstream_variables = lambda wildcards: get_run_mode_variables(wildcards.run_mode)
    threads: 2
    script:
        'analysis/automated_analysis.py'

rule create_automated_analysis_processed_data_files:
    input:
        automated_analysis_result_file
    output:
        **automated_analysis_processed_data_files
    script:
        'analysis/automated_analysis_create_processed_data_files.py'
        
rule create_automated_report:
    input:
        #star_log=star_log_file,
        **automated_analysis_processed_data_files
    output:
        automated_report
    params:
        downstream_variables = lambda wildcards: get_run_mode_variables(wildcards.run_mode),
        r_shared_scripts= repo_dir + '/analysis/shared_functions.R'
    script:
        'analysis/automated_analysis_create_report.Rmd'

rule split_final_bam:
    input:
        final_bam
    output:
        temp(split_reads_sam_files),
        split_reads_read_type,
        split_reads_strand_type
    params:
        prefix=split_reads_root
    shell:
        """
        sambamba view -F 'mapping_quality==255' -h {input} | \
        python {repo_dir}/snakemake/scripts/split_reads_by_strand_info.py \
        --prefix {params.prefix} /dev/stdin
        """

rule split_reads_sam_to_bam:
    input:
        split_reads_sam_pattern
    output:
        split_reads_bam_pattern
    threads: 2
    shell:
        "sambamba view -S -h -f bam -t {threads} -o {output} {input}"

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
