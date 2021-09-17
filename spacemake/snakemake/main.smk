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
import scanpy as sc

from spacemake.util import dge_to_sparse_adata, attach_barcode_file,\
    create_meshed_adata, parse_barcode_file
from spacemake.cmdline_utils import ProjectDF, ConfigFile

################
# Shell prefix #
################
shell.prefix('set +o pipefail; JAVA_TOOL_OPTIONS="-Xmx8g -Xss2560k" ; ')

####
# this file should contain all sample information, sample name etc.
####
# configfile should be loaded from command line

###############
# Global vars #
###############
temp_dir = config['temp_dir']
repo_dir = os.path.dirname(workflow.snakefile)
spacemake_dir = os.path.dirname(os.path.dirname(workflow.snakefile))
# create puck_data root directory from pattern
config['puck_data']['root'] = config['microscopy_out']

# set root dir where the processed_data goes
project_dir = os.path.join(config['root_dir'], 'projects/{project}')

# moved barcode_flavor assignment here so that additional samples/projects are equally processed
project_df = ProjectDF(config['project_df'], ConfigFile('config.yaml'))

#################
# DIRECTORY STR #
#################
raw_data_root = project_dir + '/raw_data'
raw_data_illumina = raw_data_root + '/illumina'
raw_data_illumina_reads = raw_data_illumina + '/reads/raw'
raw_data_illumina_reads_reversed = raw_data_illumina + '/reads/reversed'
processed_data_root = project_dir + '/processed_data/{sample}'
processed_data_illumina = processed_data_root + '/illumina'

reports_root = os.path.join(config['root_dir'], 'reports')

project_df_file = reports_root + '/project_df.csv'
sample_overview_file = reports_root + '/sample_overview.html'
sample_read_metrics_db = reports_root + '/sample_read_metrics_db.tsv'

illumina_root = project_dir + '/processed_data/{sample}/illumina'
complete_data_root = illumina_root + '/complete_data'

##############
# Demux vars #
##############
# Undetermined files pattern
# they are the output of bcl2fastq, and serve as an indicator to see if the demultiplexing has finished
demux_dir_pattern = os.path.join(config['root_dir'], 'raw_data/demultiplex_data/{demux_dir}')
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

###
# splitting the reads
###

split_reads_root = complete_data_root + '/split_reads{polyA_adapter_trimmed}/'

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
barcode_readcounts_suffix = '{polyA_adapter_trimmed}.txt.gz'
barcode_readcounts = complete_data_root + '/out_readcounts' + barcode_readcounts_suffix
strand_info = split_reads_strand_type

# united final.bam
top_barcodes_suffix = '{polyA_adapter_trimmed}.{n_beads}_beads.txt'
top_barcodes = complete_data_root + '/topBarcodes' + top_barcodes_suffix
top_barcodes_clean = complete_data_root + '/topBarcodesClean' + top_barcodes_suffix
spatial_barcodes = complete_data_root + '/spatialBarcodes.txt'

# united dgu
dge_root = complete_data_root + '/dge'
dge_out_prefix = dge_root + '/dge'
dge_out_suffix = '{dge_type}{dge_cleaned}{polyA_adapter_trimmed}{mm_included}'
dge_out = dge_out_prefix + dge_out_suffix + '.{n_beads}_beads.txt.gz'
dge_out_summary = dge_out_prefix + dge_out_suffix + '.{n_beads}_beads.summary.txt'
dge_out_h5ad = dge_out_prefix + dge_out_suffix + '.{n_beads}_beads.h5ad'
dge_out_h5ad_obs = dge_out_prefix + dge_out_suffix + '.{n_beads}_beads.obs.csv'
dge_spatial = dge_out_prefix + dge_out_suffix + '.spatial_beads.h5ad'
dge_spatial_obs = dge_out_prefix + dge_out_suffix + '.spatial_beads.obs.csv'
dge_spatial_mesh_suffix = '.spatial_beads.mesh_{spot_diameter_um}_{spot_distance_um}'
dge_spatial_mesh_prefix = dge_out_prefix + dge_out_suffix + dge_spatial_mesh_suffix
dge_spatial_mesh = dge_spatial_mesh_prefix + '.h5ad'
dge_spatial_mesh_obs = dge_spatial_mesh_prefix + '.obs.csv'

dge_types = ['.exon', '.intron', '.all', '.Reads_exon', '.Reads_intron', '.Reads_all']

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
    'var_df': '/var_df.tsv'
    }

# prepend automated_result_root
automated_analysis_processed_data_files = {key: automated_analysis_root + value for key, value in automated_analysis_processed_data_files.items()}

# in silico repo depletion
ribo_depletion_log = complete_data_root + '/ribo_depletion_log.txt'
parsed_ribo_depletion_log = complete_data_root + '/parsed_ribo_depletion_log.txt'

# #########################
#  dropseq rules and vars #
# #########################
tagged_bam = complete_data_root + '/unaligned_bc_tagged.bam'
unassigned = complete_data_root + '/unaligned_bc_unassigned.bam'

# trim smart adapter from the reads
tagged_trimmed_bam = complete_data_root + '/unaligned_bc_tagged_trimmed.bam'

# trim polyA overheang if exists
tagged_polyA_adapter_trimmed_bam = complete_data_root + '/unaligned_bc_tagged.polyA_adapter_trimmed.bam'

tagged_bam_pattern = complete_data_root + '/unaligned_bc_tagged{polyA_adapter_trimmed}.bam'

# mapped reads
star_prefix  = complete_data_root + '/star{polyA_adapter_trimmed}.'
star_log_file = star_prefix + 'Log.final.out'

# final bam file
final_bam_suffix = '/final{polyA_adapter_trimmed}'
final_bam = complete_data_root + final_bam_suffix + '.bam'
bam_mm_included_pipe_suffix = '{dge_type}{dge_cleaned}{polyA_adapter_trimmed}.mm_included.bam'
final_bam_mm_included_pipe = complete_data_root + '/final' + bam_mm_included_pipe_suffix

# index settings
star_index = 'species_data/{species}/star_index'
bt2_rRNA_index_dir = 'species_data/{species}/bt2_rRNA_index'
bt2_rRNA_index_basename = bt2_rRNA_index_dir + '/{species}_rRNA'

################################
# Final output file generation #
################################
def get_output_files(pattern, projects = [], samples = [],
                            skip_projects = [], skip_samples = [], **kwargs):
    out_files = []
    df = project_df.df

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


####################
# HELPER FUNCTIONS #
####################
include: 'scripts/snakemake_helper_functions.py'

######################### 
# INCLUDE OTHER MODULES #
#########################
include: 'downsample.smk'
include: 'dropseq.smk'
include: 'pacbio.smk' 


# global wildcard constraints
wildcard_constraints:
    sample='(?!merged_).+',
    project='(?!merged_).+',
    dge_cleaned='|.cleaned',
    dge_type = '|'.join(dge_types),
    pacbio_ext = 'fq|fastq',
    polyA_adapter_trimmed = '|.polyA_adapter_trimmed',
    mm_included = '|.mm_included',
    n_beads = '[0-9]+|spatial',
    spot_diameter_um = '[0-9]+',
    spot_distance_um = '[0-9]+'

#############
# Main rule #
#############
rule all:
    input:
        #get_final_output_files(fastqc_pattern, skip_merged = True, ext = fastqc_ext, mate = [1,2]),
        # this will also create the clean dge
        get_output_files(automated_report),
        get_output_files(qc_sheet)

#####################
# DOWNSAMPLE MODULE # 
#####################

#rule downsample:
#    input:
 #       get_output_files(downsample_saturation_analysis, projects = config['downsample']['projects'], samples = config['downsample']['samples'])

#################
# MERGE SAMPLES #
#################
include: 'merge_samples.smk'

#########
# RULES #
#########
ruleorder: link_raw_reads > link_demultiplexed_reads 

rule demultiplex_data:
    params:
        demux_barcode_mismatch=lambda wildcards: int(project_df.get_metadata('demux_barcode_mismatch', demux_dir = wildcards.demux_dir)),
        sample_sheet=lambda wildcards: project_df.get_metadata('sample_sheet', demux_dir = wildcards.demux_dir),
        output_dir= lambda wildcards: expand(demux_dir_pattern, demux_dir=wildcards.demux_dir)
    input:
        lambda wildcards: project_df.get_metadata('basecalls_dir', demux_dir = wildcards.demux_dir)
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
            demux_dir = project_df.get_metadata('demux_dir', sample_id = wildcards.sample,
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
    return([project_df.get_metadata('R'+ wildcards.mate, sample_id = wildcards.sample, project_id = wildcards.project)])

rule link_raw_reads:
    input:
        unpack(get_reads)
    output:
        raw_reads_pattern
    shell:
        """
        ln -rs {input} {output}
        """

rule zcat_pipe:
    input: "{name}.fastq.gz"
    output: temp("{name}.fastq")
    threads: 2
    shell: "unpigz --keep --processes {threads} --stdout $(readlink {input}) >> {output}"

rule reverse_first_mate:
    input:
        # these implicitly depend on the raw reads via zcat_pipes
        R1 = raw_reads_mate_1,
        R2 = raw_reads_mate_2
    params:
        bc = lambda wildcards: get_bc_preprocess_settings(wildcards)
    output:
        assigned = tagged_bam,
        unassigned = unassigned,
        bc_stats = reverse_reads_mate_1.replace(reads_suffix, ".bc_stats.tsv")
    log:
        reverse_reads_mate_1.replace(reads_suffix, ".preprocessing.log")
    threads: 4
    shell:
        "python {spacemake_dir}/preprocess.py "
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
        unpack(get_final_bam)
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
        "zcat {input} | cut -f2 | head -{wildcards.n_beads} > {output}"

rule clean_top_barcodes:
    input:
        top_barcodes
    output:
        top_barcodes_clean
    script:
        'scripts/clean_top_barcodes.py'

rule create_spatial_barcodes:
    input:
        unpack(get_puck_file)
    output:
        temp(spatial_barcodes)
    run:
        bc = parse_barcode_file(input[0])
        bc['cell_bc'] = bc.index
        bc[['cell_bc']].to_csv(output[0], header=False, index=False)

rule create_dge:
    # creates the dge. depending on if the dge has _cleaned in the end it will require the
    # topBarcodesClean.txt file or just the regular topBarcodes.txt
    input:
        unpack(get_top_barcodes),
        unpack(get_dge_input_bam)
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

rule create_h5ad_dge:
    input:
        unpack(get_puck_file),
        dge_txt = dge_out,
        dge_summary_txt = dge_out_summary
    # output here will either be n_beads=number, n_beads=spatial
    output: dge_out_h5ad, dge_out_h5ad_obs
    run:
        adata = dge_to_sparse_adata(
            input['dge_txt'],
            input['dge_summary_txt'])
        if 'barcode_file' in input.keys() and wildcards.n_beads == 'spatial':
            adata = attach_barcode_file(adata, input['barcode_file'])
        adata.write(output[0])
        adata.obs.to_csv(output[1])

rule create_mesh_spatial_dge:
    input:
        dge_spatial
    output:
        dge_spatial_mesh,
        dge_spatial_mesh_obs
    params:
        puck_data = lambda wildcards: project_df.get_puck(
            project_id = wildcards.project,
            sample_id = wildcards.sample)
    run:
        adata = sc.read(input[0])
        adata = create_meshed_adata(adata,
            width_um = float(params['puck_data']['width_um']),
            bead_diameter_um = float(params['puck_data']['spot_diameter_um']),
            spot_diameter_um = float(wildcards.spot_diameter_um),
            spot_distance_um = float(wildcards.spot_distance_um))
        adata.write(output[0])
        adata.obs.to_csv(output[1])

rule parse_ribo_log:
    input: unpack(get_ribo_depletion_log)
    output: parsed_ribo_depletion_log
    script: 'scripts/parse_ribo_log.py'

rule create_qc_sheet:
    input:
        unpack(get_qc_sheet_input_files),
        ribo_log=parsed_ribo_depletion_log
    params:
        sample_info = lambda wildcards: project_df.get_sample_info(
            wildcards.project, wildcards.sample),
        run_modes = lambda wildcards: get_run_modes_from_sample(
            wildcards.project, wildcards.sample)
    output:
        qc_sheet
    script:
        "scripts/qc_sequencing_create_sheet.Rmd"

rule run_automated_analysis:
    input:
        unpack(get_automated_analysis_dge_input)
    output:
        automated_analysis_result_file
    params:
        run_mode_variables = lambda wildcards: get_run_mode_variables(wildcards.run_mode)
    script:
        'scripts/automated_analysis.py'

rule create_automated_analysis_processed_data_files:
    input:
        automated_analysis_result_file
    output:
        **automated_analysis_processed_data_files
    script:
        'scripts/automated_analysis_create_processed_data_files.py'
        
rule create_automated_report:
    input:
        #star_log=star_log_file,
        **automated_analysis_processed_data_files
    output:
        automated_report
    params:
        run_mode_variables = lambda wildcards: get_run_mode_variables(wildcards.run_mode),
        r_shared_scripts= repo_dir + '/scripts/shared_functions.R'
    script:
        'scripts/automated_analysis_create_report.Rmd'

rule split_final_bam:
    input:
        unpack(get_final_bam)
    output:
        temp(split_reads_sam_files),
        split_reads_read_type,
        split_reads_strand_type
    params:
        prefix=split_reads_root
    shell:
        """
        sambamba view -F 'mapping_quality==255' -h {input} | \
        python {repo_dir}/scripts/split_reads_by_strand_info.py \
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

rule create_rRNA_index:
    input:
        unpack(get_rRNA_genome)
    output:
        directory(bt2_rRNA_index_dir)
    params:
        basename=bt2_rRNA_index_basename
    shell:
        """
        mkdir -p {output}

        bowtie2-build --ftabchars 12 \
                      --offrate 1 \
                      {input} \
                      {params.basename}
        """

rule map_to_rRNA:
    input:
        unpack(get_bt2_rRNA_index),
        reads=raw_reads_mate_2
    output:
        ribo_depletion_log
    params:
        species=lambda wildcards: project_df.get_metadata(
            'species', project_id = wildcards.project,
            sample_id = wildcards.sample)
    run:
        if 'index' in input.keys():
            shell("bowtie2 -x {input.index}/{params.species}_rRNA -U {input.reads} -p 20 --very-fast-local > /dev/null 2> {output}")
        else:
            shell("echo 'no_rRNA_index' > {output}")
