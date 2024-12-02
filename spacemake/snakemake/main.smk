#########
# about #
#########
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

from spacemake.preprocess.dge import dge_to_sparse_adata, attach_barcode_file,\
    parse_barcode_file, load_external_dge, attach_puck
from spacemake.spatial.util import create_meshed_adata
import spacemake.spatial.puck_collection as puck_collection
from spacemake.project_df import ProjectDF
from spacemake.config import ConfigFile
from spacemake.errors import SpacemakeError

# load root dir_variable
project_root = config['root_dir']

################
# Shell prefix #
################
shell.prefix('set +o pipefail; JAVA_TOOL_OPTIONS="-Xmx8g -Xss2560k" ; ')

####
# this file should contain all sample information, sample name etc.
####
# configfile should be loaded from command line
# populate if not exists
config['samples'] = config.get('samples', [])
config['projects'] = config.get('projects', [])

###############
# Global vars #
###############
global_tmp = config['temp_dir']
repo_dir = os.path.dirname(workflow.snakefile)
spacemake_dir = os.path.dirname(os.path.dirname(workflow.snakefile))

import logging
smk_logger = logging.getLogger("spacemake.main.smk")

#######################
# DIRECTORY STRUCTURE #
#######################
include: 'variables.py'

########################
# UNIQUE PIPELINE VARS #
########################
# set the tool script directories
dropseq_tools = config['external_bin']['dropseq_tools']

reports_dir = complete_data_root + '/reports'
smart_adapter = config['adapters']['smart']


# moved barcode_flavor assignment here so that additional samples/projects are equally processed
project_df = ProjectDF(config['project_df'], ConfigFile.from_yaml('config.yaml'))

species_genome = 'species_data/{species}/genome.fa'
species_annotation = 'species_data/{species}/annotation.gtf'

####################
# HELPER FUNCTIONS #
####################
include: 'scripts/snakemake_helper_functions.py'

_module_output_hooks = []
def register_module_output_hook(hook, module="built-in"):
    _module_output_hooks.append( (hook, module) )

def get_module_outputs():
    outputs = []
    for hook, module in _module_output_hooks:
        for out in hook(project_df=project_df, config=config):
            smk_logger.debug(f"output provided by '{module}' module (via '{hook.__name__}'): '{out}'")
            outputs.append(out)
    
    return outputs

######################
# CHECKPOINT PARSERS #
######################
def get_pbf_from_checkpoint_pucks(wildcards):
    checkpoint_output = checkpoints.checkpoint_pucks.get(**wildcards).output[0]
    puck_barcode_file_ids = glob_wildcards(os.path.join(checkpoint_output, "{p}.chk")).p
    puck_barcode_file_ids = list(set(puck_barcode_file_ids))
    return puck_barcode_file_ids


def checkpoint_puck_collection(wildcards):
    puck_barcode_file_ids = get_pbf_from_checkpoint_pucks(wildcards)
    out_files = get_puck_collection_stitching_input(wildcards, puck_barcode_file_ids, to_mesh=False)
    return out_files


def checkpoint_puck_collection_mesh(wildcards):
    puck_barcode_file_ids = get_pbf_from_checkpoint_pucks(wildcards)
    out_files = get_puck_collection_stitching_input(wildcards, puck_barcode_file_ids, to_mesh=True)
    return out_files


def checkpoint_puck_files(wildcards):
    puck_barcode_file_ids = get_pbf_from_checkpoint_pucks(wildcards)
    out_files = {"dge": get_all_dges(wildcards, puck_barcode_file_ids),
                 "automated_report": get_output_files(automated_report,
                  data_root_type = 'complete_data', downsampling_percentage='', 
                  puck_barcode_file_ids=puck_barcode_file_ids, check_puck_collection=True),
                 "qc_report": get_output_files(qc_sheet,
                  data_root_type = 'complete_data', downsampling_percentage='', 
                   puck_barcode_file_ids=puck_barcode_file_ids, check_puck_collection=True),}
    return out_files


def checkpoint_barcode_files(wildcards):
    checkpoint_output = checkpoints.checkpoint_barcodes.get(**wildcards).output[0]
    puck_barcode_file_ids = glob_wildcards(os.path.join(checkpoint_output, "{p}.chk")).p
    puck_barcode_file_ids = list(set(puck_barcode_file_ids))
    out_files = get_barcode_files_matching_summary_input(wildcards, puck_barcode_file_ids)
    return out_files


######################### 
# INCLUDE OTHER MODULES #
#########################
include: 'downsample.smk'
include: 'mapping.smk'
include: 'dropseq.smk'
include: 'longread.smk'

#####################
# CUSTOM USER RULES #
#####################
# top-level targets to be injected as input-dependencies of the 
# master rule should be registered by calling register_module_output_hook()
# with a function that returns a list of output files.

if "custom_rules" in config:
    custom = os.path.join(config["pwd"], config["custom_rules"])
    include: custom


# global wildcard constraints
wildcard_constraints:
    umi_cutoff = r'\d+',
    dge_cleaned=r'|\.cleaned',
    dge_type = r'|'.join(dge_types),
    pacbio_ext = r'fq|fastq|bam',
    polyA_adapter_trimmed = r'|\.polyA_adapter_trimmed',
    mm_included = r'|\.mm_included',
    n_beads = r'[0-9]+|spatial|external',
    is_external = r'|\.external',
    spot_diameter_um = r'[0-9]+',
    spot_distance_um = r'[0-9]+|hexagon',
    data_root_type = r'complete_data|downsampled_data',
    downsampling_percentage = r'\/[0-9]+|',
    puck_barcode_file_id = r'(?!puck_collection)[^.]+',
    puck_barcode_file_id_qc = r'[^.]+'

#############
# Main rule #
#############
rule run_analysis:
    input:
        # get outputs from registered hooks
        get_module_outputs(),

        # get FastQC reports
        (
            get_output_files(
                    fastqc_pattern, ext = fastqc_ext, mate=['1', '2'],
                    data_root_type = 'complete_data', downsampling_percentage = '',
                    filter_merged=True) 
                if config['with_fastqc'] else []
        ),
        # get flag for DGE (based on checkpoint, i.e., not explicitly generating files)
        get_expanded_pattern_project_sample(dge_out_done),


##############
# DOWNSAMPLE #
##############
rule downsample:
    input:
        get_output_files(downsample_saturation_analysis,
            samples = config['samples'],
            projects = config['projects'],
            puck_barcode_file_matching_type = "spatial_matching")

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
        retVal=$(bcl2fastq \
            --no-lane-splitting --fastq-compression-level=9 \
            --mask-short-adapter-reads 15 \
            --barcode-mismatch {params.demux_barcode_mismatch} \
            --output-dir {params.output_dir} \
            --sample-sheet {params.sample_sheet} \
            --runfolder-dir {input} \
            -p {threads})
        
        if [ $retVAL==0 ]; then
            echo "demux finished: $(date)" > {output}
        fi
        exit $retVal
        """

rule link_demultiplexed_reads:
    input:
        ancient(unpack(get_demux_indicator))
    output:
        raw_reads_pattern
    params:
        demux_dir = lambda wildcards: expand(demux_dir_pattern,
            demux_dir = project_df.get_metadata('demux_dir', sample_id = wildcards.sample_id,
                                     project_id = wildcards.project_id)),
        reads_folder = raw_data_illumina_reads
    shell:
        """
        mkdir -p {params.reads_folder}

        find {params.demux_dir} -type f -wholename '*/{wildcards.sample_id}/*R{wildcards.mate}*.fastq.gz' -exec ln -sr {{}} {output} \; 
        """

rule link_raw_reads:
    input:
        unpack(get_reads)
    output:
        raw_reads_pattern
    run:
        if len(input) == 1:
            # either link raw reads
            shell("ln -rs {input} {output}")
        else:
            # or append reads together
            shell("cat {input} > {output}")
            

rule zcat_pipe:
    input: "{name}.fastq.gz"
    output: temp("{name}.fastq")
    threads: 2
    shell: "unpigz --keep --processes {threads} --stdout $(readlink {input}) >> {output}"

rule tag_reads_bc_umi:
    input:
        # these implicitly depend on the raw reads via zcat_pipes
        R1 = raw_reads_mate_1,
        R2 = raw_reads_mate_2
    params:
        bc = lambda wildcards: get_bc_preprocess_settings(wildcards)
    output:
        assigned = tagged_bam,
        log = tagged_bam_log
    log:
        reverse_reads_mate_1.replace(reads_suffix, ".preprocessing.log")
    threads: max(min(workflow.cores * 0.5, 16), 1)
    shell:
        "python {spacemake_dir}/bin/fastq_to_uBAM.py "
        "--sample={wildcards.sample_id} "
        "--read1={input.R1} "
        "--read2={input.R2} "
        "--parallel={threads} "
	    "--out-bam={output.assigned} "
        "--cell='{params.bc.cell}' "
        "--UMI='{params.bc.UMI}' "
        "--bam-tags='{params.bc.bam_tags}' "
        "--log-file='{output.log}' "

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
        barcode_readcounts,
        barcode_readcounts_log
    params:
        cell_barcode_tag = lambda wildcards: get_bam_tag_names(
            project_id = wildcards.project_id,
            sample_id = wildcards.sample_id)['{cell}']
    threads: max(min(workflow.cores * 0.5, 16), 1)
    shell:
        # {dropseq_tools}/BamTagHistogram -m 32g 
        """
        python {spacemake_dir}/bin/BamTagHistogram.py \
        --input {input} \
        --output {output[0]} \
        --tag {params.cell_barcode_tag} \
        --min-count 1 \
        --log-file {output[1]} \
        """
        #READ_MQ=0

rule merge_stats_prealigned_spatial_barcodes:
    input:
        unpack(get_barcode_files),
        bc_counts = barcode_readcounts
    params:
        pbc_id = lambda wildcards: project_df.get_puck_barcode_ids_and_files(
            project_id=wildcards.project_id, sample_id=wildcards.sample_id
        )[0],
        min_threshold = lambda wildcards:
            project_df.config.get_run_mode(list(get_run_modes_from_sample(
            wildcards.project_id, wildcards.sample_id).keys())[0]).variables['spatial_barcode_min_matches']
    output:
        puck_count_prealigned_barcode_matches_summary
    # at most 50% of CPU resources allocated to finding tiles
    # with 64 cores, in 'margaret', this takes 10 minutes to match 
    # ~3,5k targets of 2 M barcodes, against a query of ~600 M reads
    threads: max(workflow.cores * 0.5, 1)
    shell:
        "python {spacemake_dir}/snakemake/scripts/n_intersect_sequences.py"
        " --query {input.bc_counts}"
        " --query-plain-skip 1"
        " --query-plain-column 1"
        " --target {input.puck_barcode_files}"
        " --target-id {params.pbc_id}"
        " --target-column 'cell_bc'"
        " --summary-output {output}"
        " --min-threshold {params.min_threshold}"
        " --n-jobs {threads}"

rule create_top_barcode_whitelist:
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

rule create_spatial_barcode_file:
    input:
        unpack(get_puck_file),
        unpack(get_all_barcode_readcounts),
        puck_barcode_files_filtered
    output:
        parsed_spatial_barcodes,
        temp(parsed_spatial_barcodes_summary)
    shell:
        "python {spacemake_dir}/snakemake/scripts/n_intersect_sequences.py"
        " --query {input.bc_readcounts}"
        " --query-plain-skip 1"
        " --query-plain-column 1"
        " --target {input.barcode_file}"
        " --target-id {wildcards.puck_barcode_file_id}"
        " --target-column 'cell_bc'"
        " --output {output[0]}"
        " --summary-output {output[1]}"
        " --n-jobs {threads}"   
        " --chunksize 10000000"

rule create_spatial_barcode_whitelist:
    input: parsed_spatial_barcodes
    output: temp(spatial_barcodes)
    run:
        bc = pd.read_csv(input[0])
        bc = bc[['cell_bc']]
        # bc = bc.append({'cell_bc': 'NNNNNNNNNNNN'}, ignore_index=True)

        # save both the whitelist and the beads in a separate file
        bc[['cell_bc']].to_csv(output[0], header=False, index=False)

rule create_filtered_bc_summary:
    input: puck_count_prealigned_barcode_matches_summary
    output: temp(puck_barcode_files_filtered)
    params:
        min_threshold = lambda wildcards:
                project_df.config.get_run_mode(list(get_run_modes_from_sample(
                wildcards.project_id, wildcards.sample_id).keys())[0]).variables['spatial_barcode_min_matches']
    run:
        df_prealigned = pd.read_csv(input[0])
        df_filtered = df_prealigned[df_prealigned['matching_ratio'] > params.min_threshold]
        df_filtered.to_csv(output[0], index=False)
        
        
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
            project_id = wildcards.project_id,
            sample_id = wildcards.sample_id)['{cell}'],
        umi_tag = lambda wildcards: get_bam_tag_names(
            project_id = wildcards.project_id,
            sample_id = wildcards.sample_id)['{UMI}']
    threads: 1
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
        TMP_DIR={global_tmp} \
        {params.dge_extra_params}
        """

rule create_h5ad_dge:
    input:
        unpack(get_puck_file),
        unpack(get_raw_dge)
    # output here will either be n_beads=number, n_beads=spatial
    output: dge_out_h5ad, dge_out_h5ad_obs
    run:
        if wildcards.is_external == '.external':
            adata = load_external_dge(input['dge'])
        else:
            adata = dge_to_sparse_adata(
                input['dge'],
                input['dge_summary'])
        # attach barcodes
        if 'barcode_file' in input.keys() and wildcards.n_beads == 'spatial':
            adata = attach_barcode_file(adata, input['barcode_file'])
            adata = attach_puck(
                adata,
                project_df.get_puck(
                    project_id = wildcards.project_id,
                    sample_id = wildcards.sample_id,
                    return_empty=True
                )
            )

        # add 'cell_bc' name to index for same format as individual pucks
        # this also ensures compatibility with qc_sequencing_create_sheet.Rmd
        adata.write(output[0])
        adata.obs.to_csv(output[1])

rule create_mesh_spatial_dge:
    input:
        dge_spatial
    output:
        dge_spatial_mesh,
        dge_spatial_mesh_obs
    params:
        puck_data = lambda wildcards: project_df.get_puck_variables(
            project_id = wildcards.project_id,
            sample_id = wildcards.sample_id),
        pbf_metrics = lambda wildcards: project_df.get_puck_barcode_file_metrics(
            project_id = wildcards.project_id,
            sample_id = wildcards.sample_id,
            puck_barcode_file_id = wildcards.puck_barcode_file_id,
            polyA_adapter_trimmed = wildcards.polyA_adapter_trimmed)
    run:
        adata = sc.read(input[0])
        if wildcards.spot_distance_um == 'hexagon':
            mesh_type = 'hexagon'
            # if hexagon, this will be ignored
            spot_distance_um = 10.0
        else:
            mesh_type = 'circle'
            spot_distance_um = float(wildcards.spot_distance_um)
        adata = create_meshed_adata(adata,
            px_by_um = params['pbf_metrics']['px_by_um'],
            bead_diameter_um = params['puck_data']['spot_diameter_um'],
            spot_diameter_um = float(wildcards.spot_diameter_um),
            spot_distance_um = spot_distance_um,
            mesh_type = mesh_type)
        adata.write(output[0])
        adata.obs.to_csv(output[1])

rule puck_collection_stitching:
    input:
        unpack(checkpoint_puck_collection),
        puck_barcode_files_summary
    output:
        dge_spatial_collection,
        dge_spatial_collection_obs
    params:
        puck_data = lambda wildcards: project_df.get_puck_variables(
                                project_id = wildcards.project_id,
                                sample_id = wildcards.sample_id),
        puck_metadata = lambda wildcards: project_df.get_puck_barcode_ids_and_files(
            project_id=wildcards.project_id, sample_id=wildcards.sample_id
        ),
    run:
        # get the compatible puck_ids for the puck_files
        puck_barcode_file_ids = get_pbf_from_checkpoint_pucks(wildcards)

        _pc = puck_collection.merge_pucks_to_collection(
            # takes all input except the dge_out_done and puck_barcode_files
            input[:-1],
            puck_barcode_file_ids,
            params['puck_data']['coordinate_system'],
            "",
            "puck_id",
        )
        
        # x_pos and y_pos to be global coordinates
        _pc.obs['x_pos'] = _pc.obsm['spatial'][..., 0]
        _pc.obs['y_pos'] = _pc.obsm['spatial'][..., 1]

        _pc.write_h5ad(output[0])
    
        # add 'cell_bc' name to index for same format as individual pucks
        # this also ensures compatibility with qc_sequencing_create_sheet.Rmd
        df = _pc.obs
        df.index = np.arange(len(df))
        df.index.name = "cell_bc"
        # only get numeric columns, to avoid problems during summarisation
        # we could implement sth like df.A.str.extract('(\d+)')
        # to avoid losing information from columns that are not numeric
        df._get_numeric_data().to_csv(output[1])


# TODO: collapse this with previous rule so we have a single point where we create the dge_spatial_collection
rule puck_collection_stitching_meshed:
    input:
        unpack(checkpoint_puck_collection_mesh),
        puck_barcode_files_summary
    output:
        dge_spatial_collection_mesh,
        dge_spatial_collection_mesh_obs
    params:
        puck_data = lambda wildcards: project_df.get_puck_variables(
                                project_id = wildcards.project_id,
                                sample_id = wildcards.sample_id),
        puck_metadata = lambda wildcards: project_df.get_puck_barcode_ids_and_files(
            project_id=wildcards.project_id, sample_id=wildcards.sample_id
        ),
    run:
        puck_barcode_file_ids = get_pbf_from_checkpoint_pucks(wildcards)

        _pc = puck_collection.merge_pucks_to_collection(
            # takes all input except the dge_out_done and puck_barcode_files
            input[:-1],
            puck_barcode_file_ids,
            params['puck_data']['coordinate_system'],
            "",
            "puck_id",
        )
        
        # x_pos and y_pos to be global coordinates
        _pc.obs['x_pos'] = _pc.obsm['spatial'][..., 0]
        _pc.obs['y_pos'] = _pc.obsm['spatial'][..., 1]

        _pc.write_h5ad(output[0])

        # add 'cell_bc' name to index for same format as individual pucks
        # this also ensures compatibility with qc_sequencing_create_sheet.Rmd
        df = _pc.obs
        df.index = np.arange(len(df))
        df.index.name = "cell_bc"

        # only get numeric columns, to avoid problems during summarisation
        # we could implement sth like df.A.str.extract('(\d+)')
        # to avoid losing information from columns that are not numeric
        df._get_numeric_data().to_csv(output[1])

rule create_qc_sheet:
    input:
        unpack(get_qc_sheet_input_files),
        ribo_log=parsed_ribo_depletion_log
    params:
        sample_info = lambda wildcards: project_df.get_sample_info(
            wildcards.project_id, wildcards.sample_id),
        puck_variables = lambda wildcards:
            project_df.get_puck_variables(wildcards.project_id, wildcards.sample_id,
                return_empty=True),
        pbf_metrics = lambda wildcards: project_df.get_puck_barcode_file_metrics(
            project_id = wildcards.project_id,
            sample_id = wildcards.sample_id,
            puck_barcode_file_id = wildcards.puck_barcode_file_id_qc,
            polyA_adapter_trimmed=wildcards.polyA_adapter_trimmed),
        is_spatial = lambda wildcards:
            project_df.is_spatial(wildcards.project_id, wildcards.sample_id,
                puck_barcode_file_id=wildcards.puck_barcode_file_id_qc),
        run_modes = lambda wildcards: get_run_modes_from_sample(
            wildcards.project_id, wildcards.sample_id)
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
        is_spatial = lambda wildcards:
            project_df.is_spatial(wildcards.project_id, wildcards.sample_id,
                puck_barcode_file_id=wildcards.puck_barcode_file_id_qc),
        run_mode_variables = lambda wildcards:
            project_df.config.get_run_mode(wildcards.run_mode).variables
    script:
        'scripts/automated_analysis.py'

rule run_novosparc_denovo:
    input:
        automated_analysis_result_file
    output:
        novosparc_denovo_h5ad
    threads: 4
    shell:
        "python {spacemake_dir}/spatial/novosparc_integration.py"
        " --single_cell_dataset {input}"
        " --output {output}"

rule run_novosparc_with_reference:
    input:
        unpack(get_novosparc_with_reference_input_files),
        sc_adata=automated_analysis_result_file
    output:
        novosparc_with_reference_h5ad
    shell:
        "python {spacemake_dir}/spatial/novosparc_integration.py"
        " --single_cell_dataset {input.sc_adata}"
        " --spatial_dataset {input.st_adata}"
        " --output {output}"

rule create_automated_analysis_processed_data_files:
    input:
        automated_analysis_result_file
    output:
        **automated_analysis_processed_data_files
    params:
        is_spatial = lambda wildcards:
            project_df.is_spatial(wildcards.project_id, wildcards.sample_id,
                puck_barcode_file_id=wildcards.puck_barcode_file_id_qc),
    script:
        'scripts/automated_analysis_create_processed_data_files.py'
        
rule create_automated_report:
    input:
        **automated_analysis_processed_data_files,
    threads: 1
    output:
        automated_report
    params:
        run_mode_variables = lambda wildcards:
            project_df.config.get_run_mode(wildcards.run_mode).variables,
        puck_variables = lambda wildcards:
            project_df.get_puck_variables(wildcards.project_id, wildcards.sample_id,
                return_empty=True),
        pbf_metrics = lambda wildcards: project_df.get_puck_barcode_file_metrics(
            project_id = wildcards.project_id,
            sample_id = wildcards.sample_id,
            puck_barcode_file_id = wildcards.puck_barcode_file_id_qc,
            polyA_adapter_trimmed=wildcards.polyA_adapter_trimmed),
        is_spatial = lambda wildcards:
            project_df.is_spatial(wildcards.project_id, wildcards.sample_id,
                puck_barcode_file_id=wildcards.puck_barcode_file_id_qc),
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

rule create_barcode_files_matching_summary:
    input:
        unpack(checkpoint_barcode_files),
        bc_summary=puck_barcode_files_filtered
    output:
        puck_barcode_files_summary
    params:
        puck_variables = lambda wildcards:
            project_df.get_puck_variables(wildcards.project_id, wildcards.sample_id,
                return_empty=True),
        run_mode_variables = lambda wildcards:
            project_df.config.get_run_mode(list(get_run_modes_from_sample(
            wildcards.project_id, wildcards.sample_id).keys())[0]).variables,
    run:
        import os
        out_df = pd.DataFrame(columns=[
            'puck_barcode_file_id',
            'puck_barcode_file',
            'parsed_barcode_file',
            'n_barcodes',
            'n_matching',
            'matching_ratio', 
            'x_pos_min_px',
            'x_pos_max_px',
            'y_pos_min_px',
            'y_pos_max_px',
            'px_by_um'])

        df_puck_barcode_files_filtered = pd.read_csv(input['bc_summary'])

        if ('puck_barcode_files' in input.keys() and
            'parsed_spatial_barcode_files' in input.keys()):
            for pbf_id, pbf, parsed_barcode_file in zip(
                    df_puck_barcode_files_filtered['puck_barcode_file_id'],
                    input['puck_barcode_files'],
                    input['parsed_spatial_barcode_files'],
                ):

                pbf_df = parse_barcode_file(pbf)
                n_barcodes = pbf_df.shape[0]
                n_matching = pd.read_csv(parsed_barcode_file).shape[0]
                matching_ratio = round(float(n_matching)/n_barcodes, 2)

                # calculate puck metrics
                x_pos_min_px = pbf_df.x_pos.min()
                x_pos_max_px = pbf_df.x_pos.max()
                y_pos_min_px = pbf_df.y_pos.min()
                y_pos_max_px = pbf_df.y_pos.max()

                px_by_um = (x_pos_max_px - x_pos_min_px) 
                px_by_um = px_by_um / params['puck_variables']['width_um']
                
                out_df = pd.concat([out_df, pd.DataFrame({
                    'puck_barcode_file_id': [pbf_id],
                    'puck_barcode_file': [pbf],
                    'parsed_barcode_file': [parsed_barcode_file],
                    'n_barcodes': [n_barcodes],
                    'n_matching': [n_matching],
                    'matching_ratio': [matching_ratio],
                    'x_pos_min_px': [x_pos_min_px],
                    'x_pos_max_px': [x_pos_max_px],
                    'y_pos_min_px': [y_pos_min_px],
                    'y_pos_max_px': [y_pos_max_px],
                    'px_by_um': [px_by_um],
                })], ignore_index=True, sort=False)

        out_df.to_csv(output[0], index=False)

###################
# Puck checkpoint #
###################
checkpoint checkpoint_barcodes:
    input:
        bc_summary_file=puck_barcode_files_filtered
    output:
        bc_summaries=temp(directory([bc_out_dir]))
    run:
        os.mkdir(output[0])
        barcodes_df = pd.read_csv(input.bc_summary_file)
        for p in barcodes_df['puck_barcode_file_id'].tolist():
            with open(output[0] + f"/{p}.chk", "w") as out:
                out.write("")

checkpoint checkpoint_pucks:
    input:
        bc_summary_file=puck_barcode_files_summary
    output:
        dge_pointers=temp(directory(dge_out_dir))
    run:
        os.mkdir(output.dge_pointers)
        barcodes_df = pd.read_csv(input.bc_summary_file)
        for p in barcodes_df['puck_barcode_file_id'].tolist():
            with open(output.dge_pointers + f"/{p}.chk", "w") as out:
                out.write("")

rule aggregate_dge:
    input:
        unpack(checkpoint_puck_files),
    output:
        temp(touch(dge_out_done))
