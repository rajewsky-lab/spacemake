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

from spacemake.preprocess import dge_to_sparse_adata, attach_barcode_file,\
    parse_barcode_file, load_external_dge, attach_puck_variables
from spacemake.spatial import create_meshed_adata
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
        for out in hook():
            print(f"output provided by '{module}' module (via '{hook.__name__}'): '{out}'")
            outputs.append(out)
    
    return outputs

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
    umi_cutoff = '\d+',
    dge_cleaned='|\.cleaned',
    dge_type = '|'.join(dge_types),
    pacbio_ext = 'fq|fastq|bam',
    polyA_adapter_trimmed = '|\.polyA_adapter_trimmed',
    mm_included = '|\.mm_included',
    n_beads = '[0-9]+|spatial|external',
    is_external = '|\.external',
    spot_diameter_um = '[0-9]+',
    spot_distance_um = '[0-9]+|hexagon',
    data_root_type = 'complete_data|downsampled_data',
    downsampling_percentage = '\/[0-9]+|',
    puck_barcode_file_id = '[^.]+'

#############
# Main rule #
#############
rule run_analysis:
    input:
        # create fastq
        unpack(
            lambda wildcards: get_output_files(
                    fastqc_pattern, ext = fastqc_ext, mate=['1', '2'],
                    data_root_type = 'complete_data', downsampling_percentage = '',
                    filter_merged=True) 
                if config['with_fastqc'] else []
        ),
        unpack(get_all_dges),
        # this will also create the clean dge
        get_output_files(automated_report, data_root_type = 'complete_data',
            downsampling_percentage='', puck_barcode_file_matching_type='spatial_matching'),
        get_output_files(qc_sheet, data_root_type = 'complete_data',
            downsampling_percentage='', run_on_external=False,
            puck_barcode_file_matching_type='spatial_matching'),
        # finally, everything registered via register_module_output_hook()
        get_module_outputs(),


rule get_whitelist_barcodes:
    input:
        get_output_files(barcode_readcounts,
            data_root_type='complete_data',
            downsampling_percentage='',
            run_on_external=False,
        ),
        get_output_files(puck_barcode_files_summary,
            data_root_type = 'complete_data',
            downsampling_percentage='', run_on_external=False)

##############
# DOWNSAMPLE #
##############
rule downsample:
    input:
        get_output_files(downsample_saturation_analysis,
            samples = config['samples'],
            projects = config['projects'])

#############
# NOVOSPARC #
#############
#rule novosparc:
#    input:
#        *get_novosparc_input_files(config)

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
        unassigned = unassigned,
        bc_stats = reverse_reads_mate_1.replace(reads_suffix, ".bc_stats.tsv")
    log:
        reverse_reads_mate_1.replace(reads_suffix, ".preprocessing.log")
    threads: 4
    shell:
        "python {spacemake_dir}/preprocess/cmdline.py "
        "--sample={wildcards.sample_id} "
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
            project_id = wildcards.project_id,
            sample_id = wildcards.sample_id)['{cell}'],
    shell:
        """
        {dropseq_tools}/BamTagHistogram -m 32g \
        I= {input} \
        O= {output} \
        TAG={params.cell_barcode_tag} \
        READ_MQ=0
        """

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
        unpack(get_all_barcode_readcounts)
    output:
        parsed_spatial_barcodes
    run:
        # load all readcounts
        bc_readcounts=[pd.read_table(bc_rc, skiprows=1,
            names=['read_n', 'cell_bc']) for bc_rc in input['bc_readcounts']]

        # join them together
        bc_readcounts = pd.concat(bc_readcounts)

        # remove duplicates
        bc_readcounts.drop_duplicates(subset='cell_bc', keep='first',
            inplace=True)

        # load barcode file and parse it
        bc = parse_barcode_file(input[0])
        bc.reset_index(level=0, inplace=True)
        # inner join to get rid of barcode without any data
        bc = pd.merge(bc, bc_readcounts, how='inner', on='cell_bc')
        bc = bc[['cell_bc', 'x_pos', 'y_pos']]

        bc.to_csv(output[0], index=False)

rule create_spatial_barcode_whitelist:
    input: 
        parsed_bc = parsed_spatial_barcodes,
        bc_summary = puck_barcode_files_summary
    output: temp(spatial_barcodes)
    params:
        run_mode_variables = lambda wildcards:
            project_df.config.get_run_mode(list(get_run_modes_from_sample(
            wildcards.project_id, wildcards.sample_id).keys())[0]).variables
    run:
        bc = pd.read_csv(input['parsed_bc'])
        bc = bc[['cell_bc']]
        bc = bc.append({'cell_bc': 'NNNNNNNNNNNN'}, ignore_index=True)

        # save both the whitelist and the beads in a separate file
        bc[['cell_bc']].to_csv(output[0], header=False, index=False)

        # update the project df (but do not dump) to only process some tiles above the threshold
        barcodes_summary = pd.read_csv(input['bc_summary'])
        barcodes_summary_filter = barcodes_summary[
            barcodes_summary.matching_ratio >= params['run_mode_variables']['spatial_barcode_min_matches']
        ]

        _puck_barcode_files = barcodes_summary_filter['puck_barcode_file'].values
        _puck_barcode_files_id = barcodes_summary_filter['puck_barcode_file_id'].values

        project_df.add_update_sample(
            action='update',
            project_id=wildcards.project_id,
            sample_id=wildcards.sample_id,
            puck_barcode_file=_puck_barcode_files,
            puck_barcode_file_id=_puck_barcode_files_id,
        )
        
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
    # at most 8 dges will be created the same time
    threads: max(workflow.cores * 0.125, 1)
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
            adata = attach_puck_variables(
                adata,
                project_df.get_puck_variables(
                    project_id = wildcards.project_id,
                    sample_id = wildcards.sample_id,
                    return_empty=True
                )
            )
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
            puck_barcode_file_id = wildcards.puck_barcode_file_id)
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
        unpack(get_puck_collection_stitching_input)
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
        _pc = puck_collection.merge_pucks_to_single_adata(
            input,
            params['puck_metadata'][0],
            params['puck_data']['coordinate_system'],
            params['puck_data']['puck_id_regex'],
            "puck_id",
            no_reset_index,
            no_transform,
            merge_output,
            join_output,
            output,
        )

        _pc.write_h5ad(output[0])
        _pc.obs.to_csv(output[1])

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
            puck_barcode_file_id = wildcards.puck_barcode_file_id),
        is_spatial = lambda wildcards:
            project_df.is_spatial(wildcards.project_id, wildcards.sample_id,
                puck_barcode_file_id=wildcards.puck_barcode_file_id),
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
                puck_barcode_file_id=wildcards.puck_barcode_file_id),
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
                puck_barcode_file_id=wildcards.puck_barcode_file_id),
    script:
        'scripts/automated_analysis_create_processed_data_files.py'
        
rule create_automated_report:
    input:
        unpack(get_parsed_puck_file),
        **automated_analysis_processed_data_files,
    # spawn at most 4 automated analyses
    threads: max(workflow.cores / 8, 1)
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
            puck_barcode_file_id = wildcards.puck_barcode_file_id),
        is_spatial = lambda wildcards:
            project_df.is_spatial(wildcards.project_id, wildcards.sample_id,
                puck_barcode_file_id=wildcards.puck_barcode_file_id),
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
        unpack(get_barcode_files_matching_summary_input)
    output:
        puck_barcode_files_summary
    params:
        pbf_ids = lambda wildcards: project_df.get_puck_barcode_ids_and_files(
            project_id = wildcards.project_id,
            sample_id = wildcards.sample_id)[0],
        puck_variables = lambda wildcards:
            project_df.get_puck_variables(wildcards.project_id, wildcards.sample_id,
                return_empty=True)
    run:
        out_df = pd.DataFrame(columns=[
            'puck_barcode_file_id',
            'n_barcodes',
            'n_matching',
            'matching_ratio', 
            'x_pos_min_px',
            'x_pos_max_px',
            'y_pos_min_px',
            'y_pos_max_px',
            'px_by_um'])

        if ('puck_barcode_files' in input.keys() and
            'parsed_spatial_barcode_files' in input.keys()):
            for pbf_id, pbf, parsed_barcode_file in zip(
                    params['pbf_ids'],
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
                
                out_df = out_df.append({
                    'puck_barcode_file_id': pbf_id,
                    'puck_barcode_file': parsed_barcode_file,
                    'n_barcodes': n_barcodes,
                    'n_matching': n_matching,
                    'matching_ratio': matching_ratio,
                    'x_pos_min_px': x_pos_min_px,
                    'x_pos_max_px': x_pos_max_px,
                    'y_pos_min_px': y_pos_min_px,
                    'y_pos_max_px': y_pos_max_px,
                    'px_by_um': px_by_um,
                }, ignore_index=True)

            out_df.to_csv(output[0], index=False)
        else:
            # save empty file
            out_df.to_csv(output[0], index=False)