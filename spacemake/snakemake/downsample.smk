#########
# about #
#########
__version__ = '0.1.0'
__author__ = ['Nikos Karaiskos', 'Tamas Ryszard Sztanka-Toth']
__licence__ = 'GPL'
__email__ = ['nikolaos.karaiskos@mdc-berlin.de', 'tamasryszard.sztanka-toth@mdc-berlin.de']

# first create downsampling for 10, 20 .. 90
downsampled_ratios = range(10,100,10)

downsample_root = illumina_root + '/downsampled_data'
downsampled_sample_root = downsample_root + '/{ratio}'
downsampled_bam = downsampled_sample_root + '/final_downsampled_{ratio}{polyA_adapter_trimmed}.bam'
downsampled_bam_mm_included_pipe = downsampled_sample_root + '/final_downsampled_{ratio}' + bam_mm_included_pipe_suffix

downsampled_readcounts = downsampled_sample_root + '/out_readcounts' + barcode_readcounts_suffix
downsampled_top_barcodes = downsampled_sample_root + '/topBarcodes' + top_barcodes_suffix
downsampled_top_barcodes_clean = downsampled_sample_root + '/topBarcodesClean' + top_barcodes_suffix

# dges
downsampled_dge_root = downsampled_sample_root + '/dge'
downsampled_dge_out_prefix = downsampled_dge_root + '/downsampled_dge'
downsampled_dge_out = downsampled_dge_out_prefix + dge_out_suffix + '.txt.gz'
downsampled_dge_out_summary = downsampled_dge_out_prefix + dge_out_suffix + '.summary.gz'

#downsample_qc_sheet = downsampled_sample_root  + '/qc_sheet_{united_sample}_{puck}_downsampled_{ratio}.html'

downsample_saturation_analysis = downsample_root + '/{project_id}_{sample_id}_{run_mode}_saturation_analysis.html'

rule downsample_bam:
    input:
        final_bam
    output:
        downsampled_bam
    params:
        downsample_dir = downsampled_sample_root
    threads: 4
    shell:
        """
        mkdir -p {params.downsample_dir}

        sambamba view -o {output} -f bam -t {threads} -s 0.{wildcards.ratio} {input}
        """

rule downsampled_filter_mm_reads:
    input:
        downsampled_bam
    output:
        pipe(downsampled_bam_mm_included_pipe)
    shell:
        """
        python {repo_dir}/scripts/filter_mm_reads.py \
            --in-bam {input} \
            --out-bam {output}
        """

rule downsample_bam_tag_histogram:
    input:
        downsampled_bam
    output:
        downsampled_readcounts
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

rule create_downsampled_top_barcodes_file:
    input:
        downsampled_readcounts
    output:
        downsampled_top_barcodes
    shell:
        "zcat {input} | cut -f2 | head -{wildcards.n_beads} > {output}"
        
rule downsampled_clean_top_barcodes:
    input:
        downsampled_top_barcodes
    output:
        downsampled_top_barcodes_clean
    script:
        'scripts/clean_top_barcodes.py'

def get_downsampled_top_barcodes(wildcards):
    if wildcards.dge_cleaned == "":
        return {"top_barcodes": downsampled_top_barcodes}
    else:
        return {'top_barcodes': downsampled_top_barcodes_clean}

def get_downsampled_mapped_final_bam(wildcards):
    if wildcards.mm_included == '.mm_included':
        return {'reads': downsampled_bam_mm_included_pipe}
    else:
        return {'reads': downsampled_bam}

rule create_downsampled_dge:
    input:
        unpack(get_downsampled_mapped_final_bam),
        unpack(get_downsampled_top_barcodes) 
    output:
        downsampled_dge=downsampled_dge_out,
        downsampled_dge_summary=downsampled_dge_out_summary
    params:
        downsampled_dge_root = downsampled_dge_root,
        downsampled_dge_extra_params = lambda wildcards: get_dge_extra_params(wildcards),
        cell_barcode_tag = lambda wildcards: get_bam_tag_names(
            project_id = wildcards.project,
            sample_id = wildcards.sample)['{cell}'],
        umi_tag = lambda wildcards: get_bam_tag_names(
            project_id = wildcards.project,
            sample_id = wildcards.sample)['{UMI}']
    shell:
        """
        mkdir -p {params.downsampled_dge_root}

        {dropseq_tools}/DigitalExpression \
        -m 16g \
        I= {input.reads}\
        O= {output.downsampled_dge} \
        SUMMARY= {output.downsampled_dge_summary} \
        CELL_BC_FILE={input.top_barcodes} \
        CELL_BARCODE_TAG={params.cell_barcode_tag} \
        MOLECULAR_BARCODE_TAG={params.umi_tag} \
        TMP_DIR={temp_dir} \
        {params.downsampled_dge_extra_params}
        """

#rule create_downsample_qc_sheet:
#    input:
#        star_log = united_star_log,
#        reads_type_out=united_reads_type_out,
#        parameters_file=united_qc_sheet_parameters_file,
#        read_counts=united_barcode_readcounts,
#        dge_all_summary = downsampled_dge_root + '/downsampled_dge_all_summary.txt'
#    output:
#        downsample_qc_sheet
#    script:
#        "analysis/qc_sequencing_create_sheet.Rmd"


def get_saturation_analysis_input(wildcards):
    # create dictionary with the right downsampling files where  the key
    files = {}

    for ratio in downsampled_ratios:
        # dge_files contains dge/summary file paths per run_mode
        dge_files = get_dges_from_project_sample(
            project_id = wildcards.project,
            sample_id = wildcards.sample,
            dge_out_pattern = downsampled_dge_out,
            dge_out_summary_pattern = downsampled_dge_out_summary,
            ratio=ratio)

        for key, file_path in dge_files.items():
            files[f'downsample_{ratio}_{key}'] = file_path

    #dge_summaries['downsampled_100'] = expand(dge_all_summary,
    #    united_project = wildcards.united_project,
    #    united_sample = wildcards.united_sample)
    #

    return files

rule create_saturation_analysis:
    input:
        unpack(get_saturation_analysis_input)
    output:
        downsample_saturation_analysis
    script:
        "scripts/saturation_analysis.Rmd"
