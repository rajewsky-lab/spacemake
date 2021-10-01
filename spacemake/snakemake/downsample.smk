#########
# about #
#########
__version__ = '0.1.0'
__author__ = ['Nikos Karaiskos', 'Tamas Ryszard Sztanka-Toth']
__licence__ = 'GPL'
__email__ = ['nikolaos.karaiskos@mdc-berlin.de', 'tamasryszard.sztanka-toth@mdc-berlin.de']

# first create downsampling for 10, 20 .. 90
downsampled_ratios = range(10,100,10)

rule downsample_bam:
    input:
        unpack(get_final_bam)
    output:
        downsampled_bam
    params:
        downsample_dir = downsampled_data_prefix,
        ratio = lambda wildcards: wildcards.downsampling_percentage[1:]
    threads: 4
    shell:
        """
        mkdir -p {params.downsample_dir}

        sambamba view -o {output} -f bam -t {threads} \
            -s 0.{params.ratio} {input}
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

def get_saturation_analysis_input(wildcards):
    # create dictionary with the right downsampling files where  the key
    files = {}

    run_modes = get_run_modes_from_sample(wildcards.project, wildcards.sample)

    for run_mode in run_modes:
        for ratio in downsampled_ratios:
        # dge_files contains dge/summary file paths per run_mode
            dge_summary = get_dge_from_run_mode(
                project_id = wildcards.project,
                sample_id = wildcards.sample,
                run_mode = run_mode,
                data_root_type = 'downsampled_data',
                downsampling_percentage = '/' + str(ratio))['dge_summary']

            files[f'downsampled_dge_summary.{run_mode}.{ratio}'] = dge_summary

        files[f'downsampled_dge_summary.{run_mode}.100'] = get_dge_from_run_mode(
            project_id = wildcards.project,
            sample_id = wildcards.sample,
            run_mode = run_mode,
            data_root_type = 'complete_data',
            downsampling_percentage = '')['dge_summary']

    return files

rule create_saturation_analysis:
    input:
        unpack(get_saturation_analysis_input)
    output:
        downsample_saturation_analysis
    params:
        sample_info = lambda wildcards: project_df.get_sample_info(
            wildcards.project, wildcards.sample),
        run_modes = lambda wildcards: get_run_modes_from_sample(
            wildcards.project, wildcards.sample)
    script:
        "scripts/saturation_analysis.Rmd"
