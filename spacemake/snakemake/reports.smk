"""
Snakemake rules for generating diagnostic metrics and plots.
"""

rule bamstats:
    """
    Gather data for read-length histograms and for cell-barcode
    UMI count knee-plots.
    """
    input: complete_data_root + "/{bamname}.bam"
    output:
        stats=stats_dir + "/{bamname}.bc_stats.tsv",
        rlen=stats_dir + "/{bamname}.readlen.tsv"
    log: log_dir + "/{bamname}.bamstats.log"
    shell:
        "python {bin_dir}/bamstats.py {input} "
        "  --sample={wildcards.sample_id} "
        "  --log-level={log_level} "
        "  --log-file={log} "
        "  --out-stats={output.stats} "
        "  --out-length={output.rlen} "


rule complexity_sampling:
    """
    Sub-sample mapped bam files to get UMI counts as function of sequencing read number
    """
    input: complete_data_root + "/{bamname}.bam"
    output: stats_dir + "/{bamname}.complexity.tsv"
    log: log_dir + "/{bamname}.complexity.log"
    shell:
        "python {bin_dir}/complexity.py {input} "
        "  --sample={wildcards.sample_id} "
        "  --log-level={log_level} "
        "  --log-file={log} "
        "  --out-tsv={output} "


def get_all_bamstats_for_sample(wildcards):
    raw = wc_fill(stats_dir  + "/unaligned_bc_tagged.readlen.tsv", wildcards),
    trimmed = wc_fill(stats_dir + "/cutadapt.csv", wildcards),
    mapped = []
    not_mapped = []
    complexity = []
    for bampath in get_all_mapped_bams(wildcards)['mapped_bams']:
        dirname, basename = os.path.split(bampath)
        mapped.append(f"{dirname}/stats/{basename.replace('.bam', '.readlen.tsv')}")
        not_mapped.append(f"{dirname}/stats/not_{basename.replace('.bam', '.readlen.tsv')}")
        complexity.append(f"{dirname}/stats/{basename.replace('.bam', '.complexity.tsv')}")

    return { 'raw': raw, 'trimmed':trimmed, 'mapped': mapped, 'not_mapped': not_mapped, 'complexity' : complexity}


overview_tsv = stats_dir + "/{sample_id}.overview.tsv"
overview_pdf = plots_dir + "/{sample_id}.overview.pdf"

rule overview:
    """
    Combine statistics from raw, adapter-trimmed, mapped and unmapped reads into one overview
    """
    input:
        unpack(get_all_bamstats_for_sample)
    output:
        tsv = overview_tsv,
        pdf = overview_pdf,
    log: log_dir + '/overview_plot.log'
    params:
        mapped_list = lambda wildcards, input: " ".join(input.mapped),
        not_mapped_list = lambda wildcards, input: " ".join(input.not_mapped),
        map_strategy = lambda wildcards: SAMPLE_MAP_STRATEGY[(wildcards.project_id, wildcards.sample_id)],
    shell:
        "python {bin_dir}/plot_overview.py "
        "  --sample={wildcards.sample_id} "
        "  --log-level={log_level} "
        "  --log-file={log} "
        "  --raw={input.raw} "
        "  --trimmed={input.trimmed} "
        "  --map-strategy='{params.map_strategy}' "
        "  --mapped {params.mapped_list} "
        "  --not-mapped {params.not_mapped_list} "
        "  --out-tsv {output.tsv} "
        "  --out-pdf {output.pdf}"

dge_legacy_params = "{dge_type}{dge_cleaned}{polyA_adapter_trimmed}{mm_included}.{n_beads}_beads_{puck_barcode_file_id}{is_external}"
# dge_type, cleaned -> out
# polyA_adapter-trimmed -> out
# mm_included -> stays
# {n_beads}_beads_ -> name the subsetting strategy -> top10k, top100k, spatial
# puck_barcode_file -> essential (tiles)

dge_metrics_plot = plots_dir + '/{sample_id}' + dge_legacy_params + '.{ref}.metrics.pdf'
dge_knee_plot = plots_dir + '/{sample_id}' + dge_legacy_params + '.{ref}.knee.pdf'

rule dge_stats:
    input: dge_out_h5ad
    output: 
        metrics=dge_metrics_plot,
        knee=dge_knee_plot
    log: log_dir + '/{sample_id}' + dge_legacy_params + '.{ref}.dge_stats_plot.log'
    shell:
        "python {bin_dir}/plot_dge_stats.py "
        "  --log-level={log_level} "
        "  --log-file={log} "
        "  --debug={log_debug} "
        "  --sample={wildcards.sample_id} "
        "  --reference={wildcards.ref} "
        "  --out-metrics={output.metrics} "
        "  --out-knee={output.knee} "
        "  {input}"

base_freqs_plot = plots_dir + '/{sample_id}.base_freqs.pdf'

rule make_base_freqs_plot:
    input: tagged_bam
    output: base_freqs_plot
    shell:
        "python {bin_dir}/plot_base_freqs.py "
        "  --log-level={log_level} "
        "  --log-file={log} "
        "  --debug={log_debug} "
        "  --sample={wildcards.sample_id} "
        "  --out-pdf={output} "
        "  {input}"

complexity_plot = plots_dir + '/{sample_id}.complexity.pdf'

rule make_complexity_plot:
    input:
        unpack(get_all_bamstats_for_sample)
    output: complexity_plot
    shell:
        "python {bin_dir}/plot_complexity.py "
        "  --log-level={log_level} "
        "  --log-file={log} "
        "  --debug={log_debug} "
        "  --sample={wildcards.sample_id} "
        "  --out-plot={output} "
        "  --cutadapt={input.trimmed} "
        "  {input.complexity}"

def get_plots():
    out_files = []
    for (project_id, sample_id), row in project_df.df.iterrows():
        # barcode & UMI freq stats
        out_files.append(wc_fill(base_freqs_plot, dotdict(project_id=project_id, sample_id=sample_id)))

        # complexity sub-sampling plot
        out_files.append(wc_fill(complexity_plot, dotdict(project_id=project_id, sample_id=sample_id)))

        # pre-processing and mapping overview
        out_files.append(wc_fill(overview_tsv, dotdict(project_id=project_id, sample_id=sample_id)))
        out_files.append(wc_fill(overview_pdf, dotdict(project_id=project_id, sample_id=sample_id)))

        # dge stats
        context = dotdict(
            project_id=project_id,
            sample_id=sample_id,
            dge_type='.all',
            dge_cleaned='',
            polyA_adapter_trimmed='.polyA_adapter_trimmed',
            mm_included='',
            n_beads=100000,
            puck_barcode_file_id='no_spatial_data',
            is_external='',
            ref='{ref}'
        )
        out_files.append(
            expand(dge_metrics_plot.format(**context), ref=REF_NAMES[(project_id, sample_id)])
        )
        out_files.append(
            expand(dge_metrics_plot.format(**context), ref=REF_NAMES[(project_id, sample_id)])
        )
    return out_files


register_module_output_hook(get_plots, "reports.smk")
