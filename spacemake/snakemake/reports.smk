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
        stats=reports_dir + "/{bamname}.bc_stats.tsv",
        rlen=reports_dir + "/{bamname}.readlen.tsv"
    shell:
        "python {repo_dir}/scripts/bamstats.py {input} --sample={wildcards.sample_id} --out-stats={output.stats} --out-length={output.rlen}"


rule complexity_sampling:
    """
    Sub-sample mapped bam files to get UMI counts as function of sequencing read number
    """
    input: complete_data_root + "/{bamname}.bam"
    output: reports_dir + "/{bamname}.complexity.tsv"
    shell:
        "python {repo_dir}/scripts/complexity.py {input} --sample={wildcards.sample_id} --out-tsv={output}"


def get_all_bamstats_for_sample(wildcards):
    raw = wc_fill(reports_dir  + "/unaligned_bc_tagged.readlen.tsv", wildcards),
    trimmed = wc_fill(reports_dir + "/cutadapt.csv", wildcards),
    mapped = []
    not_mapped = []
    complexity = []
    for bampath in get_all_mapped_bams(wildcards)['mapped_bams']:
        dirname, basename = os.path.split(bampath)
        mapped.append(f"{dirname}/reports/{basename.replace('.bam', '.readlen.tsv')}")
        not_mapped.append(f"{dirname}/reports/not_{basename.replace('.bam', '.readlen.tsv')}")
        complexity.append(f"{dirname}/reports/{basename.replace('.bam', '.complexity.tsv')}")

    return { 'raw': raw, 'trimmed':trimmed, 'mapped': mapped, 'not_mapped': not_mapped, 'complexity' : complexity}


overview_tsv = reports_dir + "/{sample_id}.overview.tsv"
overview_pdf = reports_dir + "/{sample_id}.overview.pdf"

rule overview:
    """
    Combine statistics from raw, adapter-trimmed, mapped and unmapped reads into one overview
    """
    input:
        unpack(get_all_bamstats_for_sample)
    output:
        tsv = overview_tsv,
        pdf = overview_pdf,
    params:
        mapped_list = lambda wildcards, input: " ".join(input.mapped),
        not_mapped_list = lambda wildcards, input: " ".join(input.not_mapped),
        map_strategy = lambda wildcards: SAMPLE_MAP_STRATEGY[(wildcards.project_id, wildcards.sample_id)],
    shell:
        "python {repo_dir}/scripts/overview_plot.py "
        " --raw={input.raw} "
        " --sample={wildcards.sample_id} "
        " --trimmed={input.trimmed} "
        " --map-strategy='{params.map_strategy}' "
        " --mapped {params.mapped_list} "
        " --not-mapped {params.not_mapped_list} "
        " --out-tsv {output.tsv} "
        " --out-pdf {output.pdf}"


def get_overview_reports():
    out_files = []
    for (project_id, sample_id), row in project_df.df.iterrows():
        out_files.append(wc_fill(overview_tsv, dotdict(project_id=project_id, sample_id=sample_id)))
        out_files.append(wc_fill(overview_pdf, dotdict(project_id=project_id, sample_id=sample_id)))

    return out_files

register_module_output_hook(get_overview_reports, "reports.smk")
