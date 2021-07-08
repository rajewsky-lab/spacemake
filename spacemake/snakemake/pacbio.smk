# pacbio_projects, pacbio_samples, pacbio_ext = glob_wildcards(pacbio_fq)
pacbio_projects, pacbio_samples, pacbio_ext = [], [], []
pacbio_reports = []
pacbio_rRNA_counts = []
pacbio_stats = []

byo_package = '/data/rajewsky/projects/slide_seq/.repos/pb_annotate/byo'
pythonpath_prepend = 'PYTHONPATH=$PYTHONPATH:' + byo_package
blocks_fa_file = '/data/rajewsky/projects/slide_seq/.repos/pb_annotate/blocks.fa'

for project, sample in zip(pacbio_projects, pacbio_samples):
    pacbio_reports = pacbio_reports + expand(pacbio_report, project=project, sample=sample)
    pacbio_rRNA_counts = pacbio_rRNA_counts + expand(pacbio_rRNA_out, project=project, sample=sample)
    pacbio_stats = pacbio_stats + expand(pacbio_stats_file, project=project, sample=sample)

bin = "/data/rajewsky/projects/slide_seq/.repos/pb_annotate/pb_annotate"

rule pacbio:
    input: 
        reports=pacbio_reports,
        rRNA_counts=pacbio_rRNA_counts,
        stats=pacbio_stats
    output:
        overview=pacbio_overview,
        bead_overview=pacbio_bead_overview,
        csv=pacbio_overview_csv
    shell:
        "{pythonpath_prepend} python {bin} overview --rRNA-same-place --multi-page --breakdown {output.bead_overview} --output {output.overview} --csv-out {output.csv} {input.stats}"

rRNA_index = "/data/rajewsky/indices/rRNA_hsa_bwa_0.7.17/rRNA_hsa.fa"
# SAM Flag 2308 means, "not primary", "supplementary", or "unmapped". 
# We kick these out for counting rRNA hits

rule rRNA_fq:
    input: pacbio_fq
    output: pacbio_rRNA_out
    threads: 4
    shell: "bwa mem -t {threads} -x pacbio {rRNA_index} {input} | samtools view -F 2308 /dev/stdin | wc -l > {output}"


rule annotate_fq:
    input: pacbio_fq
    output:
        report=pacbio_report,
        stats=pacbio_stats_file,
        run_summary=pacbio_run_summary
    shell: "{pythonpath_prepend} python {bin} scan --deletions --report {output.report} --blocks {blocks_fa_file} --summary {output.stats} {input}  > {output.run_summary}"
