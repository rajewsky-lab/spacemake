#########
# about #
#########
__version__ = '0.2'
__author__ = ['Marvin Jens', 'Tamas Ryszard Sztanka-Toth']
__email__ = ['marvin.jens@mdc-berlin.de', 'tamasryszard.sztanka-toth@mdc-berlin.de']

lr_root = project_dir + "/processed_data/{sample_id}/longread"
lr_cache_dir = lr_root + "/cache/"
lr_ann_dir = lr_root + "/annotation/"
lr_stats_dir = lr_root + "/stats/"
lr_report_dir = lr_root + "/reports/"
lr_examples_dir = lr_root + "/examples/"
lr_cDNA_dir = lr_root + "/cDNA/"

# targets
lr_ann = lr_ann_dir + "{sample_id}.annotation.tsv"
lr_stats = lr_stats_dir + "{sample_id}.stats.tsv"
lr_report = lr_report_dir + "{sample_id}.donuts.pdf"
lr_report_stats = lr_stats_dir + "{sample_id}.report.tsv"
lr_edits = lr_report_dir + "{sample_id}.oligo_edits.pdf"
lr_cDNA = lr_cDNA_dir + "{sample_id}.fa"
lr_cDNA_log = lr_cDNA_dir + "{sample_id}.log"
lr_cDNA_oligo_analysis = lr_cDNA_dir + "{sample_id}.oligo_analysis.csv"
lr_cDNA_bam = lr_cDNA_dir + "{sample_id}.bam"
lr_examples = lr_examples_dir + "{sample_id}.txt"

lr_overview_dir = os.path.join(config['root_dir'], 'longread_overview/')
lr_overview_pdf = lr_overview_dir + 'fidelity.pdf'
lr_overview_csv = lr_overview_dir + 'overview.csv'

LR_RAW_FILES = {}
LR_SIGNATURE = {}
LR_REPORT_STATS = []
def get_longread_output(project_df=None, config=None, **kw):
    """
    This function is called from main.smk at least once 
    to determine which output files need to be generated
    from longread longread analysis.
    We use this opportunity to populate LR_RAW_FILES
    """
    out_files = []
    for index, row in project_df.df.iterrows():
        # for run_mode in row["run_mode"]:
        #     run_mode_variables = project_df.config.get_run_mode(run_mode).variables
            if row.longreads:
                LR_REPORT_STATS.extend(
                    expand(lr_report_stats, project_id=index[0], sample_id=index[1])
                )
                out_files += \
                expand(
                    lr_report,
                    project_id=index[0],
                    sample_id=index[1],
                ) + \
                expand(
                    lr_edits,
                    project_id=index[0],
                    sample_id=index[1],
                ) + \
                expand(
                    lr_cDNA_bam,
                    project_id=index[0],
                    sample_id=index[1],
                ) + \
                expand(
                    lr_cDNA_oligo_analysis,
                    project_id=index[0],
                    sample_id=index[1],
                )

                LR_RAW_FILES[index[1]] = row.longreads
                LR_SIGNATURE[index[1]] = row.longread_signature

    # if we have any longread analysis, generate an overview plot
    if out_files:
        out_files.append(lr_overview_pdf)

    return out_files

register_module_output_hook(get_longread_output, "longread.smk")

def get_args(wc):
    args = f""" \
    --cache={lr_cache_dir} \
    --annotation-out={lr_ann_dir} \
    --stats-out={lr_stats_dir} \
    --report-out={lr_report_dir} \
    --examples-out={lr_examples_dir} \
    --sample={wc.sample_id} \
    --signature={LR_SIGNATURE[wc.sample_id]} \
    """.format(sample_id=wc.sample_id, project_id=wc.project_id)
    return args

# Use {root_dir}/longread.yaml to set intact_bead layout and other settings that only make sense for
# long reads
longread_cmd = """
python -m spacemake.longread \
    --parallel={threads} \
    --config=longread.yaml \
    {params.args} \
"""

rule map_cDNA:
    input: lr_cDNA
    output:
        bam=lr_cDNA_bam,
        tmp=temp(directory(lr_cDNA_dir + 'tmp/'))
    params:
        index = lambda wc : get_star_index(wc)['index'],
        annotation = lambda wc: get_species_genome_annotation(wc)['annotation'],
        star_prefix = lr_cDNA_dir + 'tmp/',
    threads: 64
    shell:
        """
        mkdir -p {params.star_prefix}
        STARlong \
            --runThreadN {threads} \
            --genomeDir  {params.index} \
            --genomeLoad NoSharedMemory \
            --readFilesIn {input} \
            --readFilesType Fastx \
            --outSAMtype BAM Unsorted \
            --outSAMunmapped Within \
            --outSAMattributes All \
            --outSAMprimaryFlag AllBestScore \
            --outStd BAM_Unsorted \
            --outFilterMultimapScoreRange 2 \
            --outFilterScoreMin 0 \
            --outFilterScoreMinOverLread 0 \
            --outFilterMatchNminOverLread 0 \
            --outFilterMatchNmin 30 \
            --outFilterMismatchNmax 1000 \
            --winAnchorMultimapNmax 200 \
            --seedSearchStartLmax 12 \
            --seedPerReadNmax 100000 \
            --seedPerWindowNmax 100 \
            --alignTranscriptsPerReadNmax 100000 \
            --alignTranscriptsPerWindowNmax 10000 \
            --outFileNamePrefix {output.tmp} | \
            {dropseq_tools}/TagReadWithGeneFunction \
            I=/dev/stdin \
            O={output.bam} \
            ANNOTATIONS_FILE={params.annotation}
        """       

rule cmd_alnstats:
    input:
        rules.map_cDNA.output.bam
    output:
        oligo_csv=lr_cDNA_oligo_analysis,
    params:
        out = lambda wc: lr_cDNA_dir.format(**wc),
    shell:
        "alnstats --parse-oligos --out-csv={params.out} --out-pdf={params.out} --out-png={params.out} {input}"

rule cmd_overview:
    input:
        reports=lambda wc: LR_REPORT_STATS
    output:
        pdf=lr_overview_pdf,
        csv=lr_overview_csv,
    params:
        out_path=lambda wc: lr_overview_dir.format(**wc),
        args=""
    shell: longread_cmd + " overview --output {params.out_path} {input.reports} "

rule cmd_report:
    input:
        stats=lr_stats
    output:
        donuts=lr_report,
        repstats=lr_report_stats
    params:
        args=get_args
    threads: 1
    shell: longread_cmd + " report"

rule cmd_extract:
    input: 
        fname = lambda wc: LR_RAW_FILES[wc.sample_id],
        ann = lr_ann
    output: lr_cDNA
    params:
        args=get_args
    log: lr_cDNA_log
    # params:
    #     known_barcodes = lambda wc: known_barcodes.get(wc.name,"")
    shell: longread_cmd + " extract {input.fname} 2> {log} > {output}"

rule cmd_edits:
    input: 
        fname = lambda wc: LR_RAW_FILES[wc.sample_id],
        stats = lr_stats
    output: lr_edits
    params:
        args=get_args
    threads: 1
    shell: longread_cmd + " edits {input.fname}"

rule cmd_annotate:
    input:
        fname = lambda wc: LR_RAW_FILES[wc.sample_id],
        ann = lr_ann
    output: lr_stats
    params:
        args=get_args
    threads: 1
    shell: longread_cmd + " annotate  {input.fname}"

rule cmd_align:
    input: 
        fname = lambda wc: LR_RAW_FILES[wc.sample_id]
    output: lr_ann
    params:
        args=get_args
    threads: 64
    shell: longread_cmd + " align {input.fname}"
