#########
# about #
#########
__version__ = '0.2'
__author__ = ['Marvin Jens', 'Tamas Ryszard Sztanka-Toth']
__email__ = ['marvin.jens@mdc-berlin.de', 'tamasryszard.sztanka-toth@mdc-berlin.de']

##################
# include pacbio #
##################
pb_root = project_dir + "/processed_data/{sample_id}/pacbio"
pb_cache_dir = pb_root + "/cache/"
pb_ann_dir = pb_root + "/annotation/"
pb_stats_dir = pb_root + "/stats/"
pb_report_dir = pb_root + "/reports/"
pb_examples_dir = pb_root + "/examples/"
pb_cDNA_dir = pb_root + "/cDNA/"

# targets
pb_ann = pb_ann_dir + "{sample_id}.annotation.tsv"
pb_stats = pb_stats_dir + "{sample_id}.stats.tsv"
pb_report = pb_report_dir + "{sample_id}.donuts.pdf"
pb_report_stats = pb_stats_dir + "{sample_id}.report.tsv"
pb_edits = pb_report_dir + "{sample_id}.oligo_edits.pdf"
pb_cDNA = pb_cDNA_dir + "{sample_id}.fa"
pb_cDNA_log = pb_cDNA_dir + "{sample_id}.log"
pb_cDNA_oligo_analysis = pb_cDNA_dir + "{sample_id}.oligo_analysis.csv"
pb_cDNA_bam = pb_cDNA_dir + "{sample_id}.bam"
pb_examples = pb_examples_dir + "{sample_id}.txt"

# pacbio_overview = '/data/rajewsky/projects/slide_seq/.config/pacbio_overview.pdf'
# print(config)
pb_overview_dir = os.path.join(config['root_dir'], 'pacbio_overview/')
pb_overview_csv = pb_overview_dir + 'overview.csv'
pb_overview_pdf = pb_overview_dir + 'fidelity.pdf'

# pb_run_summary = processed_data_pacbio + "/{sample}.examples.txt"
# # pb_rRNA_out = processed_data_pacbio + "/{sample}.rRNA.txt"
# pb_overview = "/data/rajewsky/projects/slide_seq/.config/pb_overview.pdf"
# pb_overview_csv = "/data/rajewsky/projects/slide_seq/.config/pb_overview.csv"

# pb_bead_overview = (
#     "/data/rajewsky/projects/slide_seq/.config/pb_bead_overview.pdf"
# )

PB_RAW_FILES = {}
PB_SIGNATURE = {}
PB_REPORT_STATS = []
def get_longread_output():
    """
    This function is called from main.smk at least once 
    to determine which output files need to be generated
    from pacbio longread analysis.
    We use this opportunity to populate PB_RAW_FILES
    """
    out_files = []
    for index, row in project_df.df.iterrows():
        # for run_mode in row["run_mode"]:
        #     run_mode_variables = project_df.config.get_run_mode(run_mode).variables
            if row.longreads:
                PB_REPORT_STATS.extend(
                    expand(pb_report_stats, project_id=index[0], sample_id=index[1])
                )
                out_files += \
                expand(
                    pb_report,
                    project_id=index[0],
                    sample_id=index[1],
                ) + \
                expand(
                    pb_edits,
                    project_id=index[0],
                    sample_id=index[1],
                ) + \
                expand(
                    pb_cDNA_bam,
                    project_id=index[0],
                    sample_id=index[1],
                ) + \
                expand(
                    pb_cDNA_oligo_analysis,
                    project_id=index[0],
                    sample_id=index[1],
                )

                PB_RAW_FILES[index[1]] = row.longreads
                PB_SIGNATURE[index[1]] = row.longread_signature

    # if we have any pacbio analysis, generate an overview plot
    if out_files:
        out_files.append(pb_overview_pdf)

    # print("PACBIO OUTPUT FILES", out_files)
    # print("PB_REPORT_STATS", PB_REPORT_STATS)
    return out_files

def get_args(wc):
    args = f""" \
    --cache={pb_cache_dir} \
    --annotation-out={pb_ann_dir} \
    --stats-out={pb_stats_dir} \
    --report-out={pb_report_dir} \
    --examples-out={pb_examples_dir} \
    --sample={wc.sample_id} \
    --signature={PB_SIGNATURE[wc.sample_id]} \
    """.format(sample_id=wc.sample_id, project_id=wc.project_id)
    return args

# Use {root_dir}/longread.yaml to set intact_bead layout and other settings that only make sense for
# pacbio
longread_cmd = """
python -m spacemake.longread \
    --parallel={threads} \
    --config=longread.yaml \
    {params.args} \
"""

rule map_cDNA:
    input: pb_cDNA
    output:
        bam=pb_cDNA_bam,
        tmp=temp(directory(pb_cDNA_dir + 'tmp/'))
    params:
        index = lambda wc : get_star_index(wc)['index'],
        annotation = lambda wc: get_species_genome_annotation(wc)['annotation'],
        star_prefix = pb_cDNA_dir + 'tmp/',
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
        oligo_csv=pb_cDNA_oligo_analysis,
    params:
        out = lambda wc: pb_cDNA_dir.format(**wc),
    shell:
        "alnstats --parse-oligos --out-csv={params.out} --out-pdf={params.out} --out-png={params.out} {input}"

rule cmd_overview:
    input:
        reports=lambda wc: PB_REPORT_STATS
    output:
        pdf=pb_overview_pdf,
        csv=pb_overview_csv,
    params:
        out_path=lambda wc: pb_overview_dir.format(**wc),
        args=""
    shell: longread_cmd + " overview --output {params.out_path} {input.reports} "

rule cmd_report:
    input:
        stats=pb_stats
    output:
        donuts=pb_report,
        repstats=pb_report_stats
    params:
        args=get_args
    threads: 1
    shell: longread_cmd + " report"

rule cmd_extract:
    input: 
        fname = lambda wc: PB_RAW_FILES[wc.sample_id],
        ann = pb_ann
    output: pb_cDNA
    params:
        args=get_args
    log: pb_cDNA_log
    # params:
    #     known_barcodes = lambda wc: known_barcodes.get(wc.name,"")
    shell: longread_cmd + " extract {input.fname} 2> {log} > {output}"

rule cmd_edits:
    input: 
        fname = lambda wc: PB_RAW_FILES[wc.sample_id],
        stats = pb_stats
    output: pb_edits
    params:
        args=get_args
    threads: 1
    shell: longread_cmd + " edits {input.fname}"

rule cmd_annotate:
    input:
        fname = lambda wc: PB_RAW_FILES[wc.sample_id],
        ann = pb_ann
    output: pb_stats
    params:
        args=get_args
    threads: 1
    shell: longread_cmd + " annotate  {input.fname}"

rule cmd_align:
    input: 
        fname = lambda wc: PB_RAW_FILES[wc.sample_id]
    output: pb_ann
    params:
        args=get_args
    threads: 64
    shell: longread_cmd + " align {input.fname}"
