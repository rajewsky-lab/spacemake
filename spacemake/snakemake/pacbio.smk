#########
# about #
#########
__version__ = '0.2'
__author__ = ['Marvin Jens', 'Tamas Ryszard Sztanka-Toth']
__email__ = ['marvin.jens@mdc-berlin.de', 'tamasryszard.sztanka-toth@mdc-berlin.de']

##################
# include pacbio #
##################
pb_root = project_dir + "/processed_data/{sample}/pacbio"
pb_cache_dir = pb_root + "/cache/"
pb_ann_dir = pb_root + "/annotation/"
pb_stats_dir = pb_root + "/stats/"
pb_report_dir = pb_root + "/reports/"

# targets
pb_ann = pb_ann_dir + "{sample}.annotation.tsv"
pb_stats = pb_stats_dir + "{sample}.stats.tsv"
pb_report = pb_report_dir + "{sample}.donuts.pdf"
pb_edits = pb_report_dir + "{sample}.oligo_edits.pdf"

# pb_run_summary = processed_data_pacbio + "/{sample}.examples.txt"
# # pb_rRNA_out = processed_data_pacbio + "/{sample}.rRNA.txt"
# pb_overview = "/data/rajewsky/projects/slide_seq/.config/pb_overview.pdf"
# pb_overview_csv = "/data/rajewsky/projects/slide_seq/.config/pb_overview.csv"

# pb_bead_overview = (
#     "/data/rajewsky/projects/slide_seq/.config/pb_bead_overview.pdf"
# )

PB_RAW_FILES = {}
PB_SIGNATURE = {}
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
                out_files += \
                expand(
                    pb_report,
                    project=index[0],
                    sample=index[1],
                ) + \
                expand(
                    pb_edits,
                    project=index[0],
                    sample=index[1],
                )
                PB_RAW_FILES[index[1]] = row.longreads
                PB_SIGNATURE[index[1]] = row.longread_signature

    print("PACBIO OUTPUT FILES", out_files)
    return out_files

def get_args(wc):
    args = f""" \
    --cache={pb_cache_dir} \
    --annotation-out={pb_ann_dir} \
    --stats-out={pb_stats_dir} \
    --report-out={pb_report_dir} \
    --sample={wc.sample} \
    --signature={PB_SIGNATURE[wc.sample]} \
    """.format(sample=wc.sample, project=wc.project)
    return args

# Use {root_dir}/longread.yaml to set intact_bead layout and other settings that only make sense for
# pacbio
longread_cmd = """
python -m spacemake.longread \
    --parallel={threads} \
    --config=longread.yaml \
    {params.args} \
    {input.fname} \
"""

rule cmd_report:
    input:
        fname=lambda wc: PB_RAW_FILES[wc.sample],
        stats=pb_stats
    output:
        donuts=pb_report
    params:
        args=get_args
    threads: 1
    shell: longread_cmd + " report"

rule cmd_edits:
    input: 
        fname = lambda wc: PB_RAW_FILES[wc.sample],
        stats = pb_stats
    output: pb_edits
    params:
        args=get_args
    threads: 1
    shell: longread_cmd + " edits"

rule cmd_annotate:
    input: 
        fname = lambda wc: PB_RAW_FILES[wc.sample],
        ann = pb_ann
    output: pb_stats
    params:
        args=get_args
    threads: 1
    shell: longread_cmd + " annotate"

rule cmd_align:
    input: 
        fname = lambda wc: PB_RAW_FILES[wc.sample]
    output: pb_ann
    params:
        args=get_args
    threads: 64
    shell: longread_cmd + " align"
