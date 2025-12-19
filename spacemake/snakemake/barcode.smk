#########
# about #
#########
__author__ = ['Marvin Jens']
__licence__ = 'GPL'
__email__ = ['marvin.jens@mdc-berlin.de']

import spacemake.snakemake.variables as smv

def get_puck_barcode_files(wc, input):

    df = pd.read_csv(input.tile_match_summary)
    print(">>> getting flowcell capture area barcodes")
    select = df["pass_threshold"] == 1
    print(df.loc[select])

    pbc =  project_df.get_puck_barcode_ids_and_files(
            project_id=wc.project_id, sample_id=wc.sample_id
        )
    print(f"pdf.get_puck_barcode_ids_and_files() -> {pbc}")
    res = " ".join(df.loc[select, 'puck_barcode_file'].tolist())
    if not res:
        res = "no_spatial_data"

    return res

rule cb_index_relevant_tiles:
    input:
        tile_match_summary=smv.puck_count_prealigned_barcode_matches_summary
    params:
        puck_barcode_files=get_puck_barcode_files
    output:
        bci=smv.capture_area_bci
    threads: 4
    run:
        if params.puck_barcode_files == "no_spatial_data":
            # we have no whitelist. So no real BCI is needed
            shell(
                "touch {output.bci}"
            )
        else:
            # the real McCoy.
            shell(
                "cat {params.puck_barcode_files} | "
                " python -m isal.igzip -dc | "
                " python -m scbamtools.bin.cb_correct index "
                "  --index {output.bci} "
            )

rule cb_correct:
    input:
        # wedge between final_bam and dropseq.smk:filter_mm_reads
        ubam=smv.ubam,
        #bam=final_bam_mm_included_pipe # this is the pipe() output of dropseq.smk:filter_mm_reads rule
        bci=smv.capture_area_bci
    output:
        match=smv.ubam_corrected,
        # nomatch=smv.ubam_nomatch
        stats=smv.ubam_correction_stats
    params:
        rel_ubam=lambda wildcards, input: os.path.basename(input.ubam)
    threads: 32
    shell:
        "python -m scbamtools.bin.cb_correct sam "
        "  --input {input.ubam} "
        "  --index {input.bci} "
        "  --bam-out {output.match} "
        "  --stats-out {output.stats}"
        "  --threads {threads} "
        "  --nomatch-out discard " #{output.nomatch}"

def get_correction_reference(wc, input):
    # print("get_puck_barcode_files")
    # print(project_df.get_puck_barcode_ids_and_files(
    #         wc.project_id, wc.sample_id
    #     )
    # )
    df = pd.read_csv(wc_fill(var.puck_barcode_files_summary, wc))
    # print(">>> getting flowcell capture area barcodes")
    # print(df)
    if len(df) == 0:
        # fallback: top 100k barcodes
        return f"zcat {input.bc_counts} | grep -v '#' | head -n 100000 | cut -f 2 | "
    else:
        files = " ".join(df['puck_barcode_file'].tolist())
        return f"cat {files} | python -m isal.igzip -dc | "

# rule cb_index_relevant_tiles:
#     input:
#         puck_summary=var.puck_barcode_files_summary,
#         bc_counts=barcode_readcounts_prealigned
#         # top=top_barcodes
#     params:
#         cb_ref=get_correction_reference
#     output:
#         bci=var.capture_area_bci
#     threads: 4
#     shell:
#         " {params.cb_ref} python -m scbamtools.bin.cb_correct index "
#         "  --index {output.bci} "

rule cb_index_corrected_sample:
    input: barcode_readcounts
    output: smv.corrected_sample_bci
    shell:
        "zcat {input} | grep -v '#' | cut -f 2 | "
        " python -m scbamtools.bin.cb_correct index "
        " --index {output} "


ruleorder: make_whitelist_for_dge > create_spatial_barcode_whitelist

# only required for DropSeqTools DigitalExpression
rule make_whitelist_for_dge:
    input:
        unpack(get_puck_file),
        bci=smv.corrected_sample_bci
    output:
        spatial_barcodes_corrected
    shell:
        "zcat {input.barcode_file} | "
        " python -m scbamtools.bin.cb_correct query "
        "  --unique "
        "  --index {input.bci} "
        "  --dist 0 "
        "  --out-mode match "
        "  --output {output} "
