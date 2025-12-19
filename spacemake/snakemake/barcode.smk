#########
# about #
#########
__author__ = ['Marvin Jens']
__licence__ = 'GPL'
__email__ = ['marvin.jens@mdc-berlin.de']

import spacemake.snakemake.variables as smv

def maybe_get_puck_count_prealigned_barcode_matches_summary(wc):
    pbc =  project_df.get_puck_barcode_ids_and_files(
            project_id=wc.project_id, sample_id=wc.sample_id
        )
    print(f"maybe_get_puck_count_prealigned_barcode_matches_summary(): pbc={pbc}")
    if len(pbc[0]) > 0:
        return smv.puck_count_prealigned_barcode_matches_summary
    else:
        # we have a no_spatial_data sample
        return 'no_spatial_data'

def get_puck_barcode_files(wc, input):
    fname = input.tile_match_summary
    res = "no_spatial_data"

    if os.path.exists(fname) and fname != "no_spatial_data":
        print(f"trying to read {fname}")
        df = pd.read_csv(fname)
        print(">>> getting flowcell capture area barcodes")
        select = df["pass_threshold"] == 1
        print(df.loc[select])
        puck_barcode_files = df.loc[select, 'puck_barcode_file'].tolist()
        if puck_barcode_files:
            res = " ".join(puck_barcode_files)
    else:
        print(f"get_puck_barcode_files(): file {fname} does not exist-> no_spatial_data")

    return res


rule place_no_spatial_data_indicator:
    output:
        no_spatial_data="no_spatial_data"
    run:
        shell(
            "touch {output.no_spatial_data}"
        )

rule cb_index_relevant_tiles:
    input:
        tile_match_summary=maybe_get_puck_count_prealigned_barcode_matches_summary
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
        ubam=smv.ubam,
        bci=smv.capture_area_bci
    output:
        match=smv.ubam_corrected,
        # nomatch=smv.ubam_nomatch
        stats=smv.ubam_correction_stats
    params:
        rel_ubam=lambda wildcards, input: os.path.basename(input.ubam)
    threads: 32
    run:
        if os.path.getsize(input.bci) == 0:
            print(f"about to link {params.rel_ubam} to {output.match} because bci is empty")
            # no spatial data -> just link input to output
            shell(
                "ln -s {params.rel_ubam} {output.match} ; "
                "touch {output.stats} "
            )
        else:
            shell(
                "python -m scbamtools.bin.cb_correct sam "
                "  --input {input.ubam} "
                "  --index {input.bci} "
                "  --bam-out {output.match} "
                "  --stats-out {output.stats}"
                "  --threads {threads} "
                "  --nomatch-out discard " #{output.nomatch}"
            )

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
