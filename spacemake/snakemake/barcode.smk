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
    # print(f"maybe_get_puck_count_prealigned_barcode_matches_summary(): pbc={pbc}")
    if len(pbc[0]) > 0:
        return smv.puck_count_prealigned_barcode_matches_summary
    else:
        # we have a no_spatial_data sample
        return 'no_spatial_data'

def get_puck_barcode_files(wc, input):
    fname = input.tile_match_summary
    res = "no_spatial_data"

    if os.path.exists(fname) and fname != "no_spatial_data":
        # print(f"trying to read {fname}")
        df = pd.read_csv(fname)
        # print(">>> getting flowcell capture area barcodes")
        select = df["pass_threshold"] == 1
        # print(df.loc[select])
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
                " python -m scbamtools.bin.cb_correct "
                "  --sample {wildcards.sample_id} "
                "  index "
                "  --index {output.bci} "
            )

rule cb_correct:
    input:
        ubam=smv.ubam,
        bci=smv.capture_area_bci
    output:
        match=smv.ubam_corrected,
        nomatch=smv.ubam_nomatch,
        stats=smv.ubam_correction_stats
    params:
        rel_ubam=lambda wildcards, input: os.path.basename(input.ubam)
    threads: 32
    run:
        pbc = project_df.get_puck_barcode_ids_and_files(
            project_id=wildcards.project_id, sample_id=wildcards.sample_id
        )
        print(pbc)
        if (len(pbc[0]) and pbc[0] != "no_spatial_data") and (os.path.getsize(input.bci) > 0):
            print(f"about to link {params.rel_ubam} to {output.match} because we do not have reference barcodes")
            # no spatial data -> just link input to output
            shell(
                "ln -s {params.rel_ubam} {output.match} ; "
                "touch {output.stats} ; touch {output.nomatch} "
            )
        else:
            shell(
                "python -m scbamtools.bin.cb_correct "
                "  --sample {wildcards.sample_id} "
                "  sam "
                "  --input {input.ubam} "
                "  --index {input.bci} "
                "  --bam-out {output.match} "
                "  --stats-out {output.stats}"
                "  --threads {threads} "
                "  --nomatch-out {output.nomatch} " #{output.nomatch}"
            )
        else:
            shell(
                "echo ln -s {params.rel_ubam} {output.match}; \n"
                "ln -s {params.rel_ubam} {output.match}; "
                "touch {output.stats}"
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

##################################
# Estimate correction gains rule #
##################################

rule cb_correct_sample:
    input:
        ubam=smv.ubam, 
        bci=smv.capture_area_bci
    output:
        stats=smv.ubam_correction_sample_stats 
    params:
        sample_size=int(config.get('ecg_sample_size', 10) * 1e6)
    threads: 32
    run:
        if os.path.getsize(input.bci) == 0:
            # no spatial data 
            shell(
                "touch {output.stats}"
            )
        else:
            shell(
                "samtools view -@6 -h {input.ubam} | head -n {params.sample_size} | "
                "python -m scbamtools.bin.cb_correct "
                "  --sample {wildcards.sample_id} "
                "  sam "
                "  --input /dev/stdin "
                "  --index {input.bci} "
                "  --bam-out /dev/null "
                "  --bam-out-mode S "
                "  --stats-out {output.stats}"
                "  --threads {threads} "
                "  --nomatch-out discard " #{output.nomatch}"
            )


rule estimate_correction_gains:
    input:
        get_output_files(ubam_correction_sample_stats,
            data_root_type = 'complete_data',
            downsampling_percentage = '',
            run_on_external=False,
            projects=config.get("projects", []),
            samples=config.get("samples", []),
            filter_merged=True
        )
    output:
        ecg="estimated_correction_gains.csv"
    run:
        # collect the estimated correction gains from the generated files
        from spacemake.snakemake.variables import ubam_correction_sample_stats

        import pandas as pd
        import numpy as np
        from scbamtools.tk import summarize_edit_stats

        # which part of the fname is project_id and sample_id?
        p_ix = ubam_correction_sample_stats.split('/').index("{project_id}")
        s_ix = ubam_correction_sample_stats.split('/').index("{sample_id}")

        # process all input files
        data = {
            'project_id': [],
            'sample_id': [],
            'estimated_correction_gain': []
        }
        for fname in input:
            data['project_id'].append(fname.split("/")[p_ix])
            data['sample_id'].append(fname.split("/")[s_ix])
            
            
            df = pd.read_csv(fname, sep='\t')
            op, (S_freq, I_freq, D_freq) = summarize_edit_stats(df)

            f = df.groupby("op")["n"].agg("sum")
            F = f / f.sum()
            all_edits = ["S", "I", "_"]
            found_edits = [e for e in all_edits if e in F.index] # intersection while preserving order
            F.loc["combined"] = F.loc[found_edits].sum()
            boost = 100 * F / F.loc["="]
            # logger.info(f"estimated correction gains for {fname}: {boost.loc['combined']:.2f} %")
            data['estimated_correction_gain'].append(np.round(boost.loc["combined"], 2))

        pd.DataFrame(data).set_index(['project_id', 'sample_id']).to_csv(output.ecg, sep='\t')
