#########
# about #
#########
__author__ = ['Marvin Jens']
__licence__ = 'GPL'
__email__ = ['marvin.jens@mdc-berlin.de']

import spacemake.snakemake.variables as var

rule cb_correct:
    input:
        # wedge between final_bam and dropseq.smk:filter_mm_reads
        ubam=var.ubam,
        #bam=final_bam_mm_included_pipe # this is the pipe() output of dropseq.smk:filter_mm_reads rule
        bci=var.capture_area_bci
    output:
        match=var.ubam_corrected,
        # nomatch=var.ubam_nomatch
        stats=var.ubam_correction_stats
    threads: 32
    shell:
        "python -m scbamtools.bin.cb_correct sam "
        "  --input {input.ubam} "
        "  --index {input.bci} "
        "  --bam-out {output.match} "
        "  --stats-out {output.stats}"
        "  --threads {threads} "
        "  --nomatch-out discard " #{output.nomatch}"

def get_puck_barcode_files(wc):
    df = pd.read_csv(wc_fill(var.puck_barcode_files_summary, wc))
    print(">>> getting flowcell capture area barcodes")
    print(df)
    return " ".join(df['puck_barcode_file'].tolist())

rule cb_index_relevant_tiles:
    input:
        var.puck_barcode_files_summary
    params:
        puck_barcode_files=get_puck_barcode_files
    output:
        bci=var.capture_area_bci
    threads: 4
    shell:
        "cat {params.puck_barcode_files} | "
        " python -m isal.igzip -dc | "
        " python -m scbamtools.bin.cb_correct index "
        "  --index {output.bci} "

rule cb_index_corrected_sample:
    input: barcode_readcounts
    output: var.corrected_sample_bci
    shell:
        "zcat {input} | grep -v '#' | cut -f 2 | "
        " python -m scbamtools.bin.cb_correct index "
        " --index {output} "


ruleorder: make_whitelist_for_dge > create_spatial_barcode_whitelist

# only required for DropSeqTools DigitalExpression
rule make_whitelist_for_dge:
    input:
        unpack(get_puck_file),
        bci=var.corrected_sample_bci
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
