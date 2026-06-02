#########
# about #
#########
__version__ = '0.1.0'
__author__ = ['Nikos Karaiskos', 'Tamas Ryszard Sztanka-Toth', 'Marvin Jens']
__licence__ = 'GPL'
__email__ = ['nikolaos.karaiskos@mdc-berlin.de', 'tamasryszard.sztanka-toth@mdc-berlin.de', 'marvin.jens@charite.de']


    # elif wildcards.dge_cleaned == "":
    #     return {"barcode_file": top_barcodes}
    # else:

# for "no_spatial_data" we still need the top100k barcodes (do we?)

rule filter_mm_reads:
    input:
        unpack(maybe_get_puck_file), # this depends directly on the fc_....txt.gz from the capture area now
        final_bam=get_final_bam,
    output:
        pipe(final_bam_mm_included_pipe)
    threads: 2
    shell:
        "samtools view -Sh -@ {threads} {input.final_bam} | "
        "python {repo_dir}/scripts/filter_mm_reads.py "
        "  --barcode-list {input.barcode_file} "
        "  --out-sam {output} "
        "  --sample {wildcards.sample_id}"
