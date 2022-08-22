#########
# about #
#########
__version__ = '0.1.0'
__author__ = ['Nikos Karaiskos', 'Tamas Ryszard Sztanka-Toth']
__licence__ = 'GPL'
__email__ = ['nikolaos.karaiskos@mdc-berlin.de', 'tamasryszard.sztanka-toth@mdc-berlin.de']

###################################################
# Snakefile containing the dropseq pipeline rules #
###################################################
rule remove_smart_adapter:
    input:
        tagged_bam
    output:
        pipe(tagged_trimmed_bam)
    params:
        reports_dir = reports_dir
    shell:
        """
        mkdir -p {params.reports_dir}

        {dropseq_tools}/TrimStartingSequence OUTPUT_SUMMARY={params.reports_dir}/remove_smart_adapter.report.txt \
            INPUT={input} \
            OUTPUT={output} \
            SEQUENCE={smart_adapter} \
            MISMATCHES=0 \
            NUM_BASES=5 \
            COMPRESSION_LEVEL=0
        """

rule remove_polyA:
    input:
        tagged_trimmed_bam
    output:
        temp(tagged_polyA_adapter_trimmed_bam)
    params:
        reports_dir = reports_dir
    shell:
        """
        {dropseq_tools}/PolyATrimmer OUTPUT_SUMMARY={params.reports_dir}/remove_polyA.report.txt \
            MISMATCHES=0 \
            INPUT={input} \
            OUTPUT={output} \
            NUM_BASES=6
        """

rule filter_mm_reads:
    input:
        unpack(get_final_bam)
    output:
        pipe(final_bam_mm_included_pipe)
    shell:
        """
        python {repo_dir}/scripts/filter_mm_reads.py \
            --in-bam {input} \
            --out-bam {output}
        """
