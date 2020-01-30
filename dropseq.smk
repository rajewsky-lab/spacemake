###################################################
# Snakefile containing the dropseq pipeline rules #
###################################################
rule merge_reads:
    input:
        R1=reverse_reads_mate_1,
        R2=reverse_reads_mate_2
    params:
        tmp_dir = dropseq_tmp_dir
    output:
        pipe(dropseq_merged_reads)
    shell:
        'java -jar {picard_tools} FastqToSam F1={input.R1} F2={input.R2} SM={wildcards.sample} O={output} TMP_DIR={params.tmp_dir}'

rule tag_cells:
    input:
        rules.merge_reads.output
    output:
        pipe(dropseq_cell_tagged)
    params:
        reports_dir = dropseq_reports_dir
    shell:
        """
        mkdir -p {params.reports_dir}

        {dropseq_tools}/TagBamWithReadSequenceExtended SUMMARY={params.reports_dir}/tag_cells.report.txt \
            BASE_RANGE=1-12 \
            BASE_QUALITY=10 \
            BARCODED_READ=1 \
            DISCARD_READ=False \
            TAG_NAME=XC \
            NUM_BASES_BELOW_QUALITY=1 \
            INPUT={input} \
            OUTPUT={output} \
            COMPRESSION_LEVEL=0
        """

rule tag_umis:
    input:
        rules.tag_cells.output
    output:
        pipe(dropseq_umi_tagged)
    params:
        reports_dir = dropseq_reports_dir
    shell:
        """
        {dropseq_tools}/TagBamWithReadSequenceExtended SUMMARY={params.reports_dir}/tag_umis.report.txt \
            BASE_RANGE=13-20 \
            BASE_QUALITY=10 \
            BARCODED_READ=1 \
            DISCARD_READ=True \
            TAG_NAME=XM \
            NUM_BASES_BELOW_QUALITY=1 \
            INPUT={input} \
            OUTPUT={output} \
            COMPRESSION_LEVEL=0
        """

rule remove_xc_tag:
    input:
        rules.tag_umis.output
    output:
        pipe(dropseq_tagged_filtered)
    shell:
        """
        {dropseq_tools}/FilterBam TAG_REJECT=XQ INPUT={input} OUTPUT={output} COMPRESSION_LEVEL=0
        """

rule remove_smart_adapter:
    input:
        rules.remove_xc_tag.output
    output:
        pipe(dropseq_tagged_filtered_trimmed)
    params:
        reports_dir = dropseq_reports_dir
    shell:
        """
        {dropseq_tools}/TrimStartingSequence OUTPUT_SUMMARY={params.reports_dir}/remove_smart_adapter.report.txt \
            SEQUENCE={smart_adapter} \
            MISMATCHES=0 \
            NUM_BASES=5 \
            INPUT={input} \
            OUTPUT={output} \
            COMPRESSION_LEVEL=0
        """

rule remove_polyA:
    input:
        rules.remove_smart_adapter.output
    output:
        temporary(dropseq_tagged_filtered_trimmed_polyA)
    params:
        reports_dir = dropseq_reports_dir
    shell:
        """
        {dropseq_tools}/PolyATrimmer OUTPUT_SUMMARY={params.reports_dir}/remove_polyA.report.txt \
            MISMATCHES=0 \
            NUM_BASES=6 \
            INPUT={input} \
            OUTPUT={output} \
        """

rule sam_to_fastq:
    input:
        rules.remove_polyA.output
    output:
        pipe(dropseq_star_input)
    params:
        tmp_dir = dropseq_tmp_dir
    shell:
        """
        java -jar {picard_tools} SamToFastq INPUT={input} FASTQ={output} TMP_DIR={params.tmp_dir}
        """

def get_star_inputs(wildcards):
    # This function will return 3 things required by STAR:
    #    - annotation (.gtf file)
    #    - genome (.fa file)
    #    - index (a directory where the STAR index is)
    species = config['illumina_runs'][wildcards.run]['samples'][wildcards.sample]['species']

    return {
        'annotation': config['knowledge']['annotations'][species],
        'genome': config['knowledge']['genomes'][species],
        'index': config['knowledge']['indices'][species]['star']
    }

rule map_reads:
    input:
        unpack(get_star_inputs),
        reads=rules.sam_to_fastq.output
    output:
        reads=temporary(dropseq_mapped_reads),
        log=star_log_file
    threads: 8
    params:
        tmp_dir = dropseq_tmp_dir,
        star_prefix = dropseq_root + '/star_'
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --genomeDir  {input.index} \
            --readFilesIn {input.reads} \
            --outFileNamePrefix {params.star_prefix}
        """

rule sort_mapped_reads:
    input:
        rules.map_reads.output.reads
    output:
        dropseq_mapped_sorted_reads 
    shell:
        """
        java -jar {picard_tools} SortSam I={input} O={output} SO=queryname
        """

def get_genome(wildcards):
    species = config['illumina_runs'][wildcards.run]['samples'][wildcards.sample]['species']

    return {
        'genome': config['knowledge']['genomes'][species]
    }

rule merge_bam:
    input:
        unpack(get_genome),
        unmapped_bam=rules.remove_polyA.output,
        mapped=rules.sort_mapped_reads.output
    output:
        dropseq_merged
    params:
        tmp_dir = dropseq_tmp_dir
    shell:
        """
        java -jar {picard_tools} MergeBamAlignment REFERENCE_SEQUENCE={input.genome} \
            UNMAPPED_BAM={input.unmapped_bam} \
            ALIGNED_BAM={input.mapped} \
            OUTPUT={output} \
            INCLUDE_SECONDARY_ALIGNMENTS=False \
            PAIRED_RUN=False \
            TMP_DIR={params.tmp_dir}
        """

def get_annotation(wildcards):
    species = config['illumina_runs'][wildcards.run]['samples'][wildcards.sample]['species']

    return {
        'annotation': config['knowledge']['annotations'][species]
    }

rule tag_read_with_gene_exon:
    input:
        unpack(get_annotation), 
        reads=rules.merge_bam.output
    output:
        temporary(dropseq_gene_exon_tagged)
    shell:
        """
        {dropseq_tools}/TagReadWithGeneExonFunction I={input.reads} O={output} ANNOTATIONS_FILE={input.annotation} TAG=GE
        """

rule detect_bead_substitution_error:
    input:
        rules.tag_read_with_gene_exon.output
    output:
        reads=temporary(dropseq_bead_substitution_cleaned),
        report=substitution_error_report
    threads: 4
    params:
        reports_dir = dropseq_reports_dir
    shell:
        """
        {dropseq_tools}/DetectBeadSubstitutionErrors I={input} O={output.reads} \
            OUTPUT_REPORT={output.report} \
            OUTPUT_SUMMARY={params.reports_dir}/detect_bead_substitution_error.summary.txt \
            NUM_THREADS={threads}
        """

rule detect_bead_synthesis_errors:
    input:
        rules.detect_bead_substitution_error.output
    output:
        reads=dropseq_mapped_clean_reads,
        summary=synthesis_stats_summary
    threads: 4
    params:
        reports_dir = dropseq_reports_dir
    shell:
         """
         {dropseq_tools}/DetectBeadSynthesisErrors I={input} O={output.reads} \
            REPORT={params.reports_dir}/detect_bead_synthesis_error.report.txt \
            OUTPUT_STATS={params.reports_dir}/detect_bead_synthesis_error.stats.txt \
            SUMMARY={output.summary} \
            PRIMER_SEQUENCE={smart_adapter} \
            NUM_THREADS={threads}
         """
rule bam_tag_histogram:
    input:
        dropseq_mapped_clean_reads
    output:
        dropseq_out_readcounts
    shell:
        """
        {dropseq_tools}/BamTagHistogram \
        I= {input} \
        O= {output}\
        TAG=XC
        """

rule create_top_barcodes_file:
    input:
        rules.bam_tag_histogram.output
    output:
        dropseq_top_barcodes
    shell:
        "zcat {input} | cut -f2 | head -60000 > {output}"

rule create_dge_exon_only:
    input:
        reads=dropseq_mapped_clean_reads,
        top_barcodes=dropseq_top_barcodes
    output:
        dge_exon_only 
    params:
        dge_root = dge_root
    shell:
        """
        mkdir -p {params.dge_root}

        {dropseq_tools}/DigitalExpression \
        I= {input.reads}\
        O= {output} \
        SUMMARY= {params.dge_root}/dge_summary.txt \
        CELL_BC_FILE={input.top_barcodes}
        """

rule create_dge_intron_only:
    input:
        reads=dropseq_mapped_clean_reads,
        top_barcodes=dropseq_top_barcodes
    output:
        dge_intron_only 
    params:
        dge_root = dge_root
    shell:
        """
        mkdir -p {params.dge_root}

        {dropseq_tools}/DigitalExpression \
        I= {input.reads}\
        O= {output} \
        SUMMARY= {params.dge_root}/dge_intron_summary.txt \
        CELL_BC_FILE={input.top_barcodes} \
        LOCUS_FUNCTION_LIST=null \
        LOCUS_FUNCTION_LIST=INTRONIC
        """

rule create_dge_exon_intron:
    input:
        reads=dropseq_mapped_clean_reads,
        top_barcodes=dropseq_top_barcodes
    output:
        dge_exon_intron
    params:
        dge_root = dge_root
    shell:
        """
        mkdir -p {params.dge_root}

        {dropseq_tools}/DigitalExpression \
        I= {input.reads}\
        O= {output} \
        SUMMARY= {params.dge_root}/dge_all_summary.txt \
        CELL_BC_FILE={input.top_barcodes} \
        LOCUS_FUNCTION_LIST=INTRONIC
        """

rule create_dgeReads_exon_only:
    input:
        reads=dropseq_mapped_clean_reads,
        top_barcodes=dropseq_top_barcodes
    output:
        dgeReads_exon_only 
    params:
        dge_root = dge_root
    shell:
        """
        mkdir -p {params.dge_root}

        {dropseq_tools}/DigitalExpression \
        I= {input.reads}\
        O= {output} \
        SUMMARY= {params.dge_root}/dgeReads_summary.txt \
        CELL_BC_FILE={input.top_barcodes}
        """

rule create_dgeReads_intron_only:
    input:
        reads=dropseq_mapped_clean_reads,
        top_barcodes=dropseq_top_barcodes
    output:
        dgeReads_intron_only 
    params:
        dge_root = dge_root
    shell:
        """
        mkdir -p {params.dge_root}

        {dropseq_tools}/DigitalExpression \
        I= {input.reads}\
        O= {output} \
        SUMMARY= {params.dge_root}/dgeReads_intron_summary.txt \
        CELL_BC_FILE={input.top_barcodes} \
        LOCUS_FUNCTION_LIST=null \
        LOCUS_FUNCTION_LIST=INTRONIC
        """

rule create_dgeReads_exon_intron:
    input:
        reads=dropseq_mapped_clean_reads,
        top_barcodes=dropseq_top_barcodes
    output:
        dgeReads_exon_intron
    params:
        dge_root = dge_root
    shell:
        """
        mkdir -p {params.dge_root}

        {dropseq_tools}/DigitalExpression \
        I= {input.reads}\
        O= {output} \
        SUMMARY= {dge_root}/dgeReads_all_summary.txt \
        CELL_BC_FILE={input.top_barcodes} \
        LOCUS_FUNCTION_LIST=INTRONIC
        """
