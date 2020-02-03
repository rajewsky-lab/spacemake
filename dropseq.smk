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
        """
        java -jar {picard_tools} FastqToSam i\
            F1={input.R1} \
            F2={input.R2} \
            SM={wildcards.sample} \
            O={output} \
            TMP_DIR={params.tmp_dir} \
            SO=queryname

        rm -rf {params.tmp_dir}
        """

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

        {dropseq_tools}/TagBamWithReadSequenceExtended SUMMARY={params.reports_dir}/tag_cells.summary.txt \
            BASE_RANGE=1-12 \
            BASE_QUALITY=10 \
            BARCODED_READ=1 \
            DISCARD_READ=false \
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
        {dropseq_tools}/TagBamWithReadSequenceExtended SUMMARY={params.reports_dir}/tag_umis.summary.txt \
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
        {dropseq_tools}/FilterBam \
            TAG_REJECT=XQ \
            INPUT={input} \
            OUTPUT={output} \
            COMPRESSION_LEVEL=0
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
            INPUT={input} \
            OUTPUT={output} \
            SEQUENCE={smart_adapter} \
            MISMATCHES=0 \
            NUM_BASES=5 \
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
        java -jar {picard_tools} SamToFastq \
            INPUT={input} \
            FASTQ={output} \
            TMP_DIR={params.tmp_dir}

        rm -rf {params.tmp_dir}
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

        rm -rf {params.tmp_dir}
        """

rule sort_mapped_reads:
    input:
        rules.map_reads.output.reads
    output:
        temporary(dropseq_mapped_sorted_reads)
    shell:
        """
        java -jar {picard_tools} SortSam \
            I={input} \
            O={output} \
            SO=queryname
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
        pipe(dropseq_merged)
    params:
        tmp_dir = dropseq_tmp_dir
    shell:
        """
        java -jar {picard_tools} MergeBamAlignment \
            REFERENCE_SEQUENCE={input.genome} \
            UNMAPPED_BAM={input.unmapped_bam} \
            ALIGNED_BAM={input.mapped} \
            OUTPUT={output} \
            INCLUDE_SECONDARY_ALIGNMENTS=false \
            PAIRED_RUN=false \
            TMP_DIR={params.tmp_dir}

        rm -rf {params.tmp_dir}
        """

def get_annotation(wildcards):
    species = config['illumina_runs'][wildcards.run]['samples'][wildcards.sample]['species']

    return {
        'annotation': config['knowledge']['annotations'][species]
    }

rule tag_read_with_gene:
    input:
        unpack(get_annotation), 
        reads=rules.merge_bam.output
    output:
        temporary(dropseq_gene_tagged)
    shell:
        """
        {dropseq_tools}/TagReadWithGeneFunction \
            I={input.reads} \
            O={output} \
            ANNOTATIONS_FILE={input.annotation}
        """

rule detect_bead_substitution_error:
    input:
        rules.tag_read_with_gene.output
    output:
        reads=temporary(dropseq_bead_substitution_cleaned),
        report=substitution_error_report
    threads: 8
    params:
        reports_dir = dropseq_reports_dir
    shell:
        """
        {dropseq_tools}/DetectBeadSubstitutionErrors \
            I={input} \
            O={output.reads} \
            OUTPUT_REPORT={output.report} \
            OUTPUT_SUMMARY={params.reports_dir}/detect_bead_substitution_error.summary.txt \
            NUM_THREADS={threads}
        """

rule detect_bead_synthesis_errors:
    input:
        dropseq_bead_substitution_cleaned
    output:
        reads=dropseq_final_bam,
        summary=synthesis_stats_summary
    threads: 8
    params:
        reports_dir = dropseq_reports_dir
    shell:
         """
         {dropseq_tools}/DetectBeadSynthesisErrors \
            I={input} \
            O={output.reads} \
            REPORT={params.reports_dir}/detect_bead_synthesis_error.report.txt \
            OUTPUT_STATS={params.reports_dir}/detect_bead_synthesis_error.stats.txt \
            SUMMARY={output.summary} \
            PRIMER_SEQUENCE={smart_adapter} \
            NUM_THREADS={threads}
         """
rule bam_tag_histogram:
    input:
        dropseq_final_bam
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
        "set +o pipefail; zcat {input} | cut -f2 | head -60000 > {output}"

def get_dge_extra_params(wildcards):
    dge_type = wildcards.dge_type
    print(dge_type)

    if dge_type == '_exon':
        return ''
    elif dge_type == '_intron':
        return "LOCUS_FUNCTION_LIST=null LOCUS_FUNCTION_LIST=INTRONIC"
    elif dge_type == '_all':
        return "LOCUS_FUNCTION_LIST=INTRONIC"
    if dge_type == 'Reads_exon':
        return "OUTPUT_READS_INSTEAD=true"
    elif dge_type == 'Reads_intron':
        return "OUTPUT_READS_INSTEAD=true LOCUS_FUNCTION_LIST=null LOCUS_FUNCTION_LIST=INTRONIC"
    elif dge_type == 'Reads_all':
        return "OUTPUT_READS_INSTEAD=true LOCUS_FUNCTION_LIST=INTRONIC"
        

rule create_dge:
    input:
        reads=dropseq_final_bam,
        top_barcodes=dropseq_top_barcodes
    output:
        dge=dge_out,
        dge_summary=dge_out_summary
    params:
        dge_root = dge_root,
        dge_extra_params = lambda wildcards: get_dge_extra_params(wildcards)     
    shell:
        """
        mkdir -p {params.dge_root}

        {dropseq_tools}/DigitalExpression \
        I= {input.reads}\
        O= {output.dge} \
        SUMMARY= {output.dge_summary} \
        CELL_BC_FILE={input.top_barcodes} \
        {params.dge_extra_params}
        """
