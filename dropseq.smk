#########
# about #
#########
__version__ = '0.1.0'
__author__ = ['Nikos Karaiskos', 'Tamas Ryszard Sztanka-Toth']
__licence__ = 'GPL'
__email__ = ['nikolaos.karaiskos@mdc-berlin.de', 'tamasryszard.sztanka-toth@mdc-berlin.de']

########################
# COMMON PIPELINE VARS #
########################

# trim smart adapter from the reads
dropseq_tagged_trimmed = dropseq_root + '/unaligned_tagged_trimmed.bam'

# trim polyA overheang if exists
dropseq_tagged_trimmed_polyA = dropseq_root + '/unaligned_tagged_trimmed_polyA.bam'

# mapped reads
dropseq_mapped_reads = dropseq_root + '/star_Aligned.sortedByCord.out.bam'
star_log_file = dropseq_root + '/star_Log.final.out'

# final dropseq bfinal dropseq bam
dropseq_final_bam = dropseq_root + '/final.bam'

# index bam file
dropseq_final_bam_ix = dropseq_final_bam + '.bai'

###################################################
# Snakefile containing the dropseq pipeline rules #
###################################################
rule remove_smart_adapter:
    input:
        dropseq_umi_tagged  # rules.remove_xc_tag.output
    output:
        pipe(dropseq_tagged_trimmed)
    params:
        reports_dir = dropseq_reports_dir
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
        dropseq_tagged_trimmed
    output:
        temporary(dropseq_tagged_trimmed_polyA)
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

rule map_reads:
    input:
        unpack(get_species_info),
        reads=dropseq_tagged_trimmed_polyA
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
            --readFilesType SAM SE \
            --readFilesCommand samtools view \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {params.star_prefix}

        rm -rf {params.tmp_dir}
        """

rule tag_read_with_gene:
    input:
        unpack(get_species_info), 
        reads=dropseq_mapped_reads
    output:
        dropseq_final_bam
    shell:
        """
        {dropseq_tools}/TagReadWithGeneFunction \
            I={input.reads} \
            O={output} \
            ANNOTATIONS_FILE={input.annotation}
        """
