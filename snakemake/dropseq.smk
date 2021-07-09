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
        tagged_bam  # rules.remove_xc_tag.output
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
        temp(tagged_polyA_trimmed_bam)
    params:
        reports_dir = reports_dir
    shell:
        """
        {dropseq_tools}/PolyATrimmer OUTPUT_SUMMARY={params.reports_dir}/remove_polyA.report.txt \
            MISMATCHES=0 \
            NUM_BASES=20 \
            INPUT={input} \
            OUTPUT={output} \
        """

rule map_reads:
    input:
        unpack(get_species_info),
        unpack(get_star_input_bam)
    output:
        temp(mapped_reads_sorted_headerless)
    log: star_log_file
    threads: 24
    params:
        tmp_dir = tmp_dir,
        star_prefix = lambda wildcards: complete_data_root + '/star' + wildcards.polyA_trimmed + '_'
    shell:
        """
        STAR \
            --genomeLoad NoSharedMemory \
            --genomeDir  {input.index} \
            --sjdbGTFfile {input.annotation} \
            --readFilesCommand samtools view \
            --readFilesIn {input.reads} \
            --readFilesType SAM SE \
            --outFileNamePrefix {params.star_prefix} \
            --outSAMprimaryFlag AllBestScore \
            --outSAMattributes All \
            --outSAMunmapped Within \
            --outStd BAM_Unsorted \
            --outSAMtype BAM Unsorted \
            --runThreadN {threads} > {output}

        rm -rf {params.tmp_dir}
        """

rule fix_star_bam_header:
    input:
        mapped_reads=mapped_reads_sorted_headerless,
        unmapped_tagged_reads=tagged_bam
    output: pipe(mapped_reads_sorted)
    shell:
        """
        python {repo_dir}/snakemake/scripts/fix_bam_header.py \
            --in-bam-star {input.mapped_reads} \
            --in-bam-tagged {input.unmapped_tagged_reads} \
            --out-bam {output}
        """

rule tag_read_with_gene:
    input:
        unpack(get_species_info), 
        reads=mapped_reads_sorted
    output:
        final_bam
    shell:
        """
        {dropseq_tools}/TagReadWithGeneFunction \
            I={input.reads} \
            O={output} \
            ANNOTATIONS_FILE={input.annotation}
        """

rule filter_mm_reads:
    input:
        final_bam
    output:
        pipe(final_bam_mm_filtered)
    shell:
        """
        python {repo_dir}/snakemake/scripts/filter_mm_reads.py \
            --in-bam {input} \
            --out-bam {output}
        """
