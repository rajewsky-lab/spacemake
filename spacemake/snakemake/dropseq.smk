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

rule create_star_index:
    input:
        unpack(get_species_genome_annotation)
    output:
        directory(star_index)
    threads: max(workflow.cores * 0.25, 8)
    shell:
        """
        mkdir -p {output} 

        STAR --runMode genomeGenerate \
             --runThreadN {threads} \
             --genomeDir {output} \
             --genomeFastaFiles {input.genome} \
             --sjdbGTFfile {input.annotation}
        """

rule map_reads_final_bam:
    input:
        unpack(get_species_genome_annotation),
        unpack(get_star_input_bam),
        unpack(get_star_index),
        tagged_bam=tagged_bam
    output:
        star_log_file,
        final_bam=final_bam
    threads: 8
    params:
        tmp_dir = star_tmp_dir,
        star_prefix = star_prefix
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
            --limitOutSJcollapsed 5000000 \
            --runThreadN {threads} | \
            python {repo_dir}/scripts/fix_bam_header.py \
                --in-bam-star /dev/stdin \
                --in-bam-tagged {input.tagged_bam} \
                --out-bam /dev/stdout | \
            {dropseq_tools}/TagReadWithGeneFunction \
                I=/dev/stdin \
                O={output.final_bam} \
                ANNOTATIONS_FILE={input.annotation}

        rm -rf {params.tmp_dir}
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
