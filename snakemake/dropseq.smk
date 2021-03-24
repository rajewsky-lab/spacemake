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

dropseq_tagged = dropseq_root + '/unaligned_bc_tagged.bam'
dropseq_unassigned = dropseq_root + '/unaligned_bc_unassigned.bam'

# filter out XC tag
dropseq_tagged_filtered = dropseq_root + '/unaligned_tagged_filtered.bam'

# trim smart adapter from the reads
dropseq_tagged_trimmed = dropseq_root + '/unaligned_tagged_trimmed.bam'

# trim polyA overheang if exists
dropseq_tagged_trimmed_polyA = dropseq_root + '/unaligned_tagged_trimmed_polyA.bam'

# mapped reads
dropseq_mapped_reads_unsorted_headerless = dropseq_root + '/star_Aligned.unsorted.headerless.out.bam'
dropseq_mapped_reads_unsorted = dropseq_root + '/star_Aligned.unsorted.out.bam'
dropseq_mapped_reads = dropseq_root + '/star_Aligned.sorted.out.bam'
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
        dropseq_tagged  # rules.remove_xc_tag.output
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
        temp(dropseq_tagged_trimmed_polyA)
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
        temp(dropseq_mapped_reads_unsorted_headerless)
    log: star_log_file
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
            --outSAMtype BAM Unsorted \
            --outStd BAM_Unsorted \
            --outFileNamePrefix {params.star_prefix} > {output}

        rm -rf {params.tmp_dir}
        """

rule fix_star_bam_header:
    input:
        mapped_reads=dropseq_mapped_reads_unsorted_headerless,
        unmapped_tagged_reads=dropseq_tagged_trimmed_polyA
    output: pipe(dropseq_mapped_reads_unsorted)
    shell:
        """
        python {repo_dir}/snakemake/scripts/fix_bam_header.py \
            --in-bam-star {input.mapped_reads} \
            --in-bam-tagged {input.unmapped_tagged_reads} \
            --out-bam {output}
        """

rule sort_dropseq_mapped_reads:
    input: dropseq_mapped_reads_unsorted
    output: pipe(dropseq_mapped_reads)
    threads: 4
    shell:  'sambamba sort -m 4G -o {output} -t {threads} {input}' 

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
