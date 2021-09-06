merged_bam = complete_data_root + '/unaligned_bc_tagged_merged.bam'
merged_raw_reads_mate_2 = raw_reads_prefix + '_merged_R2' + reads_suffix
merged_ribo_depletion_log = complete_data_root + '/ribo_depletion_log_merged.txt'

rule create_merged_bam:
    input:
        unpack(get_files_to_merge_snakemake(tagged_bam))
    output:
        merged_bam
    shell:
        "samtools merge -o {output} {input}"

rule create_merged_raw_r2:
    input:
        unpack(get_files_to_merge_snakemake(raw_reads_mate_2))
    output:
        merged_raw_reads_mate_2
    shell:
        "cat {input} > {output}"

rule create_merged_ribo_log:
    input:
        unpack(get_files_to_merge_snakemake(ribo_depletion_log))
    output:
        merged_ribo_depletion_log
    shell:
        "cat {input} > {output}"
