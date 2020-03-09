merged_dir = config['root_dir'] + '/projects/merged_{merged_name}'

sample_tagged_bam = merged_dir + '/{sample}_tagged.bam'
merged_bam = merged_dir + '/merged.bam'

rule tag_final_bam:
    input:
        unpack(get_dropseq_final_bam)
    output:
        temporary(sample_tagged_bam)
    threads: 4
    shell:
        """sambamba view -t {threads} -h {input} | \
            awk -v suffix=.{wildcards.sample} 'BEGIN{{OFS=FS="\t"}} /^@/ {{print $0; next}} $12=$12 suffix{{print $0}}' | \
            sambamba view -t {threads} -S /dev/stdin -o {output} -f bam"""

rule create_merged_bam:
    input:
        unpack(get_merged_inputs)
    output:
        merged_bam
    threads: 4
    shell:
        "sambamba merge -t {threads} {output} {input}"
