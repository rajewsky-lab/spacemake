FQ, = glob_wildcards("source/{fq}.fq")
BAM, = glob_wildcards("source/lima.{bam}.bam")

reports = expand("reports/{fq}.report.pdf", fq=FQ) +  expand("reports/{bam}.report.pdf", bam=BAM)
rRNA_counts = expand("rRNA/{fq}.txt", fq=FQ) + expand("rRNA/{bam}.txt", bam=BAM)

rule all:
    input: reports, rRNA_counts
    shell:
        "python {bin} overview"

bin = "/data/rajewsky/projects/slide_seq/.repos/pb_annotate/pb_annotate"

rRNA_index = "/data/rajewsky/indices/rRNA_hsa_bwa_0.7.17/rRNA_hsa.fa"
# SAM Flag 2308 means, "not primary", "supplementary", or "unmapped". 
# We kick these out for counting rRNA hits

rule rRNA_fq:
    input: "source/{fq}.fq"
    output: "rRNA/{fq}.txt"
    shell: "bwa mem -x pacbio {rRNA_index} {input} | samtools view -F 2308 /dev/stdin | wc -l > {output}"

rule rRNA_bam:
    input: "source/lima.{bam}.bam"
    output: "rRNA/{bam}.txt"
    shell: "samtools fastq {input} | bwa mem -x pacbio {rRNA_index} - | samtools view -F 2308 /dev/stdin | wc -l > {output}"


rule annotate_fq:
    input: "source/{fq}.fq"
    output: "reports/{fq}.report.pdf"
    shell: "python {bin} scan --deletions --report {output} --summary stats/{wildcards.fq}.summary.tsv {input}  > {wildcards.fq}.txt"

rule annotate_bam:
    input: "source/lima.{bam}.bam"
    output: "reports/{bam}.report.pdf"
    shell: "python {bin} scan --deletions --report {output} --summary stats/{wildcards.bam}.summary.tsv {input}  > {wildcards.bam}.txt"
