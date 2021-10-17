configfile: 'config.yaml'

spaceranger_out_id = 'sr_out-{sample}-{run_type}'

spaceranger_outs = [
    spaceranger_out_id + '/outs/web_summary.html'
]

raw_reads = 'data/reads/raw/{sample_id}_S{S}_L002_R{R}_001.fastq.gz'
linked_reads = 'data/reads/linked/{sample}_S{S}_L002_R{R}_001.fastq.gz'

spaceranger_script = 'spaceranger-1.2.0/spaceranger'

linked_reads_root = 'data/reads/linked/'
raw_reads_root = 'data/reads/raw/'

run_types = ['exon', 'exon_intron']

rule all:
    input:
        expand(spaceranger_outs, sample = config['samples'].keys(), run_type = run_types)

def get_raw_reads(wildcards):
    sample_id = config['samples'][wildcards.sample]['id']

    return expand(raw_reads, sample_id = sample_id, S= wildcards.S, R = wildcards.R)

rule link_raw_reads:
    input:
        unpack(get_raw_reads)
    output:
        linked_reads
    shell:
        "ln -sr {input} {output}"

def get_spaceranger_inputs(wildcards):
    S = config['samples'][wildcards.sample]['S']
    img = config['samples'][wildcards.sample]['img']
    sample_id = config['samples'][wildcards.sample]['id'] 

    return {
        'reads': expand(raw_reads, sample_id = sample_id, S=S, R=[1,2]),
        'img': img }

def get_refdata(wildcards):
    if wildcards.run_type == 'exon':
        return 'refdata-mm10-M23'
    elif wildcards.run_type == 'exon_intron':
        return 'refdata-pre-mm10-M23'

rule run_spaceranger_counts:
    input:
        unpack(get_spaceranger_inputs)
    output:
        spaceranger_outs
    params:
        area = lambda wildcards: config['samples'][wildcards.sample]['area'],
        sample_id = lambda wildcards: config['samples'][wildcards.sample]['id'],
        refdata = lambda wildcards: get_refdata(wildcards),
        run_id = spaceranger_out_id
    wildcard_constraints:
       run_type='|'.join(run_types) 
    threads: 8
    shell:
        # first we remove the directory, otherwise space ranger is gonna fail
        # the directory is created by snakemake, by default. but after creation
        # spaceranger thinks that it has already run...
        """
        rm -rf {params.run_id}
        {spaceranger_script} count --id={params.run_id} \
            --transcriptome={params.refdata} \
            --fastqs={raw_reads_root} \
            --sample={params.sample_id} \
            --image={input.img} \
            --localcores={threads} \
            --localmem=64 \
            --unknown-slide \
            --reorient-images
        """
