annotation_file = os.path.join(config['root_dir'],
    config['annotation_file_pattern'])
genome_file = os.path.join(config['root_dir'],
    config['genome_file_pattern'])

rule all:
    input:
        expand(annotation_file, species = config['species'],
               data_type = 'annotation'),
        expand(genome_file, species = config['species'],
               data_type = 'genome')

rule unzip:
    input:
        '{filename}.gz'
    output:
        '{filename}'
    shell: "unpigz {input}"

def get_url(wildcards):
    return config[wildcards.species + '_' + wildcards.data_type + '_url']

rule download_species_annotation:
    output:
        annotation_file.replace('.gtf', '.gtf.gz')
    params:
        url = lambda wildcards: get_url(wildcards)
    shell:
        "wget -O {output} {params.url}"

rule download_species_genome:
    output:
        genome_file.replace('.fa', '.fa.gz')
    params:
        url = lambda wildcards: get_url(wildcards)
    shell:
        "wget -O {output} {params.url}"
