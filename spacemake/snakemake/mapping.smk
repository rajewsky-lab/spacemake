"""
This module implements the mapping-strategy feature. This gives the freedom to define 
- for each sample - how reads should be mapped in detail: which mapper to use (currently
STAR or bowtie2) against which reference should be mapped (each species can now have arbitrary
references, each with its own sequence/FASTA file and optional annotation) and annotation-tagged.

The module takes over after pre-processing is done and hands over a "final.bam" or equivalent
to all downstream steps (DGE generation, qc_sheets, ...).
"""

# This is where all the python functions and string constants
# live
from spacemake.map_strategy import *

register_module_output_hook(get_mapped_BAM_output, "mapping.smk")

# TODO: port remaining python code to map_strategy.py
# to expose it to enable coverage analysis and unit-testing
# possibly best way is to turn map_strategy into a class-instance that can be set-up with 
# project_df and config and make all these utility functions methods.
# that would also remove the need for map_data which is already an attempt to mitigate global
# variables.

def is_paired_end(project_id, sample_id):
    row = project_df.df.loc[(project_id, sample_id)].to_dict()
    #print(f"{project_id} {sample_id} -> {row['barcode_flavor']}")
    flavor_d = project_df.config.get_variable("barcode_flavors", name=row['barcode_flavor'])
    pe = flavor_d.get('paired_end', 'single-end')

    return pe != 'single-end'

def get_species_reference_info(species, ref):
    refs_d = project_df.config.get_variable("species", name=species)
    return refs_d[ref]

def get_reference_attr_from_wc(wc, attr):
    d = get_species_reference_info(wc.species, wc.ref_name)
    return d[attr]

def get_map_rule(wc):
    output = wc_fill(mapped_bam, wc)
    # print(f"get_map_rules output={output}")
    return map_data['MAP_RULES_LKUP'][output]

def get_map_inputs(wc, mapper="STAR"):
    # print("get_map_inputs()")
    wc = dotdict(wc.items())
    wc.mapper = mapper
    mr = get_map_rule(wc)
    d = {
        'bam' : mr.input_path,
        'index_file' : mr.map_index_file
    }
    if hasattr(mr, "ann_final_compiled_target") and mr.ann_final_compiled_target:
        d['annotation'] = mr.ann_final_compiled_target

    # print(d)
    # print("done.")
    return d

def get_map_params(wc, output, mapper="STAR"):
    # we are expecting an input stream of - possibly uncompressed - BAM
    # this stream will get directed to either compressed BAM directly, or first
    # annotated and then written to compressed BAM

    # print("get_map_params()")
    wc = dotdict(wc.items())
    wc.mapper = mapper
    mr = get_map_rule(wc)
    annotation_cmd = f"| samtools view --threads=4 -bh /dev/stdin > {output.bam}"
    # this is a stub for "no annotation tagging"
    if hasattr(mr, "ann_final"):
        ann = mr.ann_final
        if ann and ann.lower().endswith(".gtf"):
            ann_log = wc_fill(log_dir, wc) + f"/{mr.ref_name}.{mapper}.annotator.log"
            stats_out = wc_fill(stats_dir, wc) + f"/{mr.ref_name}.{mapper}.annotator.tsv"
            annotation_cmd = (
                f"| python {bin_dir}/annotator.py "
                f"  --sample={wc.sample_id} "
                f"  --log-level={log_level} "
                f"  --log-file={ann_log} "
                f"  --debug={log_debug} "
                f"  tag "
                f"  --bam-in=/dev/stdin "
                f"  --bam-out={mr.out_path} "
                f"  --stats-out={stats_out} "
                f"  --compiled={mr.ann_final_compiled} "
                f"| samtools view --threads=4 -bh --no-PG > {output.bam} "
            )

    d = {
        'annotation_cmd' : annotation_cmd,
        'annotation' : mr.ann_final,
        'index' : mr.map_index_param,
        'flags' : mr.map_flags,
    }
    # print("done..")
    return d


##############################################################################
#### Snakemake rules for mapping, symlinks, and mapping-index generation #####
##############################################################################

ruleorder:
    map_reads_bowtie2 > map_reads_STAR > symlinks

ruleorder:    
    symlink_final_log > map_reads_STAR

rule symlinks:
    input: lambda wc: map_data['BAM_SYMLINKS'].get(wc_fill(linked_bam, wc),f"NO_BAM_SYMLINKS_for_{wc_fill(linked_bam, wc)}")
    output: linked_bam
    params:
        rel_input=lambda wildcards, input: os.path.basename(input[0])
    shell:
        "ln -s {params.rel_input} {output}"

rule symlink_final_log:
    input: lambda wc: map_data['STAR_FINAL_LOG_SYMLINKS'][wc_fill(star_log_file, wc)]
    output: star_log_file
    params:
        rel_input=lambda wildcards, input: os.path.basename(input[0])
    shell:
        "ln -s logs/{params.rel_input} {output}"

rule map_reads_bowtie2:
    input:
        # bam=lambda wc: BAM_DEP_LKUP[wc_fill(bt2_mapped_bam, wc)],
        # index=lambda wc: BAM_IDX_LKUP[wc_fill(bt2_mapped_bam, wc)],
        unpack(lambda wc: get_map_inputs(wc, mapper='bowtie2')),
    output:
        bam=bt2_mapped_bam,
        ubam=bt2_unmapped_bam, # maybe_temporary()
    log: 
        bt2 = bt2_mapped_log,
        hdr = bt2_mapped_log.replace('.log', '.splice_bam_header.log')
    params:
        auto = lambda wc, output: get_map_params(wc, output, mapper='bowtie2'),
        PE=lambda wc: "--align-paired-reads --no-mixed" if is_paired_end(wc.project_id, wc.sample_id) else "",
    threads: 32 
    shell:
        # 1) decompress unmapped reads from existing BAM
        "samtools view -f 4 {input.bam}"
        # 2) re-convert SAM into uncompressed BAM.
        #     Somehow needed for bowtie2 2.4.5. If we don't do this
        #     bowtie2 just says "0 reads"
        " "
        "| samtools view --no-PG --threads=2 -Sbu" 
        " "
        # 3) align reads with bowtie2, *preserving the original BAM tags*
        "| bowtie2 -p {threads} --reorder --mm"
        "  -x {params.auto[index]} -b /dev/stdin --preserve-tags"
        "  {params.auto[flags]} {params.PE} 2> {log.bt2}"
        " "
        # fix the BAM header to accurately reflect the entire history of processing via PG records.
        "| python {bin_dir}/splice_bam_header.py"
        "  --in-ubam {input.bam} "
        "  --log-level={log_level} "
        "  --log-file={log.hdr} "
        " "
        "| tee >( samtools view -F 4 --threads=2 -buh {params.auto[annotation_cmd]} ) "
        "| samtools view -f 4 --threads=4 -bh > {output.ubam}"


# TODO: unify these two functions and get rid of the params in parse_ribo_log rule below.
def get_ribo_log(wc):
    "used in params: which allows to make this purely optional w/o creating fake output"
    ribo_bam = bt2_mapped_bam.format(project_id=wc.project_id, sample_id=wc.sample_id, ref_name='rRNA')
    if ribo_bam in map_data['MAP_RULES_LKUP'] or ribo_bam in map_data['BAM_SYMLINKS']:
        log = bt2_mapped_log.format(sample_id=wc.sample_id, project_id=wc.project_id, ref_name="rRNA")
        return log
    else:
        return "no_rRNA_index"

def get_ribo_log_input(wc):
    rrna_bam = bt2_mapped_bam.format(sample_id=wc.sample_id, project_id=wc.project_id, ref_name="rRNA")
    if rrna_bam in map_data['MAP_RULES_LKUP']:
        # we plan to map against rRNA. This is the correct dependency:
        log = bt2_mapped_log.format(sample_id=wc.sample_id, project_id=wc.project_id, ref_name="rRNA")
    
    else:
        # no rRNA reference available. Default to the stub
        log = []

    return log

rule parse_ribo_log:
    input: get_ribo_log_input
    output: parsed_ribo_depletion_log
    params:
        ribo_log = get_ribo_log
    script: 'scripts/parse_ribo_log.py'


rule map_reads_STAR:
    input: 
        # bam=lambda wc: BAM_DEP_LKUP.get(wc_fill(star_mapped_bam, wc), f"can't_find_bam_{wc}"),
        # index=lambda wc: BAM_IDX_LKUP.get(wc_fill(star_mapped_bam, wc), f"can't find_idx_{wc}"),
        unpack(get_map_inputs)
        # bam=lambda wc: BAM_DEP_LKUP.get(wc_fill(star_mapped_bam, wc), f"can't_find_bam_{wc}"),
        # index=lambda wc: BAM_IDX_LKUP.get(wc_fill(star_mapped_bam, wc), f"can't find_idx_{wc}"),
    output:
        bam=star_mapped_bam,
        ubam=star_unmapped_bam, # maybe_temporary(
        tmp=temp(directory(star_prefix + "_STARgenome"))
    threads: 16 # bottleneck is annotation! We could push to 32 on murphy
    log:
        star=star_target_log_file,
        hdr=star_target_log_file.replace(".log", ".splice_bam_header.log")
    params:
        auto=get_map_params,
        PE=lambda wc: "PE" if is_paired_end(wc.project_id, wc.sample_id) else "SE",
        # annotation_cmd=lambda wildcards, output: get_annotation_command(output.bam),
        # annotation=lambda wilcards, output: BAM_ANN_LKUP.get(output.bam, "can't_find_annotation"),
        # flags=lambda wildcards, output: BAM_MAP_FLAGS_LKUP.get(output.bam, "can't_find_flags"),
        tmp_dir = star_tmp_dir,
        star_prefix = star_prefix
    shell:
        "STAR {params.auto[flags]}"
        "  --genomeDir {params.auto[index]}"
        "  --readFilesIn {input.bam}"
        "  --readFilesCommand samtools view -f 4"
        "  --readFilesType SAM {params.PE}"
        "  --sjdbGTFfile {params.auto[annotation]}"
        "  --outFileNamePrefix {params.star_prefix}"
        "  --runThreadN {threads}"
        " "
        "| python {bin_dir}/splice_bam_header.py"

        "  --in-ubam {input.bam}"
        "  --log-level={log_level} "
        "  --log-file={log.hdr} "
        "  --debug={log_debug} "
        " "
        "| tee >( samtools view -F 4 --threads=2 -buh {params.auto[annotation_cmd]} ) "
        "| samtools view -f 4 --threads=4 -bh > {output.ubam}"
        " "
        "; rm -rf {params.tmp_dir}"
        "; rm {params.star_prefix}Log.out"
        "; rm {params.star_prefix}Log.std.out"
        "; rm {params.star_prefix}Log.progress.out"
        "; mv {params.star_prefix}Log.final.out {log.star}"

## Automatic index generation requires
# a) that the reference sequence has already been provided (by user or another rule)
# b) that the index is to be placed in the default location

rule prepare_species_reference_sequence:
	input:
		lambda wc: get_reference_attr_from_wc(wc, 'sequence')
	output:
		species_reference_sequence
	run:
		if input[0].endswith('.fa.gz'):
			shell('unpigz -c {input} > {output}')
		else:
			shell('ln -sr {input} {output}')


rule prepare_species_reference_annotation:
	input:
		lambda wc: get_reference_attr_from_wc(wc, 'annotation')
	output:
		species_reference_annotation
	run:
		if input[0].endswith('.gtf.gz'):
			shell('unpigz -c {input} > {output}')
		else:
			shell('ln -sr {input} {output}')

# TODO: transition to species_reference_file and map_index_param
# and get rid of INDEX_FASTA_LKUP
rule create_bowtie2_index:
    input:
        species_reference_sequence
    output:
        bt2_index_file
    log:
        bt2_index_log
    params:
        auto = lambda wc: map_data['INDEX_FASTA_LKUP'][wc_fill(bt2_index_file, wc)]
    shell:
        "mkdir -p {params.auto[map_index]} \n"
        "bowtie2-build --ftabchars 12 "
        "  --offrate 1 "
        "  {params.auto[ref_path]} "
        "  {params.auto[map_index_param]} "
        " &> {log} "

rule create_star_index:
    input:
        sequence=species_reference_sequence,
        annotation=species_reference_annotation
    output:
        index_dir=directory(star_index),
        index_file=star_index_file
    log: star_index_log
    threads: max(workflow.cores * 0.25, 8)
    shell:
        """
        mkdir -p {output.index_dir} 
        STAR --runMode genomeGenerate \
             --runThreadN {threads} \
             --genomeDir {output.index_dir} \
             --genomeFastaFiles {input.sequence} \
             --sjdbGTFfile {input.annotation} &> {log}
        """

rule compile_annotation:
    input: species_reference_annotation
    output: 
        target = species_reference_annotation_compiled_target,
        path = directory(species_reference_annotation_compiled)
    log: species_reference_annotation_compiled + '/annotator_build.log'
    shell:
        "python {bin_dir}/annotator.py "
        "  --log-file={log} "
        "  --log-level={log_level} "
        "  --debug={log_debug} "
        "  build "
        "  --gtf={input} "
        "  --compiled={output.path}"

