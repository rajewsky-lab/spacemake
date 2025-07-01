import tempfile
import uuid

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

# used to fetch all info needed to create a BAM file
MAP_RULES_LKUP = map_data['MAP_RULES_LKUP']
#   key: path of output BAM file
#   value: dotdict with plenty of attributes

# used for symlink name to source mapping
BAM_SYMLINKS = map_data['BAM_SYMLINKS']

# used for automated mapping index generation
INDEX_FASTA_LKUP = map_data['INDEX_FASTA_LKUP']
#   key: bt2_index_file or star_index_file
#   value: dotdict with 
#       .ref_path: path to genome FASTA
#       .ann_path: path to annotation file (GTF) [optional]

# needed for later stages of SPACEMAKE which require one "star.Log.final.out" file.
STAR_FINAL_LOG_SYMLINKS = map_data['STAR_FINAL_LOG_SYMLINKS']

register_module_output_hook(get_mapped_BAM_output, "mapping.smk")
#####################################
#### snakemake string templates #####
#####################################

# The entry point: pre-processed, unmapped reads in a uBAM
ubam_input = "unaligned_bc_tagged.polyA_adapter_trimmed"
# this must be local and not have .bam appended!
# basically, if ubam_input were used as {target} in linked_bam it should eval to 
# "unaligned_bc_tagged.polyA_adapter_trimmed.bam"

# The end-point. Usually annotation-tagged genome alignments
final_target = "final.polyA_adapter_trimmed"
# ubam_input = "unaligned_bc_tagged{polyA_adapter_trimmed}"
# final_target = "final{polyA_adapter_trimmed}"

default_BT2_MAP_FLAGS = (
    " --local"
    " -L 10 -D 30 -R 30"
    " --ignore-quals"
    " --score-min=L,0,1.5" # require 75% of perfect match (2=base match)
)
# original rRNA mapping code used --very-fast-local and that was that.

default_STAR_MAP_FLAGS = (
    # before shared memory
    # " --genomeLoad NoSharedMemory"
    # with shared memory
    " --genomeLoad LoadAndKeep"
    " --limitBAMsortRAM 5000000000"
    " --outSAMprimaryFlag AllBestScore"
    " --outSAMattributes All"
    " --outSAMunmapped Within"
    " --outStd BAM_Unsorted"
    " --outSAMtype BAM Unsorted"
    " --limitOutSJcollapsed 5000000"
)

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

def get_star_flag(flag, default_strategy="STAR:genome:final"):
    out_files = []

    for index, row in project_df.df.iterrows():
        map_strategy = getattr(row, "map_strategy", default_strategy)
        map_rules, _ = mapstr_to_targets(map_strategy, left=ubam_input, final=final_target)
        is_merged = project_df.get_metadata(
            "is_merged", project_id=index[0], sample_id=index[1]
        )
        if is_merged:
            continue

        for mr in map_rules:
            if mr.mapper == "STAR":
                out_files += expand(flag, species=row.species, ref_name=mr.ref_name)

    out_files = sorted(set(out_files))
    # print(f"get_star_flag -> out_files={out_files}")
    return out_files

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
    wc = dotdict(wc.items())
    wc.mapper = mapper
    mr = get_map_rule(wc)
    d = {
        'bam' : mr.input_path,
        'index_file' : mr.map_index_file,
    }
    # if mapper == 'STAR':
    #     d['index_loaded'] = load_index_if_needed(expand(star_index_loaded, species=mr.species, ref_name=mr.ref_name))
    if hasattr(mr, "ann_final"):
        d['annotation'] = mr.ann_final

    return d

def get_map_params(wc, output, mapper="STAR"):
    wc = dotdict(wc.items())
    wc.mapper = mapper
    mr = get_map_rule(wc)
    ref = species_reference_sequence.format(species=mr.species, ref_name=mr.ref_name)
    annotation_cmd = f"samtools view --no-PG -F4 --threads=4 -Ch -T {ref} /dev/stdin > {output.bam}"
    # this is a stub for "no annotation tagging"
    if hasattr(mr, "ann_final"):
        ann = mr.ann_final
        if ann and ann.lower().endswith(".gtf"):
            annotation_cmd =  (
                "samtools view --no-PG -F4 --threads=2 -buh | "
                f"{dropseq_tools}/TagReadWithGeneFunction I=/dev/stdin O=/dev/stdout COMPRESSION_LEVEL=0 ANNOTATIONS_FILE={mr.ann_final} | "
                f"samtools view --no-PG --threads=4 -T {ref} -C /dev/stdin -o {mr.out_path}"
            )

    return {
        'annotation_cmd' : annotation_cmd,
        'annotation' : mr.ann_final,
        'index' : mr.map_index_param,
        'flags' : mr.map_flags,
    }

##############################################################################
#### Snakemake rules for mapping, symlinks, and mapping-index generation #####
##############################################################################

ruleorder:
    map_reads_bowtie2 > map_reads_STAR > symlinks

ruleorder:    
    symlink_final_log > map_reads_STAR

rule symlinks:
    input: lambda wc: BAM_SYMLINKS.get(wc_fill(linked_bam, wc),f"NO_BAM_SYMLINKS_for_{wc_fill(linked_bam, wc)}")
    output: linked_bam
    params:
        rel_input=lambda wildcards, input: os.path.basename(input[0])
    shell:
        "ln -s {params.rel_input} {output}"

rule symlink_final_log:
    input: lambda wc: STAR_FINAL_LOG_SYMLINKS[wc_fill(star_log_file, wc)]
    output: star_log_file
    params:
        rel_input=lambda wildcards, input: os.path.basename(input[0])
    shell:
        "ln -s {params.rel_input} {output}"

rule map_reads_bowtie2:
    input:
        # bam=lambda wc: BAM_DEP_LKUP[wc_fill(bt2_mapped_bam, wc)],
        # index=lambda wc: BAM_IDX_LKUP[wc_fill(bt2_mapped_bam, wc)],
        unpack(lambda wc: get_map_inputs(wc, mapper='bowtie2')),
    output:
        bam=bt2_mapped_bam,
        ubam=bt2_unmapped_bam
    log: bt2_mapped_bam + ".log"
    params:
        auto = lambda wc, output: get_map_params(wc, output, mapper='bowtie2'),
    threads: 32 
    shell:
        # 1) decompress unmapped reads from existing BAM/CRAM
        "samtools view -f 4 {input.bam} "
        # 2) re-convert SAM into uncompressed BAM.
        #     Somehow needed for bowtie2 2.4.5. If we don't do this
        #     bowtie2 just says "0 reads"
        " "
        "| samtools view --no-PG --threads=2 -Sbu " 
        " "
        # 3) align reads with bowtie2, *preserving the original BAM tags*
        "| bowtie2 -p {threads} --reorder --mm "
        "  -x {params.auto[index]} -b /dev/stdin --preserve-tags "
        "  {params.auto[flags]} 2> {log} "
        " "
        # fix the BAM header to accurately reflect the entire history of processing via PG records.
        "| python {repo_dir}/scripts/splice_bam_header.py "
        "  --in-ubam {input.bam}"
        " "
        "| tee >( {params.auto[annotation_cmd]} ) "
        "| samtools view -f 4 --threads=4 -Ch --no-PG > {output.ubam}"
       
        # "sambamba sort -t {threads} -m 8G --tmpdir=/tmp/tmp.{wildcards.name} -l 6 -o {output} /dev/stdin "


# TODO: unify these two functions and get rid of the params in parse_ribo_log rule below.
def get_ribo_log(wc):
    "used in params: which allows to make this purely optional w/o creating fake output"
    ribo_bam = bt2_mapped_bam.format(project_id=wc.project_id, sample_id=wc.sample_id, ref_name='rRNA')
    if ribo_bam in MAP_RULES_LKUP or ribo_bam in BAM_SYMLINKS:
        return bt2_rRNA_log.format(project_id=wc.project_id, sample_id=wc.sample_id)
    else:
        return "no_rRNA_index"

def get_ribo_log_input(wc):
    rrna_bam = bt2_mapped_bam.format(sample_id=wc.sample_id, project_id=wc.project_id, ref_name="rRNA")
    if rrna_bam in MAP_RULES_LKUP:
        # we plan to map against rRNA. This is the correct dependency:
        log = rrna_bam + '.log'
    
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
        unpack(get_map_inputs),
    output:
        bam=star_mapped_bam,
        ubam=star_unmapped_bam,
        log=star_target_log_file,
    threads: 16 # bottleneck is annotation! We could push to 32 on murphy
    params:
        auto=get_map_params,
        tmp_dir = star_tmp_dir,
        star_prefix = star_prefix
    shell:
        "STAR {params.auto[flags]}"
        " --genomeDir {params.auto[index]}"
        " --genomeLoad LoadAndKeep"
        " --readFilesIn {input.bam}"
        " --readFilesCommand samtools view -f 4"
        " --readFilesType SAM SE"
        # this needs to be removed for memory sharing
        # " --sjdbGTFfile {params.auto[annotation]}"
        " --outFileNamePrefix {params.star_prefix}"
        " --runThreadN {threads}"
        " "
        "| python {repo_dir}/scripts/splice_bam_header.py"
        " --in-ubam {input.bam}"
        " "
        "| tee >( {params.auto[annotation_cmd]} ) "
        "| samtools view -f 4 --threads=4 -Ch > {output.ubam}"
        " "
        "; rm -rf {params.tmp_dir}"


## Automatic index generation requires
# a) that the reference sequence has already been provided (by user or another rule)
# b) that the index is to be placed in the default location

from spacemake.map_strategy import get_all_mapped_bams
import spacemake.snakemake.variables as smv

def get_unmapped_cleanup_cmd(wc):
    key = (wc.project_id, wc.sample_id)
    bams = map_data["ALL_UNMAPPED_BAMS"]
    keep = map_data["BAM_IS_NOT_TEMP"]
    slated = sorted(bams[key] - keep)
    # print(f"bams={bams[key]}")
    # print(f"keep={keep}")
    # print(f"slated for removal: {slated}")
    if slated:
        return f"rm {' '.join(slated)} && "
    else:
        return ""

rule remove_unmapped:
    input:
        unpack(get_all_mapped_bams)
    output:
        smv.unmapped_removed_flag
    params:
        cleanup_cmd=get_unmapped_cleanup_cmd
    shell:
        "{params.cleanup_cmd}touch {output}"


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
    params:
        auto = lambda wc: INDEX_FASTA_LKUP[wc_fill(bt2_index_file, wc)]
    shell:
        """
        mkdir -p {params.auto[map_index]}
        bowtie2-build --ftabchars 12 \
                      --offrate 1 \
                      {params.auto[ref_path]} \
                      {params.auto[map_index_param]}
        """

rule create_star_index:
    input:
        sequence=species_reference_sequence,
        annotation=species_reference_annotation
    output:
        index_dir=directory(star_index),
        index_file=star_index_file
    threads: max(workflow.cores * 0.25, 8)
    shell:
        """
        mkdir -p {output.index_dir} 
        STAR --runMode genomeGenerate \
             --runThreadN {threads} \
             --genomeDir {output.index_dir} \
             --genomeFastaFiles {input.sequence} \
             --sjdbGTFfile {input.annotation}
        """


from spacemake.map_strategy import map_data
# print(map_data["STAR_INDICES"])

def get_unload_flags(wildcards, groupid):
    flags = []
    for (species, ref_name) in map_data["STAR_INDICES"].keys():
        flags.extend(expand(star_index_unloaded, groupid=groupid, species=species, ref_name=ref_name))

    return flags

rule unload_genome_flag:
    input: get_unload_flags

rule unload_STAR_indices_group:
    input: 
        bams = ancient(get_mapped_BAM_output(project_df))
        #get_unload_input 
        # index_dir=star_index, # we put last so it is accessible
    output:
        temp(touch(star_index_unloaded)),
    params:
        # f_locked=lambda wc: expand(star_index_locked, ref_name=wc.ref_name, species=wc.species),
        # f_locked_current=lambda wc: expand(star_index_locked_current, ref_name=wc.ref_name, species=wc.species),
        log_dir=lambda wc: expand(star_index_log_location, ref_name=wc.ref_name, species=wc.species),
        index_dir = lambda wc: map_data["STAR_INDICES"][(wc.species, wc.ref_name)]
    shell:
        """
        STAR --genomeLoad Remove --genomeDir {params.index_dir} --outFileNamePrefix {params.log_dir}/ || echo "Could not remove genome from shared memory for {output}"
        """
