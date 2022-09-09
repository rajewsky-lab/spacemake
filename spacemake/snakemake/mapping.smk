"""
This module implements the mapping-strategy feature. This gives the freedom to define 
- for each sample - how reads should be mapped in detail: which mapper to use (currently
STAR or bowtie2) against which reference should be mapped (each species can now have arbitrary
references, each with its own sequence/FASTA file and optional annotation) and annotation-tagged.

The module takes over after pre-processing is done and hands over a "final.bam" or equivalent
to all downstream steps (DGE generation, qc_sheets, ...).
"""

#####################################
#### snakemake string templates #####
#####################################
# (may eventually move to variables.py)

# The entry point: pre-processed, unmapped reads in a uBAM
ubam_input = "unaligned_bc_tagged.polyA_adapter_trimmed"
# this must be local and not have .bam appended!
# basically, if ubam_input were used as {target} in linked_bam it should eval to 
# "unaligned_bc_tagged.polyA_adapter_trimmed.bam"

# The end-point. Usually annotation-tagged genome alignments
final_target = "final.polyA_adapter_trimmed"
# ubam_input = "unaligned_bc_tagged{polyA_adapter_trimmed}"
# final_target = "final{polyA_adapter_trimmed}"

# patterns for auto-generated BAM file names and symlinks
linked_bam = complete_data_root + "/{link_name}.bam"
mapped_bam = complete_data_root + "/{ref_name}.{mapper}.bam"
star_mapped_bam = complete_data_root + "/{ref_name}.STAR.bam"
bt2_mapped_bam = complete_data_root + "/{ref_name}.bowtie2.bam"

# special log file used for rRNA "ribo depletion" stats
bt2_rRNA_log = complete_data_root + "/rRNA.bowtie2.bam.log"

# default places for mapping indices, unless specified differently in the config.yaml
star_index = 'species_data/{species}/{ref_name}/star_index'
star_index_param = star_index
star_index_file = star_index + '/SAindex'

bt2_index = 'species_data/{species}/{ref_name}/bt2_index'
bt2_index_param = bt2_index + '/{ref_name}'
bt2_index_file = bt2_index_param + '.1.bt2'

species_reference_sequence = 'species_data/{species}/{ref_name}/sequence.fa'
species_reference_annotation = 'species_data/{species}/{ref_name}/annotation.gtf'

# used to fetch all info needed to create a BAM file
MAP_RULES_LKUP = {}
#   key: path of output BAM file
#   value: dotdict with plenty of attributes

# used for symlink name to source mapping
BAM_SYMLINKS = {}

# used for automated mapping index generation
INDEX_FASTA_LKUP = {}
#   key: bt2_index_file or star_index_file
#   value: dotdict with 
#       .ref_path: path to genome FASTA
#       .ann_path: path to annotation file (GTF) [optional]

# needed for later stages of SPACEMAKE which require one "star.Log.final.out" file.
STAR_FINAL_LOG_SYMLINKS = {}
# We create a symlink from the mapping that was accompanied by the 'final' identifier

default_BT2_MAP_FLAGS = (
    " --local"
    " -L 10 -D 30 -R 30"
    " --ignore-quals"
    " --score-min=L,0,1.5" # require 75% of perfect match (2=base match)
)
# original rRNA mapping code used --very-fast-local and that was that.

default_STAR_MAP_FLAGS = (
    " --genomeLoad NoSharedMemory"
    " --outSAMprimaryFlag AllBestScore"
    " --outSAMattributes All"
    " --outSAMunmapped Within"
    " --outStd BAM_Unsorted"
    " --outSAMtype BAM Unsorted"
    " --limitOutSJcollapsed 5000000"
)


######################################################
##### Utility functions specific for this module #####
######################################################

def wc_fill(x, wc):
    """
    Versatile helper function that can render any string template used in the mapping module,
    either filling in from a Wildcards object, or from a dotdict.
    """
    return x.format(
        sample_id=getattr(wc, "sample_id", "NO_sample_id"),
        project_id=getattr(wc, "project_id", "NO_project_id"),
        species=getattr(wc, "species", "NO_species"),
        ref_name=getattr(wc, "ref_name", "NO_ref_name"),
        mapper=getattr(wc, "mapper", "NO_mapper"),
        link_name=getattr(wc, "link_name", "NO_link_name"),
        polyA_adapter_trimmed=getattr(wc, "polyA_adapter_trimmed", "")
    )

def mapstr_to_targets(mapstr, left="uBAM", final="final"):
    """
    Converts a mapping strategy provided as a string into a series of map rules, which translate into 
    BAM names and their dependencies. Downstream, rule matching is guided by the convention of the 
    strategy-defined BAM filenames ending in "STAR.bam" or "bowtie2.bam".

    Examples:

        bowtie2:rRNA->STAR:genome:final

    The input to the first mappings is always the uBAM - the CB and UMI tagged, pre-processed, 
    but unmapped reads.
    map-rules are composed of two or three paramters, separated by ':'.
    The first two parameters for the mapping are <mapper>:<reference>. The target BAM will have 
    the name <reference>.<mapper>.bam.
    Optionally, a triplet can be used <mapper>:<reference>:<symlink> where the presence of <symlink> 
    indicates that the BAM file should additionally be made accessible under the name <symlink>.bam, a useful 
    shorthand or common hook expected by other stages of SPACEMAKE. A special symlink name is "final"
    which is required by downstream stages of SPACEMAKE. If no "final" is specified, the last map-rule
    automatically is selected and symlinked as "final".
    Note, due to the special importance of "final", it may get modified to contain other flags of the run-mode
    that are presently essential for SPACEMAKE to have in the file name ("final.bam" may in fact be 
    "final.polyA_adapter_trimmed.bam")

    The example above is going to create
        
        (1) rRNA.bowtie2.bam
        
            using bowtie2 and the index associated with the "rRNA" reference


        (2) genome.STAR.bam 
        
            using STAR on the *unmapped* reads from BAM (1)


        (3) final..bam
        
            a symlink pointing to the actual BAM created in (2).

    Note that one BAM must be designated 'final.bam', or the last BAM file created will be selected as final.
    (used as input to downstream processing rules for DGE creation, etc.)

    NOTE: Parallel mappings can be implemented by using commata:

        bowtie2:rRNA:rRNA,STAR:genome:final

        This rule differs from the first example because it will align the unmapped reads from the uBAM
        in parallel to the rRNA reference and to the genome. In this way the same reads can match to both
        indices.

    NOTE: Gene tagging will be applied automatically if annotation data were provided for the associated 
    reference index (by using spacemake config add_species --annotation=... )
    """
    def process(token, left):
        """
        based on the left-hand side name (left) and the rule encoded in token,
        create one mapping rule and up to one symlink rule.
        """

        parts = token.split(":")
        mr = dotdict()
        lr = None

        if len(parts) == 2:
            mapper, ref = parts
            link_name = None
        elif len(parts) == 3:
            mapper, ref, link_name = parts
            link_name = link_name.replace("final", final)
        else:
            raise ValueError(f"map_strategy contains a map-rule with unexpected number of parameters: {parts}")

        mr.input_name = left
        mr.mapper = mapper
        mr.ref_name = ref
        mr.out_name = f"{ref}.{mapper}"

        if link_name:
            lr = dotdict(link_src=mr.out_name, link_name=link_name, ref_name=mr.ref_name)

        return mr, lr


    map_rules = []
    link_rules = []

    chain = mapstr.split("->")
    final_link = None

    while chain:
        # print("chain:", chain)
        right = chain.pop(0)
        # print(f"left='{left}' right='{right}'")
        if left == right:
            continue

        for r in right.split(","):
            mr, lr = process(r, left=left)
            map_rules.append(mr)
            
            if lr:
                link_rules.append(lr)
                # check if we have a "final" mapping
                if final in lr.link_name:
                    final_link = lr

        left = mr.out_name

    if not final_link:
        # we need to manufacture a "final" link_rule by taking the last mapping
        last = map_rules[-1]
        link_rules.append(dotdict(link_src=last.out_name, link_name=final, ref_name=last.ref_name))

    return map_rules, link_rules


def get_mapped_BAM_output(default_strategy="STAR:genome:final"):
    """
    This function is called from main.smk at least once 
    to determine which output BAM files need to be generated and 
    to parse the map_strategy into rules and dependencies.
    """
    out_files = []

    for index, row in project_df.df.iterrows():
        # rules so far are "local" to each sample. Here we create full paths, merging with
        # project_id/sample_id from the project_df

        map_strategy = getattr(row, "map_strategy", default_strategy)
        map_rules, link_rules = mapstr_to_targets(map_strategy, left=ubam_input, final=final_target)

        species_d = project_df.config.get_variable("species", name=row.species)
        is_merged = project_df.get_metadata(
            "is_merged", project_id=index[0], sample_id=index[1]
        )
        if is_merged:
            continue

        for mr in map_rules:
            mr.project_id = index[0]
            mr.sample_id = index[1]
            mr.species = row.species
            
            mr.out_path = wc_fill(mapped_bam, mr)
            
            mr.link_name = mr.input_name
            mr.input_path = wc_fill(linked_bam, mr)

            mr.ref_path = species_d[mr.ref_name]["sequence"]
            mr.ann_path = species_d[mr.ref_name].get("annotation", None)
            if mr.ann_path:
                mr.ann_final = wc_fill(species_reference_annotation, mr)
            else:
                mr.ann_final = []

            default_STAR_INDEX = wc_fill(star_index, mr)
            default_BT2_INDEX = wc_fill(bt2_index_param, mr)
            if mr.mapper == "bowtie2":
                mr.map_flags = species_d[mr.ref_name].get("BT2_flags", default_BT2_MAP_FLAGS)
                mr.map_index_param = species_d[mr.ref_name].get("BT2_index", default_BT2_INDEX) # the parameter passed on to the mapper
                mr.map_index = os.path.dirname(mr.map_index_param) # the index_dir
                mr.map_index_file = mr.map_index_param + ".1.bt2" # file present if the index is actually there

            elif mr.mapper == "STAR":
                mr.map_flags = species_d[mr.ref_name].get("STAR_flags", default_STAR_MAP_FLAGS)
                mr.map_index = species_d[mr.ref_name].get("index_dir", default_STAR_INDEX)
                mr.map_index_param = mr.map_index
                mr.map_index_file = mr.map_index + "/SAindex"

            MAP_RULES_LKUP[mr.out_path] = mr
            INDEX_FASTA_LKUP[mr.map_index_file] = mr
            #out_files.append(mr.out_path)

        # process all symlink rules
        for lr in link_rules:
            lr.link_path = linked_bam.format(project_id=index[0], sample_id=index[1], link_name=lr.link_name)
            lr.src_path = linked_bam.format(project_id=index[0], sample_id=index[1], link_name=lr.link_src)
            BAM_SYMLINKS[lr.link_path] = lr.src_path

            if lr.link_name == final_target:
                final_log_name = star_log_file.format(project_id=index[0], sample_id=index[1])
                final_log = star_target_log_file.format(ref_name=lr.ref_name, project_id=index[0], sample_id=index[1])
                # print("STAR_FINAL_LOG_SYMLINKS preparation", target, src, final_log_name, "->", final_log)
                STAR_FINAL_LOG_SYMLINKS[final_log_name] = final_log
                
                out_files.append(lr.link_path)

    # for k,v in sorted(MAP_RULES_LKUP.items()):
    #     print(f"map_rules for '{k}'")
    #     print(v)

    # print("BAM_SYMLINKS")
    # for k, v in BAM_SYMLINKS.items():
    #     print(f"    output={k} <- source={v}")

    # print("out_files", out_files)
    return out_files

register_module_output_hook(get_mapped_BAM_output, "mapping.smk")

#############################################
#### utility functions used by the rules ####
#############################################

def get_species_reference_info(species, ref):
    refs_d = project_df.config.get_variable("species", name=species)
    return refs_d[ref]

def get_reference_attr_from_wc(wc, attr):
    d = get_species_reference_info(wc.species, wc.ref_name)
    return d[attr]

def get_map_rule(wc):
    output = wc_fill(mapped_bam, wc)
    # print(f"get_map_rules output={output}")
    return MAP_RULES_LKUP[output]

def get_map_inputs(wc, mapper="STAR"):
    wc = dotdict(wc.items())
    wc.mapper = mapper
    mr = get_map_rule(wc)
    d = {
        'bam' : mr.input_path,
        'index_file' : mr.map_index_file
    }
    if hasattr(mr, "ann_final"):
        d['annotation'] = mr.ann_final

    return d

def get_map_params(wc, output, mapper="STAR"):
    wc = dotdict(wc.items())
    wc.mapper = mapper
    mr = get_map_rule(wc)
    annotation_cmd = f"| samtools view --threads=4 -bh /dev/stdin > {output}"
    # this is a stub for "no annotation tagging"
    if hasattr(mr, "ann_final"):
        ann = mr.ann_final
        if ann and ann.lower().endswith(".gtf"):
            tagging_cmd =  "| {dropseq_tools}/TagReadWithGeneFunction I=/dev/stdin O={mr.out_path} ANNOTATIONS_FILE={mr.ann_final}"
            annotation_cmd = tagging_cmd.format(dropseq_tools=dropseq_tools, mr=mr)

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
        bam=bt2_mapped_bam
    log: bt2_mapped_bam + ".log"
    params:
        auto = lambda wc, output: get_map_params(wc, output, mapper='bowtie2'),
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
        "  {params.auto[flags]} 2> {log}"
        " "
        # fix the BAM header to accurately reflect the entire history of processing via PG records.
        "| python {repo_dir}/scripts/splice_bam_header.py"
        "  --in-ubam {input.bam}"
        " "
        "{params.auto[annotation_cmd]}"
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
        # bam=lambda wc: BAM_DEP_LKUP.get(wc_fill(star_mapped_bam, wc), f"can't_find_bam_{wc}"),
        # index=lambda wc: BAM_IDX_LKUP.get(wc_fill(star_mapped_bam, wc), f"can't find_idx_{wc}"),
        unpack(get_map_inputs)
        # bam=lambda wc: BAM_DEP_LKUP.get(wc_fill(star_mapped_bam, wc), f"can't_find_bam_{wc}"),
        # index=lambda wc: BAM_IDX_LKUP.get(wc_fill(star_mapped_bam, wc), f"can't find_idx_{wc}"),
    output:
        bam=star_mapped_bam,
        log=star_target_log_file,
    threads: 16 # bottleneck is annotation! We could push to 32 on murphy
    params:
        auto=get_map_params,
        # annotation_cmd=lambda wildcards, output: get_annotation_command(output.bam),
        # annotation=lambda wilcards, output: BAM_ANN_LKUP.get(output.bam, "can't_find_annotation"),
        # flags=lambda wildcards, output: BAM_MAP_FLAGS_LKUP.get(output.bam, "can't_find_flags"),
        tmp_dir = star_tmp_dir,
        star_prefix = star_prefix
    shell:
        "STAR {params.auto[flags]}"
        " --genomeDir {params.auto[index]}"
        " --readFilesIn {input.bam}"
        " --readFilesCommand samtools view -f 4"
        " --readFilesType SAM SE"
        " --sjdbGTFfile {params.auto[annotation]}"
        " --outFileNamePrefix {params.star_prefix}"
        " --runThreadN {threads}"
        " "
        "| python {repo_dir}/scripts/splice_bam_header.py"
        " --in-ubam {input.bam}"
        " "
        "{params.auto[annotation_cmd]}"
        " "
        "; rm -rf {params.tmp_dir}"


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
