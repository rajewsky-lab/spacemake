import os
from collections import defaultdict
from spacemake.errors import *
from spacemake.util import dotdict, wc_fill
from spacemake.snakemake.variables import (
    complete_data_root,
    log_dir,
    star_log_file,
    star_target_log_file,
)

map_data = {
    # used to fetch all info needed to create a BAM file
    "MAP_RULES_LKUP": {},
    #   key: path of output BAM file
    #   value: dotdict with plenty of attributes
    # these support lookup by (project_id, sample_id)
    "ANNOTATED_BAMS": defaultdict(set),
    "REF_NAMES": defaultdict(set),
    "ALL_BAMS": defaultdict(set),
    "BAM_IS_NOT_TEMP": set(),
    "SAMPLE_MAP_STRATEGY": {},
    # used for symlink name to source mapping
    "BAM_SYMLINKS": {},
    "BAM_UNMAPPED_KEEP": set(),
    # used for automated mapping index generation
    "INDEX_FASTA_LKUP": {},
    #   key: bt2_index_file or star_index_file
    #   value: dotdict with
    #       .ref_path: path to genome FASTA
    #       .ann_path: path to annotation file (GTF) [optional]
    # needed for later stages of SPACEMAKE which require one "star.Log.final.out" file.
    "STAR_FINAL_LOG_SYMLINKS": {},
    # We create a symlink from the mapping that was accompanied by the 'final' identifier
}

#####################################
#### snakemake string templates #####
#####################################

ubam_input = "unaligned_bc_tagged.polyA_adapter_trimmed"
# this must be local and not have .bam appended!
# basically, if ubam_input were used as {target} in linked_bam it should eval to
# "unaligned_bc_tagged.polyA_adapter_trimmed.bam"

# The end-point. Usually annotation-tagged genome alignments
final_target = "final.polyA_adapter_trimmed"

# patterns for auto-generated BAM file names and symlinks
linked_bam = complete_data_root + "/{link_name}.bam"
mapped_bam = complete_data_root + "/{ref_name}.{mapper}.bam"
unmapped_bam = complete_data_root + "/not_{ref_name}.{mapper}.bam"
star_mapped_bam = complete_data_root + "/{ref_name}.STAR.bam"
star_unmapped_bam = complete_data_root + "/not_{ref_name}.STAR.bam"
bt2_mapped_bam = complete_data_root + "/{ref_name}.bowtie2.bam"
bt2_unmapped_bam = complete_data_root + "/not_{ref_name}.bowtie2.bam"

bt2_mapped_log = log_dir + "/{ref_name}.bowtie2.log"

# special log file used for rRNA "ribo depletion" stats
bt2_rRNA_log = log_dir + "/rRNA.bowtie2.bam.log"

# default places for mapping indices, unless specified differently in the config.yaml
star_index = "species_data/{species}/{ref_name}/star_index"
star_index_log = star_index + ".log"

star_index_param = star_index
star_index_file = star_index + "/SAindex"

bt2_index = "species_data/{species}/{ref_name}/bt2_index"
bt2_index_param = bt2_index + "/{ref_name}"
bt2_index_file = bt2_index_param + ".1.bt2"
bt2_index_log = bt2_index_param + ".log"
species_reference_sequence = "species_data/{species}/{ref_name}/sequence.fa"
star_idx_service = "{species}.{ref_name}.STAR_index_loaded"
species_reference_annotation = "species_data/{species}/{ref_name}/annotation.gtf"
species_reference_annotation_compiled = (
    "species_data/{species}/{ref_name}/compiled_annotation"
)
species_reference_annotation_compiled_target = (
    "species_data/{species}/{ref_name}/compiled_annotation/non_overlapping.csv"
)

#########################
### Default settings  ###
#########################

default_BT2_MAP_FLAGS = (
    " --local"
    " -L 10 -D 30 -R 30"
    " --ignore-quals"
    " --score-min=L,0,1.5"  # require 75% of perfect match (2=base match)
    " -k 10"
)
# original rRNA mapping code used --very-fast-local and that was that.

default_STAR_MAP_FLAGS = (
    # " --genomeLoad NoSharedMemory"
    " --outSAMprimaryFlag AllBestScore"
    " --outSAMattributes Standard"
    " --outSAMunmapped Within"
    " --outStd BAM_Unsorted"
    " --outSAMtype BAM Unsorted"
    " --limitOutSJcollapsed 5000000"
)

default_counting_flavor_with_annotation = "default"
default_counting_flavor_no_annotation = "custom_index"

#############################################################
##### Utility functions specific for the mapping module #####
#############################################################


def maybe_temporary(bam):
    return bam
    # DOES not currently work, as decision is made at parse time, before
    # any map-rules have been evaluated. Needs more thought.

    # # print("these are the NOT TEMP bam files", sorted(BAM_IS_NOT_TEMP))
    # if bam in map_data['BAM_IS_NOT_TEMP']:
    #     # print(bam, "is not temp!")
    #     return bam
    # else:
    #     # print(bam, "is TEMP!!!1!!11einself")
    #     return temp(bam)


def get_annotated_bams(wc):
    # print(wc)
    return {
        "annotated_bams": sorted(
            map_data["ANNOTATED_BAMS"][(wc.project_id, wc.sample_id)]
        )
    }


def get_all_mapped_bams(wc):
    # print(wc)
    return {"mapped_bams": sorted(map_data["ALL_BAMS"][(wc.project_id, wc.sample_id)])}


def get_count_flavor_str(wc):
    cflavors = []
    for bam in get_all_mapped_bams(wc)["mapped_bams"]:
        mr = map_data["MAP_RULES_LKUP"][bam]
        cflavors.append(f"{mr.ref_name}@{mr.cflavor}")

    return ",".join(cflavors)


##########################
### Core functionality ###
##########################


def validate_mapstr(mapstr, config={}, species=None):
    # print(f"validate_mapstr(mapstr={mapstr}, config={config}, species={species})")
    if species and config:
        species_d = config.get_variable("species", name=species)
    else:
        species_d = {}

    def process(token):
        """
        based on the left-hand side name (left) and the rule encoded in token,
        create one mapping rule and up to one symlink rule.
        """

        parts = token.split(":")
        if len(parts) == 2:
            mapper, ref = parts
            link_name = None
        elif len(parts) == 3:
            mapper, ref, link_name = parts
        else:
            raise ConfigVariableError(
                f"map_strategy contains a map-rule with unexpected number of parameters: {parts}"
            )

        # old format! We need to swap the order
        if ref.split("@")[0] in ["bowtie2", "STAR"]:
            ref, mapper = mapper, ref

        if species_d:
            if not ref in species_d:
                raise ConfigVariableNotFoundError(
                    ref,
                    f"reference name '{ref}' is not among the refs registered for {species}: {sorted(species_d.keys())}. Entire map-rule: {mapstr}",
                )

        if "@" in mapper:
            # we have a counting-flavor directive
            mapper, cflavor = mapper.split("@")
            if config:
                config.assert_variable("quant", cflavor)
        else:
            cflavor = None

        if cflavor:
            cfstr = f"@{cflavor}"
        else:
            cfstr = ""
        if link_name:
            lnkstr = f":{link_name}"
        else:
            lnkstr = ""

        return f"{mapper}{cfstr}:{ref}{lnkstr}"

    chain_rebuild = []
    for tokens in mapstr.split("->"):
        rebuilds = [process(r) for r in tokens.split(",")]
        chain_rebuild.append(",".join(rebuilds))

    return "->".join(chain_rebuild)


def mapstr_to_targets(mapstr, left="uBAM", final="final"):
    """
    Converts a mapping strategy provided as a string into a series of map rules, which translate into
    BAM names and their dependencies. Downstream, rule matching is guided by the convention of the
    strategy-defined BAM filenames ending in "STAR.bam" or "bowtie2.bam".

    Examples:

        rRNA:bowtie2->genome:STAR:final

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


        (3) final.bam

            a symlink pointing to the actual BAM created in (2).

    Note that one BAM must be designated 'final.bam', or the last BAM file created will be selected as final.
    (used as input to downstream processing rules for DGE creation, etc.)

    NOTE: Parallel mappings can be implemented by using commata:

        bowtie2:rRNA,STAR:genome:final

        This rule differs from the first example because it will align the unmapped reads from the uBAM
        in parallel to the rRNA reference and to the genome. In this way the same reads can match to both
        indices.

    NOTE: Gene tagging will be applied automatically if annotation data were provided for the associated
    reference index (by using spacemake config add_species --annotation=... )
    """
    mapstr = validate_mapstr(mapstr)

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
            raise ValueError(
                f"map_strategy contains a map-rule with unexpected number of parameters: {parts}"
            )

        if "@" in mapper:
            # we have a counting-flavor directive
            mapper, cflavor = mapper.split("@")
        else:
            cflavor = "auto"

        mr.input_name = left
        mr.mapper = mapper
        mr.cflavor = cflavor
        mr.ref_name = ref
        mr.out_name = f"{ref}.{mapper}"
        mr.keep_unmapped = False

        if link_name:
            lr = dotdict(
                link_src=mr.out_name, link_name=link_name, ref_name=mr.ref_name
            )

        return mr, lr

    map_rules = []
    link_rules = []

    chain = mapstr.split("->")
    final_link = None

    while chain:
        right = chain.pop(0)
        if left == right:
            continue

        last_mr = []
        for r in right.split(","):
            mr, lr = process(r, left=left)
            map_rules.append(mr)
            last_mr.append(mr)

            if lr:
                link_rules.append(lr)
                # check if we have a "final" mapping
                if final in lr.link_name:
                    final_link = lr

        # left = f"not_{mr.out_name}"
        left = mr.out_name

    for mr in last_mr:
        mr.keep_unmapped = True

    if not final_link:
        # we need to manufacture a "final" link_rule by taking the last mapping
        last = map_rules[-1]
        link_rules.append(
            dotdict(link_src=last.out_name, link_name=final, ref_name=last.ref_name)
        )

    # print("generated the following map_rules")
    # for m in map_rules:
    #     print(m)
    return map_rules, link_rules


def get_mapped_BAM_output(
    project_df=None, config=None, default_strategy="genome:STAR:final"
):
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
        map_data["SAMPLE_MAP_STRATEGY"][index] = map_strategy

        map_rules, link_rules = mapstr_to_targets(
            map_strategy, left=ubam_input, final=final_target
        )

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
            mr.out_unmapped_path = wc_fill(unmapped_bam, mr)
            if mr.keep_unmapped:
                map_data["BAM_IS_NOT_TEMP"].add(mr.out_unmapped_path)

            mr.link_name = mr.input_name
            mr.input_path = wc_fill(linked_bam, mr)

            mr.ref_path = species_d[mr.ref_name]["sequence"]
            mr.ann_path = species_d[mr.ref_name].get("annotation", None)
            map_data["ALL_BAMS"][(index[0], index[1])].add(mr.out_path)

            if mr.cflavor == "auto":
                # if we have annotation use the actual default, which works for complex
                # gene models
                # if we do NOT have annotation, assume custom_index (rRNA, spikes, CITE-seq tags etc.)
                mr.cflavor = (
                    default_counting_flavor_with_annotation
                    if mr.ann_path
                    else default_counting_flavor_no_annotation
                )

            if mr.ann_path:
                mr.ann_final = wc_fill(species_reference_annotation, mr)
                mr.ann_final_compiled = wc_fill(
                    species_reference_annotation_compiled, mr
                )
                mr.ann_final_compiled_target = wc_fill(
                    species_reference_annotation_compiled_target, mr
                )

                # keep track of all annotated BAM files we are going to create
                # for subsequent counting into DGE matrices/h5ad
                map_data["ANNOTATED_BAMS"][(index[0], index[1])].add(mr.out_path)
                map_data["REF_NAMES"][(index[0], index[1])].add(mr.ref_name)
            else:
                mr.ann_final = []

            default_STAR_INDEX = wc_fill(star_index, mr)
            default_BT2_INDEX = wc_fill(bt2_index_param, mr)
            if mr.mapper == "bowtie2":
                mr.map_index_param = species_d[mr.ref_name].get(
                    "BT2_index", default_BT2_INDEX
                )  # the parameter passed on to the mapper
                mr.map_index = os.path.dirname(mr.map_index_param)  # the index_dir
                mr.map_index_file = (
                    mr.map_index_param + ".1.bt2"
                )  # file present if the index is actually there
                mr.map_flags = species_d[mr.ref_name].get(
                    "BT2_flags", default_BT2_MAP_FLAGS
                )

            elif mr.mapper == "STAR":
                mr.map_index = species_d[mr.ref_name].get(
                    "index_dir", default_STAR_INDEX
                )
                mr.map_index_param = mr.map_index
                mr.map_index_file = mr.map_index + "/SAindex"
                mr.star_idx_service = star_idx_service.format(**mr)
                mr.map_flags = species_d[mr.ref_name].get(
                    "STAR_flags", default_STAR_MAP_FLAGS
                )

            map_data["MAP_RULES_LKUP"][mr.out_path] = mr
            map_data["INDEX_FASTA_LKUP"][mr.map_index_file] = mr
            # out_files.append(mr.out_path)

        # by convention we keep the last unmapped BAM file in a chain,
        # the others are temporary
        map_data["BAM_UNMAPPED_KEEP"].add(mr.out_unmapped_path)

        # process all symlink rules
        for lr in link_rules:
            lr.link_path = linked_bam.format(
                project_id=index[0], sample_id=index[1], link_name=lr.link_name
            )
            lr.src_path = linked_bam.format(
                project_id=index[0], sample_id=index[1], link_name=lr.link_src
            )
            map_data["BAM_SYMLINKS"][lr.link_path] = lr.src_path

            if lr.link_name == final_target:
                final_log_name = star_log_file.format(
                    project_id=index[0], sample_id=index[1]
                )
                final_log = star_target_log_file.format(
                    ref_name=lr.ref_name, project_id=index[0], sample_id=index[1]
                )
                # print("STAR_FINAL_LOG_SYMLINKS preparation", final_target, final_log_name, "->", final_log)
                map_data["STAR_FINAL_LOG_SYMLINKS"][final_log_name] = final_log

                out_files.append(lr.link_path)

    # for k,v in sorted(map_data['MAP_RULES_LKUP'].items()):
    #     print(f"map_rules for '{k}'")
    #     print(v)

    # print("BAM_SYMLINKS")
    # for k, v in map_data['BAM_SYMLINKS'].items():
    #     print(f"    output={k} <- source={v}")

    # print("out_files", out_files)
    return out_files
