tagging_cmd =  "{dropseq_tools}/TagReadWithGeneFunction I=/dev/stdin O={output.out_bam} ANNOTATIONS_FILE={params.ref_ann}"
bam_base = project_dir + "/processed_data/{sample_id}/illumina/complete_data/"
ubam_input = "unaligned_bc_tagged"
mapped_bam = bam_base + "{target}.bam"

BAM_MAP_LKUP = {}
BAM_REF_LKUP = {}
BAM_ANN_LKUP = {}
BAM_DEP_LKUP = {}
BAM_SYMLINKS = {}

BT2_MAPPINGS = []
STAR_MAPPINGS = []

def mapstr_to_targets(mapstr, inbam="uBAM", outbam="final"):
    """
    Converts a mapping strategy encoding, provided as a string, into a series of BAM file names
    and their dependencies. Rules are assigned by naming convention.
    Examples:

        [uBAM->]bowtie2:rRNA->STAR:genome:final

    uBAM stands for the CB and UMI tagged, pre-processed but unmapped reads.
    The '->' operator always extracts the unmapped fraction of reads from the bam file
    on the left.
    On the right of the operator is either a single name (of a BAM file),
    or <mapper>:<reference>. Optionally, a triplet can be used <mapper>:<reference>:<bam>
    where the presence of <bam> indicates that a BAM file with the assigned reads should be
    kept.

    If not otherwise specified, the last bam is by convention named or symlinked as 'final.bam'
    This file will be the input to downstream processing rules.

    [TODO:] Parallel and nested mappings can be implemented by parentheses and commata (implement via recursion):

        [uBAM->]bowtie2:rRNA->rRNA,STAR:genome->final

    NOTE: Gene tagging will be applied if annotation rules are provided for the associated mapping index
    """
    from collections import defaultdict
    # TODO: first tokenize by parenthesis. Optional
    targets = set()
    dep_lkup = {}
    ref_lkup = {}
    ann_lkup = {}
    mapper_lkup = {}
    chain_lengths = {}
    symlinks = {}

    def process(token, inbam):
        parts = token.split(":")
        if len(parts) == 1:
            target = parts[0]
            if target != inbam:
                dep_lkup[target] = inbam
                symlinks[target] = inbam

            return target
        else:
            print("parts", parts)
            if len(parts) == 2:
                mapper, ref = parts
                outname = f"{inbam}.{ref}.{mapper}"
            else:
                mapper, ref, linkname = parts
                outname = f"{inbam}.{ref}.{mapper}"
                symlinks[linkname] = outname

            # if ref.endswith("!"):
            #     ref = ref[:-1]
            ann_lkup[outname] = ref
            ref_lkup[outname] = ref
            mapper_lkup[outname] = mapper
            dep_lkup[outname] = inbam

            return outname

    # assume no ','
    left = inbam
    n = 0
    chain = [inbam] + mapstr.split("->") + ["final"]
    # print("  chain", chain)
    while chain:
        right = chain.pop(0)
        if left == right:
            continue

        for r in right.split(","):
            target = process(r, inbam=left)
            targets.add(target)

        left = target

    return sorted(targets), dep_lkup, ref_lkup, mapper_lkup, ann_lkup, symlinks


def get_mapped_BAM_output(default_strategy="STAR:genome:final"):
    """
    This function is called from main.smk at least once 
    to determine which output files need to be generated.
    """
    out_files = []

    for index, row in project_df.df.iterrows():
        mapstr = default_strategy
        if hasattr(row, "map_strategy") and row.map_strategy:
            mapstr = row.map_strategy

        print(index, mapstr)
        targets, dep_lkup, ref_lkup, mapper_lkup, ann_lkup, symlinks = mapstr_to_targets(mapstr, inbam=ubam_input)
        for target in targets:
            _target = mapped_bam.format(project_id=index[0], sample_id=index[1], target=target)
            BAM_DEP_LKUP[_target] = mapped_bam.format(project_id=index[0], sample_id=index[1], target=dep_lkup.get(target, None))
            BAM_MAP_LKUP[_target] = mapper_lkup.get(target, None)
            species_d = project_df.config.get_variable("species", name=row.species)
            ref = ref_lkup.get(target, None)
            if ref:
                BAM_REF_LKUP[_target] = species_d[ref]["sequence"]
                BAM_ANN_LKUP[_target] = species_d[ref]["annotation"]

            out_files.append(_target)

        for target, src in symlinks.items():
            _target = mapped_bam.format(project_id=index[0], sample_id=index[1], target=target)
            _src = mapped_bam.format(project_id=index[0], sample_id=index[1], target=src)
            BAM_SYMLINKS[_target] = _src

            out_files.append(_target)

        # SYMLINKS[mapped_bam.format(project_id=index[0], sample_id=index[1], target="uBAM")] = ubam_input.format(project_id=index[0], sample_id=index[1], target="uBAM")
        # out_files += expand(mapped_bam, project_id=index[0], sample_id=index[1], target="final")

    print("ALL TOP-LEVEL TARGETS")
    for o in out_files:
        print(
            f"  {o}\n"
            f"      input={BAM_DEP_LKUP.get(o, None)}"
            f"      mapper={BAM_MAP_LKUP.get(o, None)}"
            f"      ref={BAM_REF_LKUP.get(o, None)}"
            f"      ann={BAM_ANN_LKUP.get(o, None)}"
        )

    return out_files


def get_annotation(wc):
    return ""

def get_mapping_index(wc):
    return "index.bt2"

def get_mapping_parameters(wc):
    return "--local --ignore-quals --score-min=L,0,1.5 -L 10 -D 30 -R 30"

ruleorder: 
    map_reads_bowtie2 > map_reads_STAR > symlinks

# rule all:
#     input: get_mapped_BAM_output()

rule symlinks:
    input: lambda wc: BAM_SYMLINKS.get(f"{wc.name}.bam", 'na')
    output: "{name}.bam"
    shell:
        "ln -s {input} {output}"

rule map_reads_bowtie2:
    input:
        in_bam=lambda wc: BAM_DEP_LKUP[wc.name + '.bowtie2.bam']
    output:
        out_bam="{name}.bowtie2.bam"
    params:
        ref_ann=lambda wc: BAM_ANN_LKUP.get(wc.name + '.bowtie2.bam', 'na'),
        map_index=get_mapping_index,
        map_flags=get_mapping_parameters,
    threads: 32
    shell:
        "samtools fastq -f 4 {input.in_bam} | "
        "bowtie2 -p {threads} --reorder --mm {params.map_flags} -x {params.map_index} -b {input.in_bam} | "
        "samtools view --threads=2 -Sbuh /dev/stdin |"
        "{dropseq_tools}/TagReadWithGeneFunction I=/dev/stdin O={output.out_bam} ANNOTATIONS_FILE={params.ref_ann}"
        # "sambamba sort -t {threads} -m 8G --tmpdir=/tmp/tmp.{wildcards.name} -l 6 -o {output} /dev/stdin "

rule map_reads_STAR:
    input:
        in_bam=lambda wc: BAM_DEP_LKUP[wc.name + '.STAR.bam']
    output:
        out_bam="{name}.STAR.bam"
    params:
        ref_ann=get_annotation,
        map_idx=get_mapping_index,
        map_flags=get_mapping_parameters,
    threads: 32
    shell:
        "samtools fastq -f 4 {input.in_bam} | "
        " STAR <FILL IN STUFF> |"
        " sambamba sort -t {threads} -m 8G --tmpdir=/tmp/tmp.{wildcards.name} -l 6 -o {output} /dev/stdin "



# test_str = [
#     "uBAM->bowtie2:contaminants->STAR:genome!",
#     "uBAM->bowtie2:contaminants->cont->STAR:genome!->final",
#     "bowtie2:rRNA->rRNA,STAR:genome!->final",
# ]
# for mapstr in test_str:
#     print(">>>parsing", mapstr)
#     (
#         targets,
#         dep_lkup,
#         genome_lkup,
#         mapper_lkup,
#         tagging_lkup,
#     ) = mapstr_to_targets(mapstr)
#     for target in targets:
#         print(f"  {target} <- ({dep_lkup.get(target, None)})")
#         print(
#             f"    {mapper_lkup.get(target, None)}:{genome_lkup.get(target, None)} tagging={tagging_lkup.get(target, None)}"
#         )

# rule map_reads_bowtie2:
#     input:
#         unpack(get_species_genome_annotation),
#         unpack(get_star_input_bam),
#         unpack(get_bowtie2_index),
#         tagged_bam=tagged_bam
#     output:
#         star_log_file,
#         final_bam=final_bam
#     threads: 8
#     params:
#         tmp_dir = star_tmp_dir,
#         star_prefix = star_prefix
#     shell:
#         """
#         STAR \
#             --genomeLoad NoSharedMemory \
#             --genomeDir  {input.index} \
#             --sjdbGTFfile {input.annotation} \
#             --readFilesCommand samtools view \
#             --readFilesIn {input.reads} \
#             --readFilesType SAM SE \
#             --outFileNamePrefix {params.star_prefix} \
#             --outSAMprimaryFlag AllBestScore \
#             --outSAMattributes All \
#             --outSAMunmapped Within \
#             --outStd BAM_Unsorted \
#             --outSAMtype BAM Unsorted \
#             --limitOutSJcollapsed 5000000 \
#             --runThreadN {threads} | \
#             python {repo_dir}/scripts/fix_bam_header.py \
#                 --in-bam-star /dev/stdin \
#                 --in-bam-tagged {input.tagged_bam} \
#                 --out-bam /dev/stdout | \
#             {dropseq_tools}/TagReadWithGeneFunction \
#                 I=/dev/stdin \
#                 O={output.final_bam} \
#                 ANNOTATIONS_FILE={input.annotation}

#         rm -rf {params.tmp_dir}
#         """
