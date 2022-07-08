ubam_input = "unaligned_bc_tagged.polyA_adapter_trimmed" # this must be local and not have .bam appended!
# basically, if ubam_input were used as {target} in mapped.bam it should eval to tagged_polyA_adapter_trimmed_bam
mapped_bam = complete_data_root + "/{target}.bam"
star_mapped_bam = complete_data_root + "/{target}.STAR.bam"
bt2_mapped_bam = complete_data_root + "/{target}.bowtie2.bam"

def wc_fill(x, wc):
    return x.format(sample_id=wc.sample_id, project_id=wc.project_id, target=wc.target)

# These lookup dictionaries use the absolute paths of BAM files as keys
BAM_DEP_LKUP = {} # input on which a mapped BAM will depend (usually a uBAM or BAM of previous stage)
BAM_MAP_LKUP = {} # which mapper to use (STAR or bowtie2)
BAM_IDX_LKUP = {} # mapping index which is needed to create a BAM
BAM_REF_LKUP = {} # reference sequence
BAM_ANN_LKUP = {} # annotation file used for tagging 
BAM_SYMLINKS = {} # contains the source file to symlink to

BAM_MAP_FLAGS_LKUP = {}

default_BT2_MAP_FLAGS = (
    " --local"
    " -L 10 -D 30 -R 30"
    " --ignore-quals"
    " --score-min=L,0,1.5" # require 75% of perfect match (2=base match)
)

default_STAR_MAP_FLAGS = (
    " --genomeLoad NoSharedMemory"
    " --outSAMprimaryFlag AllBestScore"
    " --outSAMattributes All"
    " --outSAMunmapped Within"
    " --outStd BAM_Unsorted"
    " --outSAMtype BAM Unsorted"
    " --limitOutSJcollapsed 5000000"
)

def mapstr_to_targets(mapstr, inbam="uBAM", outbam="final"):
    """
    Converts a mapping strategy provided as a string into a series of BAM file names
    and their dependencies. Rule matching is otherwise ensured by the convention of the 
    target BAM filename ending in "STAR.bam" or "bowtie2.bam".

    Examples:

        [uBAM->]bowtie2:rRNA->STAR:genome:final

    uBAM stands for the CB and UMI tagged, pre-processed, but unmapped reads.
    The '->' operator always extracts the unmapped fraction of reads from the bam file
    on the left.
    On the right of the operator are the main parameters for the mapping rule <mapper>:<reference>.
    The target BAM will have the name <input>.<reference>.<mapper>.bam.
    Optionally, a triplet can be used <mapper>:<reference>:<symlink> where the presence of <symlink> 
    indicates that the BAM file with the assigned reads should be accessible under a useful 
    shorthand name as well ("final.bam" for instance).

    The example above is going to create
        
        (1) uBAM.rRNA.bowtie2.bam
        
            using bowtie2 and the index associated with the "rRNA" reference


        (2) uBAM.rRNA.bowtie2.genome.STAR.bam 
        
            using STAR on the *unmapped* reads from BAM (1)


        (3) final.bam
        
            a symlink pointing to BAM (2).

    If not otherwise specified, the last bam is by convention named or symlinked as 'final.bam'
    This file will be the input to downstream processing rules (for DGE creation).

    NOTE: Parallel mappings can be implemented by using commata:

        [uBAM->]bowtie2:rRNA:rRNA,STAR:genome:final

        This rule differs from the first example because it will align the unmapped reads from the uBAM
        in parallel to the rRNA reference and to the genome. In this way the same reads can match to both
        indices.

    NOTE: Gene tagging will be applied if annotation rules are provided for the associated mapping index
    """
    from collections import defaultdict

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
            else:
                symlinks[target] = inbam

            return target
        else:
            if len(parts) == 2:
                mapper, ref = parts
                outname = f"{inbam}.{ref}.{mapper}"
            else:
                mapper, ref, linkname = parts
                outname = f"{inbam}.{ref}.{mapper}"
                symlinks[linkname] = outname

            ann_lkup[outname] = ref
            ref_lkup[outname] = ref
            mapper_lkup[outname] = mapper
            dep_lkup[outname] = inbam

            return outname

    left = inbam
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

        # print(index, mapstr)
        targets, dep_lkup, ref_lkup, mapper_lkup, ann_lkup, symlinks = mapstr_to_targets(mapstr, inbam=ubam_input)
        for target in targets:
            _target = mapped_bam.format(project_id=index[0], sample_id=index[1], target=target)
            # print(f"TARGET: target={target} -> _target={_target}")
            BAM_DEP_LKUP[_target] = mapped_bam.format(project_id=index[0], sample_id=index[1], target=dep_lkup.get(target, None))
            mapper = mapper_lkup.get(target, None)
            BAM_MAP_LKUP[_target] = mapper
            species_d = project_df.config.get_variable("species", name=row.species)
            ref = ref_lkup.get(target, None)
            if ref:
                BAM_REF_LKUP[_target] = species_d[ref]["sequence"]
                BAM_ANN_LKUP[_target] = species_d[ref].get("annotation", None)
                if mapper == "bowtie2":
                    BAM_MAP_FLAGS_LKUP[_target] = species_d[ref].get("bt2_flags", default_BT2_MAP_FLAGS)
                    BAM_IDX_LKUP[_target] = species_d[ref].get("bt2_index", BAM_REF_LKUP[_target])
                else:
                    BAM_MAP_FLAGS_LKUP[_target] = species_d[ref].get("STAR_flags", default_STAR_MAP_FLAGS)
                    BAM_IDX_LKUP[_target] = species_d[ref].get("index_dir", "NO_STAR_index_IN_config.yaml")
                
            out_files.append(_target)

        for target, src in symlinks.items():
            _target = mapped_bam.format(project_id=index[0], sample_id=index[1], target=target)
            _src = mapped_bam.format(project_id=index[0], sample_id=index[1], target=src)
            BAM_SYMLINKS[_target] = _src

            out_files.append(_target)

    print("ALL TOP-LEVEL TARGETS")
    for o in out_files:
        print(
            f"   {o}\n"
            f"      input={BAM_DEP_LKUP.get(o, None)}\n"
            f"      mapper={BAM_MAP_LKUP.get(o, None)}\n"
            f"      ref={BAM_REF_LKUP.get(o, None)}\n"
            f"      ann={BAM_ANN_LKUP.get(o, None)}\n"
        )

    return out_files


def get_annotation_command(output):
    tagging_cmd =  "| {dropseq_tools}/TagReadWithGeneFunction I=/dev/stdin O={output.bam} ANNOTATIONS_FILE={ann_gtf}"
    ann = BAM_ANN_LKUP.get(output, None)
    if ann and ann.lower().endswith(".gtf"):
        return tagging_cmd.format(dropseq_tools=dropseq_tools, output=output, ann_gtf=ann)
    else:
        return f"| samtools view --threads=2 -bh /dev/stdin > {output.bam}"

    # TODO: add other means of annotation

# def get_symlink_source(wc):
#     output = complete_data_root + f"/{wc.name}.bam"
#     if not output in BAM_SYMLINKS:
#         for k,v in sorted(BAM_SYMLINKS.items()):
#             print(f"{k} <- {v}")
#             if not v in BAM_DEP_LKUP:
#                 print("!!! NOT IN BAM_DEP_LKUP")
#                 for k in sorted(BAM_DEP_LKUP.keys()):
#                     print(k)

#     return BAM_SYMLINKS.get(output, f'NA {output} NA')

# ruleorder: 
#     map_reads_bowtie2 > map_reads_STAR > symlinks

rule symlinks:
    input: lambda wc: BAM_SYMLINKS[wc_fill(mapped_bam, wc)]
    output: mapped_bam
    shell:
        "ln -s {input} {output}"

rule map_reads_bowtie2:
    input:
        bam=lambda wc: BAM_DEP_LKUP[wc_fill(bt2_mapped_bam, wc)],
        index=lambda wc: BAM_IDX_LKUP[wc_fill(bt2_mapped_bam, wc)],
    output:
        bam=bt2_mapped_bam
    log: bt2_mapped_bam + ".log"
    params:
        annotation_cmd=lambda wildcards, output: get_annotation_command(output),
        flags=lambda wildcards, output: BAM_MAP_FLAGS_LKUP[output.bam],
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
        "  -x {input.index} -b /dev/stdin --preserve-tags"
        "  {params.flags} 2> {log}"
        " "
        # fix the BAM header to accurately reflect the entire history of processing via PG records.
        "| python {repo_dir}/scripts/splice_bam_header.py"
        "  --in-ubam {input.bam}"
        " "
        "{params.annotation_cmd}"
        # "sambamba sort -t {threads} -m 8G --tmpdir=/tmp/tmp.{wildcards.name} -l 6 -o {output} /dev/stdin "


rule map_reads_STAR:
    input: 
        bam=lambda wc: BAM_DEP_LKUP.get(wc_fill(star_mapped_bam, wc), f"can't_find_bam_{wc}"),
        index=lambda wc: BAM_IDX_LKUP.get(wc_fill(star_mapped_bam, wc), f"can't find_idx_{wc}"),
    output:
        bam=star_mapped_bam
    threads: 32
    params:
        annotation_cmd=lambda wildcards, output: get_annotation_command(output),
        annotation=lambda wilcards, output: BAM_ANN_LKUP[output.bam],
        flags=lambda wildcards, output: BAM_MAP_FLAGS_LKUP[output.bam],
        tmp_dir = star_tmp_dir,
        star_prefix = star_prefix
    shell:
        "STAR {params.flags}"
        " --genomeDir {input.index}"
        " --readFilesIn {input.bam}"
        " --readFilesCommand 'samtools view -f 4'"
        " --readFilesType SAM SE"
        " --sjdbGTFfile {params.annotation}"
        " --outFileNamePrefix {params.star_prefix}"
        " --runThreadN {threads}"
        " "
        "| python {repo_dir}/scripts/splice_bam_header.py"
        " --in-ubam {input.bam}"
        " "
        "{params.annotation_cmd}"
        " "
        "rm -rf {params.tmp_dir}"
