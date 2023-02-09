from spacemake.contrib import __license__, __version__, __author__, __email__

ruleorder: link_raw_reads > link_demultiplexed_reads 

rule demultiplex_data:
    params:
        demux_barcode_mismatch=lambda wildcards: int(project_df.get_metadata('demux_barcode_mismatch', demux_dir = wildcards.demux_dir)),
        sample_sheet=lambda wildcards: project_df.get_metadata('sample_sheet', demux_dir = wildcards.demux_dir),
        output_dir= lambda wildcards: expand(demux_dir_pattern, demux_dir=wildcards.demux_dir)
    input:
        lambda wildcards: project_df.get_metadata('basecalls_dir', demux_dir = wildcards.demux_dir)
    output:
        demux_indicator
    threads: 8
    shell:
        """
        retVal=$(bcl2fastq \
            --no-lane-splitting --fastq-compression-level=9 \
            --mask-short-adapter-reads 15 \
            --barcode-mismatch {params.demux_barcode_mismatch} \
            --output-dir {params.output_dir} \
            --sample-sheet {params.sample_sheet} \
            --runfolder-dir {input} \
            -p {threads})
        
        if [ $retVAL==0 ]; then
            echo "demux finished: $(date)" > {output}
        fi
        exit $retVal
        """


rule link_demultiplexed_reads:
    input:
        ancient(unpack(get_demux_indicator))
    output:
        raw_reads_pattern
    params:
        demux_dir = lambda wildcards: expand(demux_dir_pattern,
            demux_dir = project_df.get_metadata('demux_dir', sample_id = wildcards.sample_id,
                                     project_id = wildcards.project_id)),
        reads_folder = raw_data_illumina_reads
    shell:
        """
        mkdir -p {params.reads_folder}

        find {params.demux_dir} -type f -wholename '*/{wildcards.sample_id}/*R{wildcards.mate}*.fastq.gz' -exec ln -sr {{}} {output} \; 
        """


rule link_raw_reads:
    input:
        unpack(get_reads)
    output:
        raw_reads_pattern
    run:
        if len(input) == 1:
            # either link raw reads
            shell("ln -rs {input} {output}")
        else:
            # or append reads together
            shell("cat {input} > {output}")
            

rule zcat_pipe:
    input: "{name}.fastq.gz"
    output: pipe("{name}.fastq")
    threads: 1
    shell: "zcat {input} >> {output}"


rule tag_reads_bc_umi:
    input:
        # these implicitly depend on the raw reads via zcat_pipes
        ## TODO: make R1 optional for bulk samples
        # R1 = raw_reads_mate_1,
        # R2 = raw_reads_mate_2
        ancient(unpack(get_linked_reads))
    params:
        bc_params = get_bc_preprocess_settings,
        read1 = lambda wc: get_linked_reads(wc).get("R1", None),
        read2 = lambda wc: get_linked_reads(wc).get("R2", None),
        reads = lambda wc: get_linked_reads(wc).get("reads", None)
    output:
        ubam = tagged_bam,
        bc_stats = stats_dir + "/{sample_id}.bc_stats.tsv"
        # bc_counts = barcode_readcounts
    log:
        log_dir + '/fastq_to_uBAM.log'
    threads: max(1, min(workflow.cores - 2, 8)) # reserve two cores for the input zcat pipes
    shell:
        "python {bin_dir}/fastq_to_uBAM.py "
        "  --sample={wildcards.sample_id} "
        "  --log-level={log_level}"
        "  --log-file={log} "
        "  --debug={log_debug} "
        "  --read1={params.read1} "
        "  --read2={params.read2} "
        "  --matrix={params.reads} "
        "  --parallel={threads} "
        "  --save-stats={output.bc_stats} "
        # "  --save-cell-barcodes={output.bc_counts} "
        "  --out-bam=/dev/stdout "
        "  {params.bc_params} "
        " | samtools view -bh --threads=4 /dev/stdin > {output.ubam} "


rule trim_adapters_polyA:
    input:
        tagged_bam
    output:
        trimmed = tagged_polyA_adapter_trimmed_bam, # temp(
        stats = stats_dir + "/cutadapt.csv",
    log: log_dir + '/cutadapt.log'
    params:
        adapter_flavor = lambda wildcards: project_df.get_metadata(
            "adapter_flavor", 
            project_id=wildcards.project_id,
            sample_id=wildcards.sample_id
        )
    threads: 12
    shell:
        "python {bin_dir}/cutadapt_bam.py {input} "
        "  --sample={wildcards.sample_id} "
        "  --log-level={log_level} "
        "  --log-file={log} "
        "  --debug={log_debug} "
        "  --config={spacemake_config} "
        "  --bam-out={output.trimmed} --bam-out-mode=b "
        "  --stats-out={output.stats} "
        "  --adapter-flavor={params.adapter_flavor} "
        "  --threads-write=4 " # scales up to 4 bc we use compression
        "  --threads-work={threads} "

rule run_fastqc:
    input:
        # we need to use raw reads here, as later during "reversing" we do the umi
        # extraction, and barcode identification (for the combinatorial barcoding)
        # in order for R1 to have the same pattern
        raw_reads_pattern
    output:
        fastqc_pattern
    params:
        output_dir = fastqc_root 
    threads: 4
    shell:
        """
        mkdir -p {params.output_dir}

        fastqc -t {threads} -o {params.output_dir} {input}
        """
# rule filter_mm_reads:
#     input:
#         unpack(get_final_bam)
#     output:
#         pipe(final_bam_mm_included_pipe)
#     shell:
#         """
#         python {repo_dir}/scripts/filter_mm_reads.py \
#             --in-bam {input} \
#             --out-bam {output}
#         """
