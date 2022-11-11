#########
# about #
#########
__version__ = '0.1.0'
__author__ = ['Nikos Karaiskos', 'Tamas Ryszard Sztanka-Toth']
__licence__ = 'GPL'
__email__ = ['nikolaos.karaiskos@mdc-berlin.de', 'tamasryszard.sztanka-toth@mdc-berlin.de']

###################################################
# Snakefile containing the dropseq pipeline rules #
###################################################
rule remove_polyA:
    input:
        tagged_bam
    output:
        trimmed = tagged_polyA_adapter_trimmed_bam, # temp(
        stats = reports_dir + "/cutadapt.csv"
    params:
        reports_dir = reports_dir,
        adap3 = os.path.join(spacemake_dir, "data/adapters.fa"),
        # adap5 = os.path.join(spacemake_dir, "data/adapters5.fa"),
    threads: 10
    shell:
        "python {repo_dir}/scripts/cutadapt_bam.py {input} "
        " --bam-out={output.trimmed} --bam-out-mode=b "
        " --stats-out={output.stats} "
        " --adapters-right={params.adap3} "
        # " --adapters-left={params.adap5} "
        " --threads-write=4 " # scales up to 4 bc we use compression
        " --threads-work={threads} "

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
