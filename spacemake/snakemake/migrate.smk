__version__ = "0.1"
__author__ = ["Nikos Karaiskos"]
__license__ = "GPL"

from spacemake.migrate import (
    check_if_all_files_exist,
    convert_bam_to_cram,
    rename_log_files,
    remove_bam_files,
    update_version_in_config
    )

include: "scripts/snakemake_helper_functions.py"
include: "variables.py"

project_folders = config["project_folders"]
projects = config["projects"]
samples = config["samples"]

bam_cleanup_markers = [
    f"{folder}/bam_files_removed.txt" for folder in config["project_folders"]
]

rule all:
    input:
        config["targets"],
        "version_update_done.txt"

rule convert_bam_to_cram:
    output:
        "{project_folder}/final.polyA_adapter_trimmed.cram"
    threads: 4
    run:
        project_folder = wildcards.project_folder
        project_id = project_folder.split("/")[1]
        sample_id = project_folder.split("/")[3]

        if check_if_all_files_exist(project_id, sample_id, "cram"):
            print(f"CRAM {output[0]} already exists. Skipping conversion.")
        else:
            convert_bam_to_cram(project_id, sample_id, threads)

        rename_log_files(project_id, sample_id)

rule remove_bam_files:
    input:
        "{project_folder}/final.polyA_adapter_trimmed.cram"
    output:
        "{project_folder}/bam_files_removed.txt"
    run:
        remove_bam_files(wildcards.project_folder, output[0])

rule update_config:
    input:
        bam_cleanup_markers,
        expand("{project_folder}/.unmapped_removed", project_folder=project_folders)
    output:
        touch("version_update_done.txt")
    run:
        print("Updating version in config.yaml ...")
        update_version_in_config()

rule generate_dummy_files:
    input:
        bam_cleanup_markers
    output:
        "{project_folder}/.unmapped_removed"
    shell:
        "touch {output}"
