__version__ = "0.1"
__author__ = ["Nikos Karaiskos"]
__license__ = "GPL"

from spacemake.project_df import ProjectDF
from spacemake.config import ConfigFile
from spacemake.errors import SpacemakeError

include: "migrate.py"
include: "scripts/snakemake_helper_functions.py"
include: "variables.py"

pdf = ProjectDF(config["project_df"], ConfigFile.from_yaml("config.yaml"))

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
            print("CRAM already exists. Skipping conversion.")
        else:
            convert_bam_to_cram(project_id, sample_id, threads)

rule remove_bam_files:
    input:
        "{project_folder}/final.polyA_adapter_trimmed.cram"
    output:
        "{project_folder}/bam_files_removed.txt"
    run:
        remove_bam_files(wildcards.project_folder, output[0])

rule update_version_in_config:
    input:
        bam_cleanup_markers  # <- just a plain list
    output:
        touch("version_update_done.txt")
    run:
        print("Updating version in config.yaml ...")






# rule convert_bam_to_cram:
#     output:
#         "{project_folder}/{project_id}/processed_data/{sample_id}/illumina/complete_data/final.polyA_adapter_trimmed.cram"
#     run:
#         print("Beginning migration ...")

#         project_id = wildcards.project_id
#         sample_id = wildcards.sample_id
#         project_folder = wildcards.project_folder

#         if check_if_all_files_exist(project_id, sample_id, "cram"):
#             print("CRAM already exists. Checking for BAM files to remove...")
#             bam_files = find_bam_files(project_folder)
#             if bam_files:
#                 remove_bam_files(project_folder)
#             else:
#                 print("No BAMs found.")
#         else:
#             convert_bam_to_cram(project_id, sample_id, threads)
#             remove_bam_files(project_folder)
