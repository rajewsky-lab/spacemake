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
triplets = config["project_sample_folder_triplets"]

rule all:
    input:
        config["targets"]

rule convert_bam_to_cram:
    output:
        "{project_folder}/{project_id}/processed_data/{sample_id}/illumina/complete_data/final.polyA_adapter_trimmed.cram"
    run:
        print("Beginning migration ...")

        project_id = wildcards.project_id
        sample_id = wildcards.sample_id
        project_folder = wildcards.project_folder

        if check_if_all_files_exist(project_id, sample_id, "cram"):
            print("CRAM already exists. Checking for BAM files to remove...")
            bam_files = find_bam_files(project_folder)
            if bam_files:
                remove_bam_files(project_folder)
            else:
                print("No BAMs found.")
        else:
            convert_bam_to_cram(project_id, sample_id, threads)
            remove_bam_files(project_folder)
