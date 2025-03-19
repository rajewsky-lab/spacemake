import os
import subprocess
import time
import yaml

from spacemake.project_df import get_global_ProjectDF
from spacemake.util import sync_timestamps

import sys


def find_bam_files(folder):
    """
    Finds all .bam files in the given folder and checks if any of them is a symlink.
    
    Returns a list of tuples of type (str, bool), e.g. ('bam_file', False)
    """
    if not os.path.isdir(folder):
        raise ValueError(f"The provided path {folder} is not a valid directory.")
    
    # Find files and check for symlinks
    bam_files = [f for f in os.listdir(folder) if f.endswith('.bam')]
    bam_file_paths = [os.path.join(folder, f) for f in bam_files]
    are_symlinks = [os.path.islink(bam_file_path) for bam_file_path in bam_file_paths]

    return list(zip(bam_file_paths, are_symlinks))


def get_map_strategy_sequences(project_id, sample_id):
    """
    Returns a dictionary of reference_types and their location, e.g. {rRNA : /path/to/disk/sequence.fa}
    """
    pdf = get_global_ProjectDF()

    map_strategy = pdf.get_sample_info(project_id, sample_id)['map_strategy']
    sequence_type = [mapping.split(':')[1] for mapping in map_strategy.split('->')]

    with open("config.yaml") as yamlfile:
        cf = yaml.safe_load(yamlfile.read())
    sample_species = pdf.get_sample_info(project_id, sample_id)['species']

    reference_type = {st : cf['species'][sample_species][st]['sequence'] for st in sequence_type}

    return reference_type


def check_if_all_files_exist(project_id, sample_id, file_type):
    """
    Checks if all required files exist for the map_strategy of the sample.

    Args:
        file_type (str): The type of file to check. Must be either 'bam' or 'cram'

    Returns:
        bool: True if all required files exist. False otherwise.
    """
    project_folder = os.path.join('projects', project_id, 'processed_data', sample_id, 'illumina', 'complete_data')

    pdf = get_global_ProjectDF()
    map_strategy = pdf.get_sample_info(project_id, sample_id)['map_strategy']

    aligner = [mapping.split(':')[0] for mapping in map_strategy.split('->')]
    sequence_type = [mapping.split(':')[1] for mapping in map_strategy.split('->')]

    files_expected = [f"{x}.{y}.{file_type}" for x, y in zip(sequence_type, aligner)]

    all_files_exist = True
    for file in files_expected:
        file_path = os.path.join(project_folder, file)
        if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
            all_files_exist = False
            print('file', file_path, 'does not exist or is empty.')
            break

    return all_files_exist


def convert_bam_to_cram(project_id, sample_id, threads=4):
    """
    Converts all BAM files to CRAM and updates the timestamps to those of the
    original files. Symbolic links are treated as such.
    """
    species_sequences = get_map_strategy_sequences(project_id, sample_id)

    if check_if_all_files_exist(project_id, sample_id, 'bam'):
        print('All required BAM files exist on disk -- will continue with CRAM conversion.')
    else:
        print('Not all required BAM files were found on disk -- will not continue with CRAM conversion.')
        return

    project_folder = os.path.join('projects', project_id, 'processed_data',
                                  sample_id, 'illumina', 'complete_data')    
    bam_files = find_bam_files(project_folder)
    
    for idx in range(len(bam_files)):
        bam_filename, bam_file_is_symlink = bam_files[idx]
        bam_filename_prefix = bam_filename.rsplit('.', 1)[0]
        cram_filename = bam_filename_prefix + '.cram'

        if os.path.exists(cram_filename) and os.path.getsize(cram_filename) > 0:
            print('CRAM file', cram_filename, 'already exists. Skipping conversion.')
            continue
        elif os.path.exists(cram_filename) and os.path.getsize(cram_filename) == 0:
            print('WARNING: CRAM file', cram_filename, 'already exists, but file size is 0. Possible permission error.')

        if 'unaligned' in bam_filename:
            # Special case for now
            print('Converting', bam_filename, 'to',
                  os.path.join(project_folder, 'unaligned_bc_tagged.polyA_adapter_trimmed.cram'),
                  '...', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            subprocess.run(
                [
                    "samtools", "view",
                    "-T", species_sequences['genome'],
                    "-C",
                    "--threads", str(threads),
                    "-o",  os.path.join(project_folder, 'unaligned_bc_tagged.polyA_adapter_trimmed.cram'),
                    bam_filename
                ])

            continue
        
        if bam_file_is_symlink:
            true_bam_filename = os.readlink(bam_filename)
            true_bam_filename_prefix = true_bam_filename.rsplit('.', 1)[0]
            try:
                os.symlink(true_bam_filename_prefix + '.cram', cram_filename)
            except FileExistsError:
                print('CRAM symlink already exists.')
        else:
            print('Converting', bam_filename, 'to', cram_filename, 
            '...', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

            for ref_type in species_sequences:
                if ref_type in bam_filename:
                    ref_sequence = species_sequences[ref_type]
                    break

            subprocess.run(
                [
                    "samtools", "view",
                    "-T", ref_sequence,
                    "-C",
                    "--threads", str(threads),
                    "-o", cram_filename,
                    bam_filename
                ])

        sync_timestamps(bam_filename, cram_filename)


def remove_bam_files(project_folder):
    file_sizes_cram = [
        os.path.getsize(os.path.join(project_folder, f))
        for f in os.listdir(project_folder)
        if f.endswith(".cram") and os.path.isfile(os.path.join(project_folder, f))
        ]
    
    file_sizes_bam = [
        os.path.getsize(os.path.join(project_folder, f))
        for f in os.listdir(project_folder)
        if f.endswith(".bam") and os.path.isfile(os.path.join(project_folder, f))
        ]
    
    # Remove BAM files
    while True:
        response = input("This action will permanently delete the BAM files from the sample. "
        "Are you sure that all necessary CRAM files have been created? [y/n]: ").strip().lower()
        if response in ['y', 'yes']:
            # Remove files
            bam_files = find_bam_files(project_folder)
            for bam_file in bam_files:
                os.remove(bam_file[0])
            print("BAM files deleted.")
            print("Total disk space saved through the migration:", 
                  f"{round((sum(file_sizes_bam)-sum(file_sizes_cram))/(1024*1024)):,}", "MB")
            return False
        elif response in ['n', 'no']:
            print('Deletion aborted. Please run spacemake migrate')
            return False
        else:
            print("Please enter 'y' or 'n'.")

