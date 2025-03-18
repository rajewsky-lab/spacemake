import os
import subprocess
import time
import yaml

from spacemake.project_df import get_global_ProjectDF
from spacemake.util import sync_timestamps


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
     

def convert_bam_to_cram(project_id, sample_id, threads=4):
    """
    Converts all BAM files to CRAM and updates the timestamps to those of the
    original files. Symbolic links are treated as such.
    """
    species_sequences = get_map_strategy_sequences(project_id, sample_id)

    project_folder = os.path.join('projects', project_id, 'processed_data',
                                  sample_id, 'illumina', 'complete_data')    
    bam_files = find_bam_files(project_folder)
    
    for idx in range(len(bam_files)):
        bam_filename, bam_file_is_symlink = bam_files[idx]
        bam_filename_prefix = bam_filename.rsplit('.', 1)[0]
        cram_filename = bam_filename_prefix + '.cram'

        if os.path.exists(cram_filename):
            print('CRAM file', cram_filename, 'already exists. Skipping conversion.')
            continue

        if 'unaligned' in bam_filename:
            continue
        
        if bam_file_is_symlink:
            true_bam_filename = os.readlink(bam_filename)
            true_bam_filename_prefix = true_bam_filename.rsplit('.', 1)[0]
            os.symlink(true_bam_filename_prefix + '.cram', cram_filename)
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
                ]
            )

        sync_timestamps(bam_filename, cram_filename)


def remove_files(project_folder):
    # - BAM files (only if CRAMs are present)
    bam_files = find_bam_files(project_folder)

    # - unaligned.bam

    # remove tiles

