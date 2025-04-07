import os
import subprocess
import sys
import time
import yaml

from spacemake.contrib import __version__
from spacemake.project_df import get_global_ProjectDF
from spacemake.util import sync_timestamps, sync_symlink_mtime

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
    import sys
    pdf = get_global_ProjectDF()

    species = pdf.get_sample_info(project_id, sample_id)['species']
    map_strategy = pdf.get_sample_info(project_id, sample_id)['map_strategy']
    sequence_type = [mapping.split(':')[1] for mapping in map_strategy.split('->')]

    reference_type = {st : f"species_data/{species}/{st}/sequence.fa" for st in sequence_type}

    # with open("config.yaml") as yamlfile:
    #     cf = yaml.safe_load(yamlfile.read())
    

    # reference_type = {st : cf['species'][sample_species][st]['sequence'] for st in sequence_type}

    return reference_type


def check_if_all_files_exist(project_id, sample_id, file_type):
    """
    Checks if all required files exist for the map_strategy of the sample.

    Args:
        file_type (str): The type of file to check. Must be either 'bam' or 'cram'

    Returns:
        bool: True if all required files exist. False otherwise.
    """
    project_folder = os.path.join('projects', project_id, 'processed_data',
                                  sample_id, 'illumina', 'complete_data')

    pdf = get_global_ProjectDF()
    sample_is_merged = pdf.get_sample_info(project_id, sample_id)['is_merged']
    map_strategy = pdf.get_sample_info(project_id, sample_id)['map_strategy']

    aligner = [mapping.split(':')[0] for mapping in map_strategy.split('->')]
    sequence_type = [mapping.split(':')[1] for mapping in map_strategy.split('->')]

    files_expected = [f"{x}.{y}.{file_type}" for x, y in zip(sequence_type, aligner)]

    if sample_is_merged:
        files_expected = ['final.polyA_adapter_trimmed.merged.' + file_type]

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
    pdf = get_global_ProjectDF()
    sample_is_merged = pdf.get_sample_info(project_id, sample_id)['is_merged']
    
    species_sequences = get_map_strategy_sequences(project_id, sample_id)

    if check_if_all_files_exist(project_id, sample_id, 'bam'):
        print('All required BAM files exist on disk -- will continue with CRAM conversion.')
    else:
        print('Not all required BAM files were found on disk -- will not continue with CRAM conversion.')
        return

    project_folder = os.path.join('projects', project_id, 'processed_data',
                                  sample_id, 'illumina', 'complete_data')    
    bam_files = find_bam_files(project_folder)
    bam_files = sorted(bam_files, key=lambda x: ("final" not in x[0], x[1])) # sort to make sure final is always checked first
    
    # simple dictionary to keep track of which ref_type is final
    ref_type_final = {}
    
    for idx in range(len(bam_files)):
        bam_filename, bam_file_is_symlink = bam_files[idx]
        bam_filename_prefix = bam_filename.rsplit('.', 1)[0]
        cram_filename = bam_filename_prefix + '.cram'

        if os.path.exists(cram_filename) and os.path.getsize(cram_filename) > 0:
            print('CRAM file', cram_filename, 'already exists. Skipping conversion.')
            continue
        elif os.path.exists(cram_filename) and os.path.getsize(cram_filename) == 0:
            print('WARNING: CRAM file', cram_filename, 'already exists, but file size is 0. Possible permission error.')

        # Special case for now
        if 'unaligned' in bam_filename:
            cram_filename = os.path.join(project_folder, 'unaligned_bc_tagged.polyA_adapter_trimmed.cram')
            if os.path.exists(cram_filename) and os.path.getsize(cram_filename) > 0:
                print('CRAM file', cram_filename, 'already exists. Skipping conversion.')
                continue
            else:
                print('Converting', bam_filename, 'to', cram_filename,
                    '...', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
                subprocess.run(
                    [
                        "samtools", "view",
                        "-T", species_sequences['genome'],
                        "-C",
                        "--threads", str(threads),
                        "-o", cram_filename,
                        bam_filename
                    ])
                sync_timestamps(bam_filename, cram_filename)
                continue
        
        if bam_file_is_symlink:
            true_bam_filename = os.readlink(bam_filename)
            true_bam_filename_prefix = true_bam_filename.rsplit('.', 1)[0]

            # Keep track of what's final to split the not_ files later.
            if 'final' in bam_filename:
                ref_type_final[true_bam_filename.rsplit('.')[0]] = True

            try:
                os.symlink(true_bam_filename_prefix + '.cram', cram_filename)
                sync_symlink_mtime(bam_filename, cram_filename)

            except FileExistsError:
                print('CRAM symlink already exists.')

        else:
            print('Converting', bam_filename, 'to', cram_filename, 
            '...', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

            for ref_type in species_sequences:
                if sample_is_merged:
                    ref_sequence = species_sequences['genome']
                    break
                if ref_type in bam_filename:
                    ref_sequence = species_sequences[ref_type]
                    break
            
            # warn if reference is gzip-compressed
            if ref_sequence.endswith('.gz'):
                print(f"ERROR: Reference file '{ref_sequence}' is gzip-compressed (.gz), which samtools cannot use for CRAM conversion.", file=sys.stderr)
                print("Please use an uncompressed FASTA (.fa) or bgzip-compressed version with proper indexing (.fai and .gzi).", file=sys.stderr)
                sys.exit(1)

            if ref_type in ref_type_final:
                #split files appropriately
                subprocess.run(
                    [
                        "samtools", "view",
                        "-T", ref_sequence,
                        "-C",
                        "-F 4",
                        "--threads", str(threads),
                        "-o", cram_filename,
                        bam_filename
                    ])
                
                directory, filename = os.path.split(cram_filename)
                not_filename = 'not_' + filename

                subprocess.run(
                    [
                        "samtools", "view",
                        "-T", ref_sequence,
                        "-C",
                        "-f 4",
                        "--threads", str(threads),
                        "-o", os.path.join(directory, not_filename),
                        bam_filename
                    ])
                
                sync_timestamps(bam_filename, os.path.join(directory, not_filename))

            else:
                subprocess.run(
                    [
                        "samtools", "view",
                        "-T", ref_sequence,
                        "-C",
                        "-F 4",
                        "--threads", str(threads),
                        "-o", cram_filename,
                        bam_filename
                    ])

                if sample_is_merged:
                    subprocess.run(
                        [
                            "touch", os.path.join(project_folder, 'final.polyA_adapter_trimmed.cram')
                        ])
                
            sync_timestamps(bam_filename, cram_filename)


def rename_log_files(project_id, sample_id):
    """
    Rename any .bam.log files (created by bowtie) so that qc_sheets and other reports
    generation downstream does not fail.
    """
    project_folder = os.path.join('projects', project_id, 'processed_data',
                                  sample_id, 'illumina', 'complete_data')    
    
    log_files = [f for f in os.listdir(project_folder) if f.endswith('.bam.log')]

    for log_filename in log_files:
        log_filename_prefix = log_filename.rsplit('.', 2)[0]
        new_log_filename = log_filename_prefix + '.cram.log'
        os.rename(os.path.join(project_folder, log_filename),
                  os.path.join(project_folder, new_log_filename))


def remove_bam_files(project_folder, output_file_path):
    bam_files = find_bam_files(project_folder)

    cram_files = [
        os.path.join(project_folder, f)
        for f in os.listdir(project_folder)
        if f.endswith(".cram") and os.path.isfile(os.path.join(project_folder, f))
    ]

    total_bam_size = sum(
        os.path.getsize(bam[0]) for bam in bam_files if os.path.exists(bam[0])
    )
    total_cram_size = sum(
        os.path.getsize(cram) for cram in cram_files if os.path.exists(cram)
    )

    # remove BAM files
    deleted_files = []
    for bam in bam_files:
        if os.path.exists(bam[0]):
            os.remove(bam[0])
            deleted_files.append(bam[0])

    saved_bytes = total_bam_size - total_cram_size
    saved_gb = saved_bytes / (1024 ** 3)

    # write report
    with open(output_file_path, "w") as out:
        out.write(f"BAM files removed: {len(deleted_files)}\n")
        out.write("Files:\n")
        for f in deleted_files:
            out.write(f"  - {f}\n")
        out.write(f"Total disk space saved: {saved_gb:.2f} GB\n")

    print(f"Deleted {len(deleted_files)} BAM files, saved ~{saved_gb:.2f} GB")


def update_version_in_config():
    """
    Starting with v0.9, the spacemake version is recorded inside the config.yaml
    as 'spacemake_version: '.

    For migration of earlier 0.8x -> 0.9, this entry has to be added.
    For future mirgations, the record is updated. 
    """
    # save original timestamps
    try:
        original_stat = os.stat("config.yaml")
    except FileNotFoundError:
        print("config.yaml not found.")
        return
    
    with open("config.yaml", "r") as f:
        config = yaml.safe_load(f)

    if "spacemake_version" not in config:
        config["spacemake_version"] = __version__
        with open("config.yaml", "w") as f:
            yaml.dump(config, f)

        # restore timestamps
        os.utime("config.yaml", (original_stat.st_atime, original_stat.st_mtime))

    else:
        # placeholder for future migrations
        return
