import os
import subprocess
import time


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

    return list(zip(bam_files, are_symlinks))


def sync_timestamps(original_file, new_file):
    """
    Sync the timestamps (access and modification time) of new_file with those of original_file.
    
    Args:
        original_file (str): Path to the file whose timestamps will be copied.
        new_file (str): Path to the file that will have its timestamps updated.
    """
    try:
        # Get the access time and modification time from original_file
        source_times = os.stat(original_file)

        # Set the same access and modification time for new_file
        os.utime(new_file, (source_times.st_atime, source_times.st_mtime))

        print(f"File timestamps of {new_file} set to match {original_file}.")
    except FileNotFoundError:
        print(f"Error: One or both of the files '{original_file}' or '{new_file}' do not exist.")
    except Exception as e:
        print(f"An error occurred: {e}")


def convert_bam_to_cram(ref_sequence, project_folder, threads=4):
    bam_files = find_bam_files(project_folder)

    for idx in range(len(bam_files)):
        bam_filename, bam_file_is_symlink = bam_files[idx]
        bam_filename_prefix = bam_filename.rsplit('.', 1)[0]
        cram_filename = bam_filename_prefix + '.cram'

        if os.path.exists(os.path.join(project_folder, cram_filename)):
            print('CRAM file', cram_filename, 'already exists. Skipping conversion.')
            continue

        # TODO: change ref sequence for genome, rRNA, phiX, custom?
        
        if bam_file_is_symlink:
            # TODO: fix timestamp for symlink 
            true_bam_filename = os.readlink(os.path.join(project_folder, bam_filename))
            true_bam_filename_prefix = true_bam_filename.rsplit('.', 1)[0]
            os.symlink(true_bam_filename_prefix + '.cram',
                       os.path.join(project_folder, bam_filename_prefix + '.cram'))
        else:
            print('Converting', bam_filename, 'to', cram_filename, 
            '...', time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
            subprocess.run(
                [
                    "samtools", "view",
                    "-T", ref_sequence,
                    "-C",
                    "--threads", str(threads),
                    "-o", os.path.join(project_folder, cram_filename),
                    os.path.join(project_folder, bam_filename)
                ]
            )

        sync_timestamps(os.path.join(project_folder, bam_filename),
                        os.path.join(project_folder, cram_filename))


def remove_files():
    # - BAM files (only if CRAMs are present)
    bam_files = find_bam_files(project_folder)

    # - unaligned.bam

    # remove tiles

