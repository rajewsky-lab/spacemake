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
    are_symlinks = [os.path.islink(bam_file) for bam_file in bam_files]

    return list(zip(bam_files, are_symlinks))

def convert_bam_to_cram(ref_sequence, project_folder, threads=4):
    bam_files = find_bam_files(project_folder)

    for idx in range(len(bam_files)):
        bam_filename, bam_file_is_symlink = bam_files[idx]
        bam_filename_prefix = bam_filename.rsplit('.', 1)[0]
        cram_filename = bam_filename_prefix + ".cram"

        # TODO: change ref sequence for genome, rRNA, phiX, custom?
        
        if bam_file_is_symlink:
            # TODO: deal with this
            continue
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

    # TODO: transfer timestamp

    return

def remove_files():
    # - BAM files (if CRAMs are present)

    # - unaligned.bam

    # remove tiles

    return

