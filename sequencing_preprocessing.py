#########
# about #
#########
__version__ = '0.1.3'
__author__ = ['Nikos Karaiskos']
__licence__ = 'GPL'
__email__ = ['nikolaos.karaiskos@mdc-berlin.de']

###########
# imports #
###########
import time
import os
import numpy as np
import re
import gzip
import subprocess
from tqdm import tqdm
import mmap
import sys

#############
# functions #
#############

def is_gzip_file(file_path):
    """Checks if file is gzipped or not."""
    try:
        gzip.GzipFile(file_path).peek(1)
        return True
    except OSError:
        return False

def open_file(file_path):
    """Extension of the standard open function in python that opens a file
    for reading regardless of being zipped or not."""
    if (is_gzip_file(file_path)):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'rt')

def estimate_num_lines(file_path):
    # read the size of the file
    full_size = os.stat(file_path).st_size
    
    # need to catch errors for small file sizes, leads to wrong estimation

    # extract the first 1,000,000 reads to estimate file size
    subprocess.call("zcat " + file_path + " | head -4000000 > temp_small.fastq", 
        shell=True)

    # zip the small file
    subprocess.call("pigz temp_small.fastq", shell=True)

    # read the size of the small file
    small_file_size = os.stat('temp_small.fastq.gz').st_size

    # remove the intermediate file
    subprocess.call("rm temp_small.fastq.gz", shell=True)

    # return the estimated number of lines of the full size
    return int(full_size * 4000000 / small_file_size)

def get_filenames(folder, pattern_list):
    """Get the names of files in a given folder that contain a list of 
    strings. It descends onto subfolders and searches recursively.
    folder -- the folder containing the data
    list   -- a list of strings to look for in the filenames"""
    filenames = []
    # r=root, d=directories, f = files
    for r, d, f in os.walk(folder):
        for file in f:
            if sum([pattern in file for pattern in pattern_list]) == len(pattern_list):
                filenames.append(os.path.join(r, file))
    return filenames

def rename_fastq_files(folder):
    """Rename all fastq files for meaningful naming: removes all SX_RY_001
    from the file names.
    folder -- the folder containing the .fastq.gz files."""

    filenames = get_filenames(folder, ['.fastq.gz'])

    # these are the new filenames
    new_filenames = []
    for file in filenames:
        if 'R1_001.fastq.gz' in file:
            new_filenames.append(re.sub('_S[0-9]+_R[0-9]_[0-9]+.fastq.gz',
                '_1.fastq.gz', file))
        elif 'R2_001.fastq.gz' in file:
            new_filenames.append(re.sub('_S[0-9]+_R[0-9]_[0-9]+.fastq.gz',
                '_2.fastq.gz', file))

    # rename the files and continue the analysis
    for file in range(len(filenames)):
        os.rename(filenames[file], new_filenames[file])        

def reverse_fastq_file(file_path):
    print ('reversing file', file_path)

    # reverse the fastq file
    idx = 1
    with open_file(file_path) as fi, open(re.sub('_1.fastq.gz', '_reversed_1.fastq', file), 'w') as fo:
        for line in tqdm(fi, total=estimate_num_lines(file_path)):
            if (idx % 4 == 2) or (idx % 4 == 0):
                line = line[:20] # read only UMI+barcode
                line = line.strip('\n')
                fo.write(line[::-1] + '\n')
            else:
                fo.write(line)
            idx += 1


########
# main #
########

if __name__ == '__main__':

    # this is the folder where the demultiplexed data is.
    folder = sys.argv[1]
    rename_fastq_files(folder)

    # get filenames of new fastq files to reverse their read1
    filenames = get_filenames(folder, ['.fastq.gz'])

    # split into read1 and read2 files
    filenames_read1 = [x for x in filenames if '_1.fastq.gz' in x]
    filenames_read2 = [x for x in filenames if '_2.fastq.gz' in x]

    print ('Reversing sequences of', len(filenames_read1), 'files')

    for file in filenames_read1:
        reverse_fastq_file(file)

    # get filenames of reversed fastq files
    filenames = get_filenames(folder, ['_reversed', 'fastq'])

    print ('\n')

    # zip read1 files. Currently requires pigz to be installed
    print ('Zipping', len(filenames), 'files')

    for file in filenames:
        print ('zipping', file)
        subprocess.call("pigz " + file, shell=True)

    print ('\n')

    # create symbolic links for the Read2 files
    for file in filenames_read2:
        print ('creating symbolic link for', file)
        subprocess.call("ln -s " + file + ' ' +
            re.sub('_2.fastq.gz', '_reversed_2.fastq.gz', file), shell=True)
