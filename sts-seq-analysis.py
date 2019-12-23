#########
# about #
#########
__version__ = '0.1'
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

#############
# functions #
#############

def is_gzip_file(filename):
    try:
        gzip.GzipFile(filename).peek(1)
        return True
    except OSError:
        return False

def open_file(filename):
    if (is_gzip_file(filename)):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'rt')

def rename_fastq_files(folder):
    """Rename all fastq files for meaningful naming: removes all SX_RY_001
    from the file names.
    folder -- the folder containing the .fastq.gz files."""

    filenames = []
    # r=root, d=directories, f = files
    for r, d, f in os.walk(folder):
        for file in f:
            if '.fastq.gz' in file:
                filenames.append(os.path.join(r, file))

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

def reverse_fastq_file(file):
    print ('reversing file', file)
    # first identify how many bases were sequenced
    with open_file(file) as fi:
        for i, line in enumerate(fi):
            if i == 1:
                read1_length = len(line.strip('\n'))
                break

    # reverse the fastq file
    idx = 1
    with open_file(file) as fi, open(re.sub('_1.fastq.gz', '_reversed_1.fastq', file), 'w') as fo:
        for line in fi:
            if idx % 4000000 == 0:
                print ('processed', int(idx/4), 'reads')
            if (idx % 4 == 2) or (idx % 4 == 4):
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
    # it will be passed as an argument to the file (or through a yaml file)
    folder = '/data/rajewsky/sequencing/slideSeq/191220_sts014_1-4_mouse_brain/'

    rename_fastq_files(folder)


    # get filenames of new fastq files to reverse their read one
    filenames = []
    # r=root, d=directories, f = files
    for r, d, f in os.walk(folder):
        for file in f:
            if '.fastq.gz' in file:
                filenames.append(os.path.join(r, file))

    # keep only read 1 files
    filenames = [x for x in filenames if '_1.fastq.gz' in x]
    
    print ('\n')
    print ('Reversing sequences of', len(filenames), 'files')

    for file in filenames:
        reverse_fastq_file(file)
        print ('\n')

    # get filenames of reversed fastq files
    filenames = []
    for r, d, f in os.walk(folder):
        for file in f:
            if '_reversed' in file and 'fastq' in file:
                filenames.append(os.path.join(r, file))

    # At this stage Read1 files are reversed, need to zip them
    print ('Zipping', len(filenames), 'files')
    print ('\n')

    for file in filenames:
        print ('zipping', file)
        subprocess.call("pigz "+file, shell=True)





    




