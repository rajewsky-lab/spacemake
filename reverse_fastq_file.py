#########
# about #
#########
__version__ = '0.1.0'
__author__ = ['Nikos Karaiskos', 'Tamas Ryszard Sztanka-Toth']
__licence__ = 'GPL'
__email__ = ['nikolaos.karaiskos@mdc-berlin.de', 'tamasryszard.sztanka-toth@mdc-berlin.de']

###########
# imports #
###########
import gzip
import subprocess
import os
from tqdm import tqdm

#############
# functions #
#############


def estimate_num_lines(file_path):
    # read the size of the file
    full_size = os.stat(file_path).st_size
    
    # need to catch errors for small file sizes, leads to wrong estimation

    # extract the first 1,000,000 reads to estimate file size
    subprocess.call("zcat " + file_path + " | head -4000000 > " + snakemake.params.tmp_file_pattern + ".fastq", shell=True)

    # zip the small file
    subprocess.call("pigz " + snakemake.params.tmp_file_pattern + ".fastq", shell=True)

    # read the size of the small file
    small_file_size = os.stat(snakemake.params.tmp_file_pattern + ".fastq.gz").st_size

    # remove the intermediate file
    subprocess.call("rm " + snakemake.params.tmp_file_pattern + ".fastq.gz", shell=True)

    # return the estimated number of lines of the full size
    return int(full_size * 4000000 / small_file_size)

def reverse_fastq_file(input_fq, output_fq):
    print ('reversing file', input_fq)

    # reverse the fastq file
    idx = 1
    with gzip.open(input_fq, 'rt') as fi, gzip.open(output_fq, 'wt') as fo:
        for line in tqdm(fi, total=estimate_num_lines(input_fq)):
            if (idx % 4 == 2) or (idx % 4 == 0):
                line = line[:20] # read only UMI+barcode
                line = line.strip('\n')
                fo.write(line[::-1] + '\n')
            else:
                fo.write(line)
            idx += 1
    print('reversing finished')

# create the directory for the reverse sequences
subprocess.call("mkdir -p " + os.path.dirname(str(snakemake.output)), shell=True)

# reverse the sequence
reverse_fastq_file(str(snakemake.input), str(snakemake.output))
