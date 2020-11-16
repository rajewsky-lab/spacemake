import gzip
import subprocess
import argparse
import os

## create argparser
parser = argparse.ArgumentParser(description='get umi from R2 and reverse R1')
parser.add_argument('--in_R1', type=str, required=True, help = 'input R1 fastq.gz file')
parser.add_argument('--in_R2', type=str, required=True, help = 'input R2 fastq.gz file')
parser.add_argument('--out_R1', type=str, required=True, help = 'output R1 fastq.gz file')
args = parser.parse_args()

def create_r1_umi_from_r2(fastq_in_r1, fastq_in_r2, fastq_out_r1):
    with gzip.open(fastq_in_r1, 'rt') as in_1, gzip.open(fastq_in_r2, 'rt') as in_2, gzip.open(fastq_out_r1, 'wt') as out_1:

        it = 0

        for line_1, line_2 in zip(in_1, in_2):
            if it % 2 == 1:
                # first 12 nt of R1 are cell barcode, first 8nt of R2 are the umi
                # but we need to reverse the cellbc
                line_1 = line_1[:12][::-1] + line_2[:8] + '\n'
            _ = out_1.write(line_1)
            it = it + 1

## call the function
subprocess.call('mkdir -p ' + os.path.dirname(str(args.out_R1)), shell=True)

print('getting UMI from %s, and reversing R1 for %s' % (args.in_R2, args.in_R1))

create_r1_umi_from_r2(args.in_R1, args.in_R2, args.out_R1) 
