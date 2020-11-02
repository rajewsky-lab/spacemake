import gzip
import subprocess

def create_r1_umi_from_r2(fastq_in_r1, fastq_in_r2, fastq_out_r2):
    with gzip.open(fastq_in_r1, 'rt') as in_1, gzip.open(fastq_in_r2, 'rt') as in_2, gzip.open(fastq_out_r1, 'wt') as out_1:
        for line_1, line_2 in zip(in_1, in_2):
            if it % 2 == 1:
                # first 12 nt of R1 are cell barcode, first 8nt of R2 are the umi
                line_1 = line_1[:12] + line_2[:8] + '\n'
            _ = out_1.write(line_1)
            it = it + 1

subprocess.call('mkdir -p ' + os.path.dirname(str(snakemake.output)), shell=True)

create_r1_umi_from_r2(snakemake.input['R1'], snakemake.input['R2'], snakemake.output) 
