import re

def parse_ribo_log(ribo_log_file):
    # before the log, there can be some perl warnings prepended. so we need to find the
    # first line outputed by bwa mem
    input_reads = 0
    aligned_reads = 0

    # ribo log summary line: first line of the summary
    first_line_regex = r'^\d+ reads; of these:$'
    first_line_found = False

    line_n = 0

    with open(ribo_log_file) as f:
        for line in f:
            stripped_line = line.strip()
            
            if stripped_line == 'no_rRNA_index':
                input_reads = -1
                aligned_reads = -1
                break

            if not first_line_found:
                if re.match(first_line_regex, stripped_line) is not None:
                    first_line_found = True
                    line_n = 0
                else:
                    # keep looking for first line
                    continue

            if line_n == 0:
                input_reads = input_reads + int(stripped_line.split(' ')[0])
            elif line_n == 3 or line_n == 4:
                aligned_reads = aligned_reads + int(stripped_line.split(' ')[0])
            # reset after the fifth line, this is needed if there are several ribolog files
            # appended one after the other. this is the case for merged samples
            elif line_n == 5:
                first_line_found = False

            line_n = line_n + 1
            
    
    if input_reads <= 0:
        return (None, None)
    else:
        return (aligned_reads, input_reads)


aligned_reads, input_reads = parse_ribo_log(snakemake.input[0])

with open(snakemake.output[0], 'w') as fo:
    fo.write(f'aligned_reads\t{aligned_reads}\n')
    fo.write(f'input_reads\t{input_reads}\n')
