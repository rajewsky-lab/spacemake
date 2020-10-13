import pysam
import difflib
import pandas as pd

optical_primer = 'GAATCACGATACGTACACCA' 

# get kmer lengths from length of primer until 4-mers
kmer_lengths = range(20, 3, -1)

# generate all posibble kmers as a list
kmer_list = [optical_primer[i:i+k] for k in kmer_lengths for i in range(len(optical_primer) - k + 1)]

# create a dataframe to store the counts
summary_table = pd.DataFrame({'count': 0}, index = list(kmer_lengths) + [0])
summary_table.index.rename('kmer_length', inplace=True)

def get_R1(read):
    return read.get_tag('XM')[::-1] + read.get_tag('XC')[::-1]

with pysam.AlignmentFile(snakemake.input[0], 'rb') as in_bam:
    with pysam.AlignmentFile(snakemake.output['tagged_bam'], 'wb', header = in_bam.header) as out_bam:
        optical_primer_len = len(optical_primer)

        # set the length of R1
        first_read = list(in_bam.head(1))[0]
        R1_len = len(get_R1(first_read))

        # for each read we need to look for the highest k-mer present
        # k-mers are ordered by their size, so first we find 20-mers etc
        for read in in_bam:
            R1 = get_R1(read)
            
            matcher = difflib.SequenceMatcher(None, optical_primer, R1)
            
            pos_optical_primer, pos_R1, kmer_len = matcher.find_longest_match(0, optical_primer_len, 0, R1_len)

            # in case kmer is shorter than 4nt, we simply output 0 for the match
            if kmer_len < 4:
                read.set_tag('XL', 0, 'i')
                read.set_tag('XK', '', 'Z')
                read.set_tag('XP', 0, 'i')

                # add one count to the 0 kmer length in the summary table
                summary_table.loc[0, 'count'] = summary_table.loc[0, 'count'] + 1
            # else, we have a longer match than 3-mer
            else:
                read.set_tag('XL', kmer_len, 'i')
                read.set_tag('XK', optical_primer[pos_optical_primer:pos_optical_primer+kmer_len], 'Z')
                read.set_tag('XP', pos_R1, 'i')

                summary_table.loc[kmer_len, 'count'] = summary_table.loc[kmer_len, 'count'] + 1

            # print the read to the output
            out_bam.write(read)

            # delete matcher object to save memory
            del matcher

summary_table.to_csv(snakemake.output['summary'])
