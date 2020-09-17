import pysam
import pandas as pd

optical_primer = 'GAATCACGATACGTACACCA' 

# get kmer lengths from length of primer until 4-mers
kmer_lengths = range(20, 3, -1)

# generate all posibble kmers as a list
kmer_list = [optical_primer[i:i+k] for k in kmer_lengths for i in range(len(optical_primer) - k + 1)]

# create a length hash table to speedup the code
kmers = {x: len(x) for x in kmer_list}

# create a dataframe to store the counts
summary_table = pd.DataFrame({'count': 0}, index = list(kmer_lengths) + [0])
summary_table.index.rename('kmer_length', inplace=True)

def get_R1(read):
    return read.get_tag('XM')[::-1] + read.get_tag('XC')[::-1]

with pysam.AlignmentFile('../scratch/test.bam', 'rb') as in_bam:
    with pysam.AlignmentFile('out/kmer.bam', 'wb', header = in_bam.header) as out_bam:
        # for each read we need to look for the highest k-mer present
        # k-mers are ordered by their size, so first we find 20-mers etc
        for read in in_bam:
            R1 = get_R1(read)
            
            for kmer in kmers.keys():
                kmer_pos = R1.find(kmer)

                # if kmer is found, it returns the position of the match
                if kmer_pos > -1:
                    kmer_found = True
                    # set the correct tags
                    read.set_tag('XL', kmers[kmer], 'i')
                    read.set_tag('XK', kmer, 'Z')
                    read.set_tag('XP', kmer_pos, 'i')
                    
                    # update the summary table
                    summary_table.loc[kmers[kmer], 'count'] = summary_table.loc[kmers[kmer], 'count'] + 1

                    # stop the iteration here
                    break

            # in case we didn't find any kmer, which means that all find commands returned -1
            if kmer_pos == -1:
                read.set_tag('XL', 0, 'i')
                read.set_tag('XK', '', 'Z')
                read.set_tag('XP', 0, 'i')

                # add one count to the 0 kmer length in the summary table
                summary_table.loc[0, 'count'] = summary_table.loc[0, 'count'] + 1

            # print the read to the output
            out_bam.write(read)

'asd'[::-1]

summary_table.to_csv('summary.csv')
