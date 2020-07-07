#!/bin/bash

# extract the read categories from a dropseq-output bam file

# take each read, split the field starting with gf:Z:, and then take whatever is there, and finally count the occurrences. in the end sort by last column

# take only 10 million reads, for faster computation
MAX_NUM_READS=10000000

sambamba view $1 | head -n $MAX_NUM_READS | grep -e 'gs:Z' |\
    # print the 2nd field (read strand) and the last field (gene strands)
    # extract the strand info from the last field
    awk '{split($NF, gsz, ":")}{print $2" "gsz[3]}' | \
    # collapse the strand information. it is currently on the second field
    awk '{  split($2, strands, ",");
            strands_l = length(strands);
            first_strand = strands[1];
            eq_count = 0;
            for(s in strands){
                if(strands[s] == first_strand){eq_count++}
            };
            if(eq_count == strands_l){print $1" "first_strand} else {print $1" amb"}
        }' |\
    sort | uniq -c | sort -nr
