#!/bin/bash

# extract the read categories from a dropseq-output bam file

# take each read, split the field starting with gf:Z:, and then take whatever is there, and finally count the occurrences. in the end sort by last column

sambamba view $1 | awk 'BEGIN {split("", count)} $5==255{for(i=1; i<=NF;++i){if($i ~ "^gf:Z:"){split($i, rtype, ":")}} if(rtype[3] in count){count[rtype[3]]++} else {count[rtype[3]] = 1} } END {for (key in count) print count[key]" "key}' | sort -rV
