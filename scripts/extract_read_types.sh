#!/bin/bash

# extract the read categories from a dropseq-output bam file

# take each read, split the field starting with gf:Z:, and then take whatever is there, and finally count the occurrences. in the end sort by last column

# select reads where mapping quality is bigger than 10 (255 means that it is certainly unique, 10 means that it is most probably unique). split the read type field. then count the occurrences of read types.
sambamba view $1 | awk 'BEGIN {split("", count)} $5>10{for(i=1; i<=NF;++i){if($i ~ "^gf:Z:"){split($i, rtype, ":")}} if(rtype[3] in count){count[rtype[3]]++} else {count[rtype[3]] = 1} } END {for (key in count) print count[key]" "key}' | \
    # for each row now we have two columns: $1 = the number of occurrences. $2 = the feature names the read maps to, separated by commas. Whenever each feature is the same per row, we collapse this information into a single feature
    # ie UTR,UTR,UTR becomes UTR for example, or CODING,CODING becomes CODING. For other cases, like CODING,INTERGENIC, we simply assign AMBIGUOUS.
    awk '{split($2, types, ","); types_l = length(types); eq_count = 0; first_field = types[1]; for (e in types){if(types[e] == first_field){eq_count++}}; if(eq_count == types_l) {print $1" "first_field} else {print $1" AMBIGUOUS"}}' | \
    # then we simply sum up the reads across features and output them in a table
    awk 'BEGIN { split("", out)} {if($2 in out) {out[$2]= out[$2] + $1} else {out[$2]=$1}} END { for (key in out){ print key, out[key]}}'

