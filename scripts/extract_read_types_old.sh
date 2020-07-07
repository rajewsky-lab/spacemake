#!/bin/bash

# old script for extracting read types

samtools view $1 |
awk '!/GE:Z:/ && $5 == "255" && match ($0, "XF:Z:") split(substr($0, RSTART+5), a, "\t") {print a[1]}' | \
awk 'BEGIN { split("INTRONIC INTERGENIC CODING UTR", keyword)
        for (i in keyword) count[keyword[i]]=0
    }
    /INTRONIC/  { count["INTRONIC"]++ }
    /INTERGENIC/  { count["INTERGENIC"]++ }
    /CODING/ {count["CODING"]++ }
    /UTR/ { count["UTR"]++ }
    END   {
        for (i in keyword) print keyword[i], count[keyword[i]]
    }'
