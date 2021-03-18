#!/usr/bin/env python3
__version__ = "0.9"
__author__ = ["Marvin Jens", ]
__license__ = "GPL"
__email__ = ['marvin.jens@mdc-berlin.de', ]

import pysam
import argparse
import os
import sys
import logging

bam_star = pysam.AlignmentFile('/dev/stdin', 'rb')
bam_tagged = pysam.AlignmentFile(sys.argv[1], 'rb', check_sq=False)

def print_header(header):
    for k, v in sorted(header.items()):
        if type(v) == dict:
            vstr = " ".join([f"{x}:{y}" for x, y in sorted(v.items())])
            print(f"@{k}:\t{vstr}")
        elif type(v) == list:
            for row in v:
                if type(row) == dict:
                    vstr = " ".join([f"{x}:{y}" for x, y in sorted(row.items())])
                else:
                    vstr = str(row)
                print(f"@{k}:\t{vstr}")
        else:
            print(f"@{k}:\t{v}")


def merge_headers(orig, star):
    merged = dict(orig)
    merged['PG'].extend(star['PG'])
    merged['SQ'] = star['SQ']
    # merged['HD']['SO'] = star['HD']['SO']  # sorted by

    return merged

star_header = bam_star.header.to_dict()
tagged_header = bam_tagged.header.to_dict()
merged_header = merge_headers(tagged_header, star_header)
# print(f"STAR header")
# print_header(star_header)

# print(f"original header")
# print_header(tagged_header)

# print("merged header")
# print_header(merged_header)

# copy input to output on stdout, just with the new header
sam_out = pysam.AlignmentFile('/dev/stdout', 'wb', header=merged_header)
for aln in bam_star.fetch(until_eof=True):
    sam_out.write(aln)
