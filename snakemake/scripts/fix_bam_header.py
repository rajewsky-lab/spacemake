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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fix .bam header of the STAR mapped output .bam')

    parser.add_argument('--in-bam-star', help='mapped star bam input')
    parser.add_argument('--in-bam-tagged', help='unmapped dropseq tagged bam')
    parser.add_argument('--out-bam', help='output bam')

    args = parser.parse_args()

    bam_star = pysam.AlignmentFile(args.in_bam_star, 'rb')
    bam_tagged = pysam.AlignmentFile(args.in_bam_tagged, 'rb', check_sq=False)


    star_header = bam_star.header.to_dict()
    tagged_header = bam_tagged.header.to_dict()
    merged_header = merge_headers(tagged_header, star_header)
    # print(f"STAR header")
    # print_header(star_header)

    # print(f"original header")
    # print_header(tagged_header)

    # print("merged header")
    # print_header(merged_header)

    # copy input to output, just with the new header
    bam_out = pysam.AlignmentFile(args.out_bam,'wb', header=merged_header)
    for aln in bam_star.fetch(until_eof=True):
        bam_out.write(aln)
