#!/usr/bin/env python3
__version__ = "0.9"
__author__ = [
    "Marvin Jens",
]
__license__ = "GPL"
__email__ = [
    "marvin.jens@mdc-berlin.de",
]

import pysam
import argparse
import os
import sys
import logging


def print_header(header):
    for k, v in sorted(header.items()):
        if type(v) == dict:
            vstr = " ".join([f"{x}:{y}" for x, y in sorted(v.items())[::-1]])
            print(f"@{k}:\t{vstr}")
        elif type(v) == list:
            for row in v:
                if type(row) == dict:
                    vstr = " ".join([f"{x}:{y}" for x, y in sorted(row.items())[::-1]])
                else:
                    vstr = str(row)
                print(f"@{k}:\t{vstr}")
        else:
            print(f"@{k}:\t{v}")


def unique_IDs(pg_list):
    from collections import defaultdict

    id_counts = defaultdict(int)

    # first, iterate over entire list and count how often each program ID is there.
    pp_list = [None]
    if len(pg_list) > 1:
        pp_list += pg_list[:-1]

    pg_new = []
    for pg, pp in zip(pg_list, pp_list):
        name = pg["ID"].split(".")[0]
        # edit in-place
        id_counts[name] += 1
        pg["ID"] = f"{name}.{id_counts[name]}"
        # id_counts[name] -= 1

        if pp:
            pname = pp["ID"].split(".")[0]
            pg["PP"] = f"{pname}.{id_counts[pname]}"

        pg_new.append(pg)

    return pg_new


def merge_headers(orig, other, enforce_RG=True):
    merged = dict(orig)
    # start from the original, including VN and RG entries...
    # connect the processing chains:
    # most recent program should be on top
    #  Previous Program (PP) of the first new output was the last Program Name (PN) in the original uBAM
    other["PG"][-1]["PP"] = orig["PG"][0]["ID"]
    merged["PG"] = unique_IDs(merged["PG"] + other["PG"])

    if "SO" in other["HD"]:
        merged["HD"]["SO"] = other["HD"]["SO"]  # keep sort-order

    # sequence identifiers should be absent from uBAM and at any rate are overwritten here
    merged["SQ"] = other["SQ"]
    if enforce_RG and not "RG" in merged or len(merged["RG"]) == 0:
        merged["RG"] = {"ID": "A", "SM": "NA"}
    # merged['HD']['SO'] = star['HD']['SO']  # sorted by

    return merged


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "STAR and bowtie2 create a new header from scratch and ignore everything upstream. "
            "This script fixes the .bam headers of such mapped output by splicing it together with "
            "the original uBAM header."
        )
    )

    parser.add_argument(
        "--in-bam",
        help="mapped star/bowtie2 bam input (default=/dev/stdin)",
        default="/dev/stdin",
    )
    parser.add_argument("--in-ubam", help="unmapped dropseq tagged bam", required=True)
    parser.add_argument(
        "--out-bam",
        help="fixed output bam (default=/dev/stdout)",
        default="/dev/stdout",
    )
    parser.add_argument("--out-mode", help="mode for output (default=b0)", default="b0")

    args = parser.parse_args()

    mbam = pysam.AlignmentFile(args.in_bam, "rb")
    ubam = pysam.AlignmentFile(args.in_ubam, "rb", check_sq=False)

    mapped_header = mbam.header.to_dict()
    ubam_header = ubam.header.to_dict()
    merged_header = merge_headers(ubam_header, mapped_header)
    # print(f"mapped BAM header")
    # print_header(mapped_header)

    # print(f"original uBAM header")
    # print_header(ubam_header)

    # print("merged header")
    # print_header(merged_header)

    # copy input to output, just with the new header
    bam_out = pysam.AlignmentFile(
        args.out_bam, f"w{args.out_mode}", header=merged_header
    )
    for aln in mbam.fetch(until_eof=True):
        bam_out.write(aln)
