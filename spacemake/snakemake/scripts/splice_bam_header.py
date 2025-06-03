#!/usr/bin/env python3
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


def read_header(source):
    hd = {
        # "HD": {},
        "SQ": [],
        "PG": [],
        "RG": [],
    }
    for line in open(source):
        if not line.startswith("@"):
            continue

        parts = line[1:].split("\t")
        top = parts[0]
        rec = {}
        for entry in parts[1:]:
            key, value = entry.rstrip().split(":", maxsplit=1)
            rec[key] = value

        if top in ["SQ", "PG", "RG"]:
            hd[top].append(rec)
        else:
            hd[top] = rec

    return hd


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


def sort_PG_list(pg_list):
    root = None
    pp = {}
    for pg in pg_list:
        if not "PP" in pg:
            # no Previous Program: we have the first entry
            root = pg
        else:
            pp[pg["PP"]] = pg

    l = [root]
    pp_id = root["ID"]
    for i in range(len(pg_list) - 1):
        pg = pp[pp_id]
        l.append(pg)
        pp_id = pg["ID"]

    return l


def unique_IDs(pg_list1, pg_list2):
    from collections import defaultdict

    id_renames = {}
    ids1 = set([pg["ID"] for pg in pg_list1])
    # print(f"pre-existing IDs: {ids1}")
    for pg in pg_list2:
        id2 = pg["ID"]
        id2_ = str(id2)
        n = 0
        while id2_ in ids1:
            # print(f"while enter: {id2_}")
            try:
                if "." in id2_:
                    n = int(id2_.split(".")[-1])
                    id2_ = id2_.rsplit(".", maxsplit=1)[0]
            except ValueError:
                pass

            n += 1
            id2_ = f"{id2_}.{n}"

        # print(f"end result: {id2} -> {id2_}")
        if id2_ != id2:
            ids1.add(id2_)
            id_renames[id2] = id2_

    for pg in pg_list2:
        pg["ID"] = id_renames.get(pg["ID"], pg["ID"])
        if "PP" in pg:
            pg["PP"] = id_renames.get(pg["PP"], pg["PP"])

    # print("suggest the following renames in pg_list2:")
    # print(id_renames)

    return pg_list2


test_pg_list1 = [
    {"ID": "fastq_to_uBAM"},
    {"ID": "bowtie2", "PP": "fastq_to_uBAM"},
]

test_pg_list2 = [
    {"ID": "bowtie2"},
]


def merge_PG_lists(orig, other):
    first = sort_PG_list(orig.copy())
    last = sort_PG_list(other.copy())

    # print("about to merge PG lists")
    # for pg in first:
    #     print("first:", pg)

    # for pg in last:
    #     print("last:", pg)

    last = unique_IDs(
        first, last
    )  # rename re-occuring IDs to make them unique if needed
    last[0]["PP"] = first[-1]["ID"]  # link the two histories
    return first + last


def merge_headers(orig, other, enforce_RG=True):
    merged = dict(orig)
    # start from the original, including VN and RG entries...
    # connect the processing chains:
    # most recent program should be on top
    #  Previous Program (PP) of the first new output was the last Program Name (PN) in the original uBAM
    # other["PG"][-1]["PP"] = orig["PG"][0]["ID"]
    merged["PG"] = merge_PG_lists(
        orig["PG"], other["PG"]
    )  # unique_IDs(merged["PG"] + other["PG"])

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
    parser.add_argument("--in-ubam", help="unmapped dropseq tagged bam")
    parser.add_argument(
        "--out-bam",
        help="fixed output bam (default=/dev/stdout)",
        default="/dev/stdout",
    )
    parser.add_argument("--out-mode", help="mode for output (default=b0)", default="b0")
    parser.add_argument("--in-header1", help="first header")
    parser.add_argument("--in-header2", help="second header")

    args = parser.parse_args()

    if args.in_bam and args.in_ubam:
        mbam = pysam.AlignmentFile(args.in_bam, "r")
        ubam = pysam.AlignmentFile(args.in_ubam, "rc", check_sq=False)
        header1 = ubam.header.to_dict()
        header2 = mbam.header.to_dict()

    if args.in_header1 and args.in_header2:
        # overwrite headers if so desired
        header1 = read_header(args.in_header1)
        header2 = read_header(args.in_header2)

    # print_header(header1)
    # print_header(header2)
    merged_header = merge_headers(header1, header2)
    # print_header(merged_header)
    # print(f"mapped BAM header")
    # print_header(header1)

    # print(f"original uBAM header")
    # print_header(ubam_header)

    # print("merged header")
    # print_header(merged_header)

    # copy input to output, just with the new header
    if args.in_bam:
        bam_out = pysam.AlignmentFile(
            args.out_bam, f"w{args.out_mode}", header=merged_header
        )
        for aln in mbam.fetch(until_eof=True):
            bam_out.write(aln)
