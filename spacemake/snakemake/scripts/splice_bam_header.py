#!/usr/bin/env python3
from spacemake.contrib import __author__, __email__, __version__
import pysam
import argparse
import os
import sys
import logging


def read_header(source):
    hd = {}
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
            if not top in hd:
                hd[top] = []
            hd[top].append(rec)
        else:
            hd[top] = rec

    return hd


def make_ID_unique(ids1, id2):
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

    return id2_


def repair_PG_lines(pg_list, fix=False):
    from collections import defaultdict

    issues = []

    no_PP = []
    w_PP = []

    ids_seen = set()
    ids_by_cmd = defaultdict(list)
    rename = {}

    for i, pg in enumerate(pg_list):
        cmd = pg["ID"].rsplit(".", maxsplit=1)[0]
        if fix:
            pg["ID"] = rename.get(pg["ID"], pg["ID"])

        if "PP" in pg:
            if fix:
                pg["PP"] = rename.get(pg["PP"], pg["PP"])

            w_PP.append(pg)
            if pg["ID"] == pg["PP"]:
                best_guess = ids_by_cmd[cmd][-1]
                issues.append(
                    f"entry ID={pg['ID']} references itself: PP={pg['PP']}. Best guess: correct PP={best_guess}"
                )
                if fix:
                    # correct entry is most likely previous run of same tools
                    pg["PP"] = best_guess

            if pg["PP"] not in ids_seen:
                issues.append(
                    f"entry ID={pg['ID']} references PP={pg['PP']}, which has not been defined yet. At a minimum, the ordering is off."
                )
        else:
            no_PP.append(pg)
            if i > 0:
                best_guess = pg_list[i - 1]["ID"]
                issues.append(
                    f"subsequent PG entry without PP tag ID={pg['ID']}. Have to guess that {no_PP[0]['ID']} is actual beginning. Best guess: correct PP={best_guess}"
                )

                if fix:
                    pg["PP"] = best_guess

        if pg["ID"] in ids_seen:
            issues.append(f"multiple use of ID {pg['ID']}.")
            if fix:
                new_ID = make_ID_unique(ids_seen, pg["ID"])
                rename[pg["ID"]] = new_ID
                pg["ID"] = new_ID

        ids_seen.add(pg["ID"])
        ids_by_cmd[cmd].append(pg["ID"])

    if fix:
        pg_list = sort_PG_list(pg_list)

    return issues, pg_list


def header_to_SAM(header):
    SAM = []
    for k, v in header.items():
        if type(v) == dict:
            vstr = "\t".join([f"{x}:{y}" for x, y in v.items()])
            SAM.append(f"@{k}\t{vstr}")
        elif type(v) == list:
            for row in v:
                if type(row) == dict:
                    keys = list(row.keys())
                    if k == "PG":
                        # override order to ensure CL is always last for max compatibility
                        keys = [x for x in ["ID", "PN", "VN", "PP", "CL"] if x in row]

                    vstr = "\t".join([f"{x}:{row[x]}" for x in keys])
                else:
                    vstr = str(row)
                SAM.append(f"@{k}\t{vstr}")
        else:
            SAM.append(f"@{k}\t{v}")

    return SAM


def print_header(header):
    print("\n".join(header_to_SAM(header)))


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


def fix_PP(repair_header, fix=False):
    """
    Due to a bug in previous versions, we may end up with ID:bowtie2.2 and PP:bowtie2.2
    in the same PG entry. This breaks samtools view (unless --no-PG is used). So, here we
    detect this very specific error and fix it.
    """
    import pysam

    bam = pysam.AlignmentFile(repair_header, "rc", check_sq=False)
    hdr = bam.header.to_dict()
    issues, pg_list = repair_PG_lines(hdr["PG"], fix=True)
    for iss in issues:
        print(f"Header issue detected: '{iss}'.")

    if len(issues) == 0:
        print("No issues detected.")
    else:
        issues_after_fix, pg_list = repair_PG_lines(hdr["PG"], fix=False)
        assert len(issues_after_fix) == 0
        print(f"Successfully generated a fixed header.")

        from tempfile import NamedTemporaryFile

        if fix:
            with NamedTemporaryFile(mode="tw") as tmp_hdr:
                tmp_hdr.write("\n".join(header_to_SAM(hdr)) + "\n")
                tmp_hdr.flush()
                # tmp_hdr.close()
                import pysam

                pysam.samtools.reheader("-P", "-i", tmp_hdr.name, args.repair_header)
                print(
                    f"Replaced original header of '{args.repair_header}' with fixed header."
                )


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
    parser.add_argument("--repair-header", default="")
    parser.add_argument("--check-header", default="")

    args = parser.parse_args()

    header1 = None
    header2 = None

    if args.repair_header:
        fix_PP(args.repair_header, fix=True)
    elif args.check_header:
        fix_PP(args.check_header, fix=False)
    else:
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
        if header1 and header2:
            merged_header = merge_headers(header1, header2)
        # print_header(merged_header)
        # print(f"mapped BAM header")
        # print_header(header1)

        # print(f"original uBAM header")
        # print_header(ubam_header)

        # print("merged header")
        # print_header(merged_header)

        # copy input to output, just with the new header
        if args.in_bam and header1 and header2:
            bam_out = pysam.AlignmentFile(
                args.out_bam, f"w{args.out_mode}", header=merged_header
            )
            for aln in mbam.fetch(until_eof=True):
                bam_out.write(aln)
