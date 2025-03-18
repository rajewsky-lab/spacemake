import pandas as pd
import argparse
import re

parser = argparse.ArgumentParser(
    description="Split .bam file to .sam files by mapped reads strand orientation"
)
parser.add_argument("file_in", metavar="in", type=str)
parser.add_argument("--prefix")

args = parser.parse_args()

prefix = args.prefix

read_type_num = {"INTERGENIC": 0, "INTRONIC": 0, "CODING": 0, "UTR": 0, "AMB": 0}

strand_type_num = {
    "minus_minus": 0,
    "minus_plus": 0,
    "plus_plus": 0,
    "plus_minus": 0,
    "plus_AMB": 0,
    "minus_AMB": 0,
}

# out_file_names = {x: prefix + x + ".sam" for x in strand_type_num.keys()}

# out_files = {x: open(out_file_names[x], "w") for x in out_file_names.keys()}


def return_collapsed(it):
    # set has exactly 1 element, meaning that all elements are the same in the list
    if len(set(it)) == 1:
        return it[0]
    else:
        return "AMB"


def get_tag_unique(line, tag="gs"):
    if m := re.search(f"{tag}\:Z\:(\S+)", line):
        gs = set(m.groups()[0].split(","))
        if len(gs) == 1:
            return gs.pop()


with open(args.file_in, "r") as fi:
    for line in fi:
        # if line is header line
        if line.startswith("@"):
            # for f in out_files.values():
            #     f.write(line)

            # # go to next iteration
            continue

        line_stripped = line.strip()

        elements = line_stripped.split()

        last = elements[-1]

        read_overlaps_gene = False

        # set gene strand
        if gene_strand := get_tag_unique(line_stripped, tag="gs"):
            # if last element begins with gs, this means that read overlaps a gene (fw or rv strands)
            read_overlaps_gene = True
            if gene_strand == "-":
                gene_strand = "minus"
            elif gene_strand == "+":
                gene_strand = "plus"

        else:
            gene_strand = "AMB"

        # print(f"gene_strand={gene_strand}")
        # set read strand
        if int(elements[1]) & 16:
            read_strand = "minus"
        else:
            read_strand = "plus"

        # get read type
        if read_overlaps_gene:
            read_type = get_tag_unique(line_stripped, tag="XF")
            # print(f"read_type={read_type}")
            # return_collapsed(elements[-3].split(":")[-1].split(","))
        else:
            # if read do not overlap a gene, it is clearly intergenic
            read_type = "INTERGENIC"

        read_type_num[read_type] = read_type_num[read_type] + 1

        strand_type = read_strand + "_" + gene_strand

        strand_type_num[strand_type] = strand_type_num[strand_type] + 1

        # # print the read to the correct split file, depending on strand orientation
        # out_files[strand_type].write(line)

with open(prefix + "read_type_num.txt", "w") as fo:
    for key, value in read_type_num.items():
        fo.write("%s %s\n" % (key, value))

with open(prefix + "strand_type_num.txt", "w") as fo:
    for key, value in strand_type_num.items():
        fo.write("%s %s\n" % (key, value))
