import os
from collections import defaultdict
import pysam
import numpy as np
import pandas as pd
import argparse
import spacemake.util as util
import logging

def parse_args():
    parser = util.make_minimal_parser(prog="bamstats.py")
    parser.add_argument("bamfile", help="input bam file")
    parser.add_argument("--out-stats", default="{bampath}/{bamname}.bc_stats.tsv", help="barcode and UMI statistics output")
    parser.add_argument("--out-length", default="{bampath}/{bamname}.readlen.tsv", help="read-length distribution output file")
    args = parser.parse_args()

    return args

logger = logging.getLogger("spacemake.bamstats")

def scan(fname, skim=0):
    UMI_counter = defaultdict(set)
    read_counter = defaultdict(int)
    read_len_hist = defaultdict(int)

    # print(fname)
    sam = util.quiet_bam_open(fname, "rb", check_sq=False, threads=2)
    last_qname = ""
    for read in util.timed_loop(sam.fetch(until_eof=True), logger, skim=skim):
        # count each read only once regardless of how many times it aligns
        if read.query_name == last_qname:
            continue

        last_qname = read.query_name

        CB = read.get_tag("CB")
        MI = read.get_tag("MI")

        UMI_counter[CB].add(MI)
        read_counter[CB] += 1

        seq = read.query_sequence
        read_len_hist[len(seq)] += 1

    # print("done processing")
    import numpy as np

    CBs = []
    UMI_counts = []
    read_counts = []
    for cb, umis in UMI_counter.items():
        CBs.append(cb)
        UMI_counts.append(len(umis))
        read_counts.append(read_counter[cb])

    CBs = np.array(CBs)
    UMI_counts = np.array(UMI_counts)
    read_counts = np.array(read_counts)

    # I = UMI_counts.argsort()[::-1]
    # for n_umis, cb in zip(UMI_counts[I], CBs[I]):
    #     print(f"{cb}\t{n_umis}")

    if len(read_len_hist):
        Lmax = np.array(list(read_len_hist.keys())).max()
    else:
        Lmax = 0

    Lhist = np.array([read_len_hist[x] for x in np.arange(Lmax + 1)])
    return CBs, UMI_counts, read_counts, Lhist


def main(args):
    logger = util.setup_logging(args, "spacemake.scripts.bamstats")
    bamname = os.path.basename(args.bamfile).split('.bam')[0]
    bampath = os.path.dirname(args.bamfile)
    sample = args.sample.format(bamname=bamname)

    CBs, UMI_counts, read_counts, Lhist = scan(args.bamfile)
    fname = util.ensure_path(args.out_stats.format(bamname=bamname, bampath=bampath))
    df = pd.DataFrame(dict(sample=sample, bam=bamname, CB=CBs, UMI=UMI_counts, reads=read_counts))
    df.to_csv(fname, sep="\t", index=None)

    fname = util.ensure_path(args.out_length.format(bamname=bamname, bampath=bampath))
    df = pd.DataFrame(dict(sample=sample, readlen=np.arange(len(Lhist)), count=Lhist))
    df.to_csv(fname, sep="\t", index=None)


if __name__ == "__main__":
    args = parse_args()
    main(args)
