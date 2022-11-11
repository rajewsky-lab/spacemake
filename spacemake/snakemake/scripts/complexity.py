import pysam
import numpy as np
import os
import argparse
import logging
from collections import defaultdict
import pandas as pd
import spacemake.util as util
from spacemake.contrib import __license__, __version__, __author__, __email__


def subsample(reads, n):
    if n == len(reads):
        sample = reads
    else:
        f = float(n) / len(reads)
        rnd = np.random.random(len(reads))
        mask = rnd <= f
        sample = list(reads[mask])

    # print(sample)
    n_uniq = len(set(sample))
    return len(sample), n_uniq


def get_raw_reads(sample):
    return int(open(f"counts/{sample}.rawcount.txt").read())


def parse_args():
    parser = argparse.ArgumentParser(
        description="assess complexity of library vs sequencing depth by sub-sampling mapped BAM files"
    )

    parser.add_argument(
        "bam",
        help="BAM file to analyze",
    )
    parser.add_argument(
        "--sample", help="sample name to use in the statistics output", default="sample"
    )
    parser.add_argument(
        "--skim",
        help="skim through the BAM by investigating only every <skim>-th record (default=0 off)",
        default=0,
        type=int,
    )
    parser.add_argument(
        "--out-tsv",
        help="save results here. default=complexity.tsv",
        default="complexity.tsv",
    )
    parser.add_argument(
        "--filter-reads",
        default="",
        help="name of a function to be used to further stratify reads (default=off)",
    )
    # parser.add_argument(
    #     "--cutadapt-stats", default="", help="path to cutadapt.csv. Used to retrieve number of raw reads"
    # )

    args = parser.parse_args()
    return args


def main(args):
    logger = logging.getLogger("spacemake.scripts.complexity")

    reads = defaultdict(list)

    logger.info(f"processing '{args.bam}'")
    sam = pysam.AlignmentFile(args.bam, "rb", check_sq=False, threads=2)
    uniq = set()
    for read in util.timed_loop(sam.fetch(until_eof=True), logger, skim=args.skim):
        if read.query_name in uniq:
            continue

        uniq.add(read.query_name)
        cb = read.get_tag("CB")
        umi = read.get_tag("MI")

        key = cb + umi
        reads["all"].append(key)

    logger.info("finished loading all reads")

    data = []
    for tagval, keys in reads.items():
        keys = np.array(keys, dtype=object)
        N_reads = len(keys)

        n_min = 10000
        if N_reads < 5 * n_min:
            continue

        n_subsamples = np.linspace(n_min, N_reads, 20)
        for n in n_subsamples:
            n, n_UMI = subsample(keys, n)
            logger.debug(
                f"sub-sampling {args.sample}\t{tagval}\t{100.0 * n/N_reads:.2f} %\t{n}\t{n_UMI}\t{n/n_UMI:.2f}"
            )
            data.append((tagval, n, n_UMI, n / n_UMI))

    df = pd.DataFrame(data, columns=["tagval", "n_reads", "n_umis", "pcr_factor"])
    df["sample"] = args.sample
    df["bamname"] = os.path.basename(args.bam).replace(".bam", "")
    df.to_csv(util.ensure_path(args.out_tsv), sep="\t")
    return df


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    args = parse_args()
    main(args)
