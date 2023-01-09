import pysam
import numpy as np
import os
import argparse
import logging
from collections import defaultdict
import pandas as pd
import spacemake.util as util
from spacemake.contrib import __license__, __version__, __author__, __email__


def parse_args():
    parser = util.make_minimal_parser(
        prog="complexity.py",
        description="assess complexity of library vs sequencing depth by sub-sampling mapped BAM files"
    )

    parser.add_argument(
        "bam",
        help="BAM file to analyze",
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
        "--chunk-size",
        type=int,
        default=10000000,
        help="number of reads to load at one time (default=10M). This saves a lot of RAM."
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

def subsample(reads, n):
    f = float(n) / len(reads)
    n_uniq = 0
    n_reads = 0
    last = None

    # rnd = np.random.random(size=len(reads))
    # mask = rnd <= f
    import random
    for r in reads:
        rnd = random.random()
        if rnd > f:
            # skip the read
            continue

        n_reads += 1
        if r != last:
            n_uniq += 1
        last = r

    return n_reads, n_uniq


def main(args):
    logger = util.setup_logging(args, "spacemake.scripts.complexity")

    logger.info(f"processing '{args.bam}'")
    sam = util.quiet_bam_open(args.bam, "rb", check_sq=False, threads=2)

    last_qname = None
    keys = []

    chunk_size = args.chunk_size

    chunks = []
    chunk = np.empty(chunk_size, dtype=np.uint64)
    chunk_free = chunk_size

    for read in util.timed_loop(sam.fetch(until_eof=True), logger, skim=args.skim):
        if read.query_name == last_qname:
            continue

        last_qname = read.query_name
        if chunk_free <= 0:
            chunks.append(chunk)
            chunk = np.empty(chunk_size, dtype=np.uint64)
            chunk_free = chunk_size

        # tags = dict(read.get_tags())
        # cb = tags['CB']
        # umi = tags['MI']
        # # Above code has less than half the reads/second! Use get_tag() instead
        # gn = tags.get('gn', 'n/a') # this might lead to overestimates of complexity. Let's ignore gene, hoping that UMI space is large enough
        cb = read.get_tag("CB")
        umi = read.get_tag("MI")
        key = cb + umi

        chunk[chunk_size - chunk_free] = hash(key)
        chunk_free -= 1

        # keys.append(hash(key)) # substitute a 64bit hash for the real strings. 
        # + saves huge amounts of RAM 
        # + makes sorting faster. 
        # - chance of hash collisions. But should be absolutely negligible as 
        #   long as 2^64 is >> number of sequenced molecules 
        #   (which should be true for the foreseeable future)

    chunks.append(chunk[:chunk_size - chunk_free])
    print(f"concatenating {len(chunks)} chunks")

    logger.debug(f"finished loading all reads. Concatenating {len(chunks)} chunks.")
    keys = np.concatenate(chunks)
    chunks = None # free some memory

    logger.debug("sorting...")
    keys.sort()
    N_reads = len(keys)

    tagval = "all"

    logger.debug("sub-sampling and unique sequence counting")

    n_min = 10000
    data = []
    if N_reads >= 5 * n_min:
        n_subsamples = np.linspace(n_min, N_reads, 20)
        for n in n_subsamples:
            n, n_UMI = subsample(keys, n)
            logger.info(
                f"sub-sampling {args.sample}\t{tagval}\t{100.0 * n/N_reads:.2f} %\t{n}\t{n_UMI}\t{n/n_UMI:.2f}"
            )
            data.append((tagval, n, n_UMI, n / n_UMI))
    else:
        n, n_UMI = subsample(keys, N_reads)
        data.append((tagval, n, n_UMI, n / n_UMI))

    df = pd.DataFrame(data, columns=["tagval", "n_reads", "n_umis", "pcr_factor"])
    df["sample"] = args.sample
    df["bamname"] = os.path.basename(args.bam).replace(".bam", "")
    df.to_csv(util.ensure_path(args.out_tsv), sep="\t")
    return df


if __name__ == "__main__":
    args = parse_args()
    main(args)
