import pandas as pd
import gzip
import itertools
import os
from collections import Counter

bases = ["A", "C", "T", "G", "N"]
kmer_len = int(snakemake.params["kmer_len"])

# in this case we look for 4 mers
kmers = ["".join(kmer) for kmer in itertools.product(bases, repeat=kmer_len)]

position_counts = None
read_len = 0
position_list = []

read_kmer_hashes = []

with gzip.open(snakemake.input[0], "rt") as fastq_in:
    line = 0
    for read in fastq_in:
        if line == 1:
            read = read.strip("\n")
            read_len = len(read)
            position_list = list(range(read_len - kmer_len + 1))
            kmer_hashes = [
                "_".join(prod)
                for prod in itertools.product(kmers, [str(x) for x in position_list])
            ]
            position_counts = pd.DataFrame(0, index=kmer_hashes, columns=["count"])

        # if line is a read
        if line % 4 == 1:
            kmer_hashes = kmer_hashes + [
                str(read[i : i + kmer_len]) + "_" + str(i) for i in position_list
            ]

        line = line + 1
        if line % 4000 == 0:
            kmer_hash_counts = Counter(kmer_hashes)
            # print(kmer_hash_counts.values())

            # update df
            position_counts.loc[kmer_hash_counts.keys(), "count"] = position_counts.loc[
                kmer_hash_counts.keys(), "count"
            ] + list(kmer_hash_counts.values())
            # print(position_counts)

            kmer_hashes = []

        if line % 4000000 == 0:
            print("%s reads processed" % (line / 4))

position_counts.index.rename("kmer_hash", inplace=True)

file_path = snakemake.output[0]
os.makedirs(os.path.dirname(file_path), exist_ok=True)

position_counts.to_csv(file_path)
