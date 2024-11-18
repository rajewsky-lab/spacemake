import logging
import os

import multiprocessing as mp
import pandas as pd
import time

from spacemake.util import message_aggregation, FASTQ_src

"""
Sequence intersection between a query file and target files.

Usage:
    python n_intersect_sequences.py \
        --query input_reads.txt \
        --target target_file1.txt target_file2.txt \
        --target-id id1 id2 \
        --summary-output output_summary.txt \
        --output target_output1.txt target_output2.txt

Author:
    Daniel León-Periñán
"""

logger_name = "spacemake.snakemake.scripts.n_intersect_sequences"
logger = logging.getLogger(logger_name)
AVAILABLE_CPUS = mp.cpu_count()


def setup_parser(parser):
    """
    Set up command-line arguments for the script.

    :param parser: Argument parser object.
    :type parser: argparse.ArgumentParser
    :returns: Updated argument parser object.
    :rtype: argparse.ArgumentParser
    """
    parser.add_argument(
        "--query",
        type=str,
        help="query reads file for lookup into 'target'(s) database",
        default=None,
    )

    parser.add_argument(
        "--query-tag",
        type=str,
        help="if input is a bam file, the specified tag will be used as the sequence",
        default="CB",
    )

    parser.add_argument(
        "--query-plain-skip",
        type=int,
        help="number of lines to skip, if target input is a plain text file",
        default=None,
    )

    parser.add_argument(
        "--query-plain-column",
        type=int,
        help="column to parse from query, if input is a plain text file",
        default=None,
    )

    parser.add_argument(
        "--query-separator",
        type=str,
        default="\t",
        help="if input is a tabular plain text file, this specifies the character for column separation",
        choices=["\t", ",", ";"],
    )

    parser.add_argument(
        "--target",
        type=str,
        nargs="+",
        help="target puck files against which search is performed",
        required=True,
    )

    parser.add_argument(
        "--output",
        nargs="+",
        type=str,
        help="where to save puck files, in the same order as the --target",
        default="",
    )

    parser.add_argument(
        "--target-id",
        type=str,
        nargs="+",
        help="target puck ids for summary output",
        required=True,
    )

    parser.add_argument(
        "--target-column",
        type=str,
        help="column name to parse from target",
        default=None,
    )

    parser.add_argument(
        "--target-separator",
        type=str,
        default="\t",
        help="if input is a tabular plain text file, this specifies the character for column separation",
        choices=["\t", ",", ";"],
    )

    parser.add_argument(
        "--summary-output",
        type=str,
        help="a summary output file containing the number of matches between query and target, per target file",
        required=False,
        default="",
    )

    parser.add_argument(
        "--min-threshold",
        type=float,
        default=0,
        help="targets with 'matching_ratio' above --min-threshold are marked as True; else, False.",
        required=False,
    )

    parser.add_argument(
        "--chunksize",
        type=int,
        default=int(100 * 1e6),
        help="number of query sequences per processing chunk (default: 100M)",
        required=False,
    )

    parser.add_argument(
        "--n-jobs",
        type=int,
        help="number of parallel jobs (up to {AVAILABLE_CPUS} in this machine) for database lookup",
        metavar=f"[0-{AVAILABLE_CPUS-1}]",
        choices=range(0, AVAILABLE_CPUS),
        default=1,
    )

    return parser


# TODO: use a partial function to pass args to not use global context
def find_matches(_target, df=None):
    """
    Find matches between query and target sequences.

    :param i: iterator index for the file argument, to compare against query sequences.
    :type i: int
    :returns: Tuple containing information about the target file, number of barcodes, number of matching barcodes, and the matching ratio.
    :rtype: tuple
    """
    global query_seqs, args

    if isinstance(_target, int):
        f_in = args.target[_target]
        df = pd.read_csv(f_in, sep=args.target_separator)
        reads = df[args.target_column]
        target = set(reads)
        if args.output != "":
            f_out = args.output[_target]
    elif isinstance(_target, set):
        f_in = None
        target = _target
        f_out = args.output[0]
    else:
        raise ValueError(
            "_target must be either integer (will use as index for 'args.target' or set of unique reads"
        )

    start = time.time()
    _intersection = target.intersection(query_seqs)
    n_matches = len(_intersection)
    print(f"queried {len(target):,} barcodes and found {n_matches:,} matches in {round(time.time() - start, 2):,} seconds")

    if args.output != "":
        df_matched = df[df[args.target_column].isin(list(_intersection))]
        df_matched = df_matched[["cell_bc", "x_pos", "y_pos"]]
        df_matched.to_csv(
            f_out, mode="a", header=not os.path.exists(f_out), index=False
        )
        print(f"saved puck file into {f_out}")

    return (f_in, len(target), n_matches)


def generate_puck_barcode_summary(df):
    df_summary = pd.DataFrame(
        {
            "x_pos_min_px": [df.x_pos.min()],
            "x_pos_max_px": [df.x_pos.max()],
            "y_pos_min_px": [df.y_pos.min()],
            "y_pos_max_px": [df.y_pos.max()],
        }
    )

    return df_summary


@message_aggregation(logger_name)
def cmdline():
    """cmdline."""
    import argparse

    global args, query_seqs

    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="calculate sequence intersection between a query and target(s) file(s)",
    )
    parser = setup_parser(parser)

    args = parser.parse_args()

    if len(args.target_id) != len(args.target):
        raise ValueError(
            f"target_id ({len(args.target_id)}) and target ({len(args.target)}) are different in size"
        )

    if args.output == "" and args.summary_output == "":
        raise ValueError("One of --output or --summary-output must have a value")

    # Create a df for the final results
    result_df = pd.DataFrame()

    # Load the query file as chunks
    query_df = pd.read_csv(
        args.query,
        sep=args.query_separator,
        skiprows=args.query_plain_skip,
        chunksize=args.chunksize,
    )

    if len(args.output) == 1:
        df = pd.read_csv(args.target[0], sep=args.target_separator).rename(
            columns={"xcoord": "x_pos", "ycoord": "y_pos"}
        )
        reads = df[args.target_column]
        target = set(reads)

    # Process in chunks
    for query_df_chunk in query_df:
        query_seqs = query_df_chunk.iloc[:, args.query_plain_column]
        logger.info(
            f"read chunk of {round(len(query_seqs)/1e6, 2)}M sequences from query file at '{args.query}'"
        )

        # get unique reads into a global context
        query_seqs = set(query_seqs)
        logger.info(f"hashed unique query reads")

        if len(args.target) == 1:
            _, n_barcodes, n_matches = find_matches(target, df)
            results = [(args.target[0], n_barcodes, n_matches)]
        elif len(args.target) > 1:
            logger.info(
                f"querying against {len(args.target)} targets with {args.n_jobs} parallel jobs"
            )
            with mp.Pool(args.n_jobs) as pool:
                results = pool.map(find_matches, range(len(args.target)))

        result_df_chunk = pd.DataFrame(
            results,
            columns=["puck_barcode_file", "n_barcodes", "n_matching"],
        )
        result_df_chunk["puck_barcode_file_id"] = args.target_id
        result_df = pd.concat([result_df, result_df_chunk])

    if len(args.output) == 1 and args.summary_output != "":
        df_summary = generate_puck_barcode_summary(df)
        result_df = (
            result_df.groupby(
                ["puck_barcode_file", "puck_barcode_file_id", "n_barcodes"]
            )
            .sum()
            .reset_index()
        )
        df_summary["puck_barcode_file"] = result_df["puck_barcode_file"]
        df_summary["puck_barcode_file_id"] = result_df["puck_barcode_file_id"]
        df_summary["n_barcodes"] = result_df["n_barcodes"]
        df_summary["n_matching"] = result_df["n_matching"]
        df_summary["matching_ratio"] = (
            df_summary["n_matching"] / df_summary["n_barcodes"]
        )
        df_summary["parsed_barcode_file"] = args.output
        df_summary.to_csv(args.summary_output, index=False)

    elif args.summary_output != "":
        # group per puck_barcode_file, and compute sum
        result_df = (
            result_df.groupby(
                ["puck_barcode_file", "puck_barcode_file_id", "n_barcodes"]
            )
            .sum()
            .reset_index()
        )
        result_df["matching_ratio"] = result_df["n_matching"] / result_df["n_barcodes"]
        result_df["pass_threshold"] = 0
        result_df["pass_threshold"][
            result_df["matching_ratio"] > args.min_threshold
        ] = 1
        result_df.to_csv(args.summary_output, index=False)


if __name__ == "__main__":
    cmdline()
