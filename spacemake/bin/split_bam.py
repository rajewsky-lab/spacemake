import os
import spacemake.util as util
import pyximport

pyximport.install()
import spacemake.cython.fast_loop as fast_loop
from spacemake.parallel import ExceptionLogging, create_named_pipes

# pyximport.install(setup_args={"include_dirs": np.get_include()}, reload_support=True)
# import spacemake.cython.fastread as fr


def bam_decompressor(bam_in, sam_out, threads=2):
    import os

    with ExceptionLogging("spacemake.parallel.split_bam") as el:
        cmd = f"samtools view -h --threads={threads} {bam_in} > {sam_out}"
        el.logger.info(f"executing '{cmd}'")
        os.system(cmd)
        el.logger.info(f"finished decompressing.")


def sam_distributor(in_pipe, worker_in_pipes, chunk_size=1000):
    with ExceptionLogging("spacemake.parallel.sam_distributor") as el:
        el.logger.info(
            f"reading from {in_pipe}, writing to {worker_in_pipes} chunk_size={chunk_size}"
        )

        n = fast_loop.distribute(
            in_pipe,
            worker_in_pipes,
            chunk_size=chunk_size,
            header_detect_func=lambda line: line.startswith("@"),
            header_broadcast=True,
        )
        el.logger.info(f"finished distributing {n} BAM records.")


def sam_CB_distributor(in_pipe, worker_in_pipes, cb_dict={b'A':0, b'C':1, b'G':2, b'T':3, b'N': 0}):
    with ExceptionLogging("spacemake.parallel.sam_CB_distributor") as el:
        el.logger.info(
            f"reading from {in_pipe}, writing to {worker_in_pipes} cb_dict={cb_dict}"
        )

        sub_size = len(list(cb_dict.keys())[0])
        n = fast_loop.distribute_by_substr(
            in_pipe,
            worker_in_pipes,
            cb_dict,
            sub_size,
            header_detect_func=lambda line: line.startswith("@"),
            header_broadcast=True,
        )
        el.logger.info(f"finished distributing {n} BAM records.")


def worker_read_sam_write_bam(in_pipe, out_pipe, mode="hSb", threads=4):
    import os

    with ExceptionLogging("spacemake.parallel.worker_read_sam_write_bam") as el:
        cmd = f"samtools view -{mode} --threads={threads} {in_pipe} > {out_pipe}"
        el.logger.info(f"executing '{cmd}'")
        os.system(cmd)

    el.logger.info(f"finished writing to '{out_pipe}'.")


## Helper functions
def make_cb_prefix_dict(L=4, n=10, include_N=True):
    lkup = {}
    for i, prefix in enumerate(util.generate_kmers(L, "ACGTN")):
        lkup[prefix.encode('ascii')] = i % n
    
    return lkup

def load_cb_routing(fname, out_folder, name, ext):
    import pandas as pd
    groups_seen = {}
    n_group = 0

    out_paths = []

    def numeric(key):
        nonlocal n_group
        if key in groups_seen:
            return groups_seen[key]
        else:
            groups_seen[key] = n_group
            n_group += 1
            out_paths.append(
                os.path.join(out_folder, f"{name}_{key}{ext}")
            )
            return groups_seen[key]

    df = pd.read_csv(fname, sep='\t').set_index('CB').sort_values('group')
    df['numeric'] = df['group'].apply(numeric)
    
    lkup = df['numeric'].to_dict()
    lkup = dict([(k.encode('ascii'), v) for k, v in lkup.items()])
    return lkup, out_paths


def store_cb_routing(lkup, fname):
    import pandas as pd
    # convert bytes to str
    d = dict([(k.decode('ascii'), v) for k, v in lkup.items()])

    # create two columns
    df = pd.DataFrame(dict(group=pd.Series(d))).reset_index(names="CB")
    df.to_csv(fname, sep='\t', index=False)

    return df


def parse_args():
    import spacemake.util as util

    parser = util.make_minimal_parser(
        "split_bam.py",
        description="break one big BAM file into multiple smaller ones",
    )
    parser.add_argument(
        "bam_in",
    )
    parser.add_argument(
        "-n", type=int, default=10, help="number of BAM files to produce (default=10)"
    )
    parser.add_argument(
        "--mode", default="rr", choices=['rr', 'cb-split', 'cb-group'], help=(
            "how would you like your BAM file split up into smaller BAMs?\n\n"
            "\trr (default)= round robin: cycle between -n <n> output files, record i is assigned to file j = i %% n\n"
            "\tcb-split = cell-barcode prefix split: read the first --cb-prefix-len <L> characters of the cell barcode "
            "\tBAM tag 'CB:Z:<cb>' and split the records based on that value among the -n <n> files as evenly as possible."
            "\tin this mode, the resulting split files are guaranteed to contain non-overlapping sets of barcodes.\n"
            "\tcb-group = cell-barcode custom routing: read the --cb-demux-list <filename> file to route sets of "
            "\tpre-defined cell barcodes into specific output files"
        )
    )
    parser.add_argument(
        "--cb-tag", default='\tCB:Z:', help="how to find the start of the cell barcode BAM tag? default='\tCB:Z:'"
    )
    parser.add_argument(
        "--cb-prefix-len", type=int, default=4,
    )
    parser.add_argument(
        "--cb-demux-table", default="", help=(
            "provide a tab-separated file with cell-barcode (or CB-prefixes) to out-file association. "
            "Expects the columns 'CB' and 'group'"
        )
    )
    parser.add_argument(
        "--cb-demux-out", default="", help=(
            "store a tab-separated file with cell-barcode (or CB-prefixes) to out-file association. "
            "only used by --mode cb-split"
        )
    )

    parser.add_argument(
        "--decompression-threads",
        type=int,
        default=2,
        help="number of threads to use for BAM decompression (default=2)",
    )

    parser.add_argument(
        "--recompression-threads",
        type=int,
        default=1,
        help="number of threads to use for BAM recompression (default=1)",
    )

    parser.add_argument("--out-folder", "-o", default=".")

    return parser.parse_args()


def split_bam_main():
    args = parse_args()
    util.setup_logging(args)
    import multiprocessing as mp

    pipe_names = ["sam_in"] + [f"in_{n}" for n in range(args.n)]

    base = os.path.dirname(args.bam_in)
    name, ext = os.path.splitext(os.path.basename(args.bam_in))


    # BAM(file) -> [bam_decompressor] -> SAM(FIFO)
    #   -> [sam_distributor] -> n x SAM(FIFO)
    #   -> n x [worker_read_sam_write_bam] -> n x BAM(file)

    out_paths = [
        os.path.join(args.out_folder, f"{name}_{n:03d}{ext}") for n in range(args.n)
    ]
    if args.mode == "rr":
        distributor = sam_distributor
        dist_kw = {}
    elif args.mode == "cb-split":
        distributor = sam_CB_distributor
        lkup = make_cb_prefix_dict(args.cb_prefix_len, args.n)
        dist_kw = {'cb_dict' : lkup}
        if args.cb_demux_out:
            store_cb_routing(lkup, args.cb_demux_out)

    elif args.mode == "cb-group":
        distributor = sam_CB_distributor
        lkup, out_paths = load_cb_routing(args.cb_demux_table, args.out_folder, name, ext)
        dist_kw = {'cb_dict': lkup}

    with create_named_pipes(pipe_names) as pipe_paths:
        p_decomp = mp.Process(
            target=bam_decompressor,
            args=(args.bam_in, pipe_paths[0]),
            kwargs=dict(threads=args.decompression_threads),
        )
        p_dist = mp.Process(
            target=distributor, args=(pipe_paths[0], pipe_paths[1:]), kwargs=dist_kw
        )
        workers = []
        for w_in, w_out in zip(pipe_paths[1:], out_paths):
            w = mp.Process(
                target=worker_read_sam_write_bam,
                args=(w_in, w_out),
                kwargs=dict(threads=args.recompression_threads),
            )
            w.start()
            workers.append(w)

        p_dist.start()
        p_decomp.start()

        p_decomp.join()
        p_dist.join()

        for w in workers:
            w.join()


if __name__ == "__main__":
    split_bam_main()
