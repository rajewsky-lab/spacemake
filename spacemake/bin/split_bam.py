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


def worker_read_sam_write_bam(in_pipe, out_pipe, mode="hSb", threads=4):
    import os

    with ExceptionLogging("spacemake.parallel.worker_read_sam_write_bam") as el:
        cmd = f"samtools view -{mode} --threads={threads} {in_pipe} > {out_pipe}"
        el.logger.info(f"executing '{cmd}'")
        os.system(cmd)

    el.logger.info(f"finished writing to '{out_pipe}'.")


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

    out_paths = [
        os.path.join(args.out_folder, f"{name}_{n:03d}{ext}") for n in range(args.n)
    ]

    # BAM(file) -> [bam_decompressor] -> SAM(FIFO)
    #   -> [sam_distributor] -> n x SAM(FIFO)
    #   -> n x [worker_read_sam_write_bam] -> n x BAM(file)

    with create_named_pipes(pipe_names) as pipe_paths:
        p_decomp = mp.Process(
            target=bam_decompressor,
            args=(args.bam_in, pipe_paths[0]),
            kwargs=dict(threads=args.decompression_threads),
        )
        p_dist = mp.Process(
            target=sam_distributor, args=(pipe_paths[0], pipe_paths[1:])
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
