import numpy as np

# import cutadapt.align
from collections import defaultdict

from spacemake.contrib import __version__, __license__, __author__, __email__
import multiprocessing as mp
from spacemake.parallel import (
    put_or_abort,
    queue_iter,
    join_with_empty_queues,
    chunkify,
    ExceptionLogging,
    log_qerr,
)
import spacemake.util as util
from time import time
import pysam
import logging
import os


def parse_cmdline():
    import argparse

    parser = argparse.ArgumentParser(
        description="trim adapters from a BAM file using cutadapt"
    )
    parser.add_argument(
        "bam_in",
        help="bam input (default=stdin)",
        default="/dev/stdin",
        # nargs="+",
    )
    parser.add_argument(
        "--bam-out",
        help="bam output (default=stdout)",
        default="/dev/stdout",
    )
    parser.add_argument(
        "--bam-out-mode",
        help="bam output mode (default=b0)",
        default="b0",
    )
    parser.add_argument(
        "--adapters-right",
        help="FASTA file with adapter sequences to trim from the right (3') end of the reads",
        default="",
    )
    parser.add_argument(
        "--adapters-left",
        help="FASTA file with adapter sequences to trim from the left (5') end of the reads",
        default="",
    )
    parser.add_argument(
        "--skim",
        help="skim through the BAM by investigating only every <skim>-th record (default=1 off)",
        default=1,
        type=int,
    )
    parser.add_argument(
        "--min-length",
        help="minimal allowed read-length left after trimming (default=18)",
        type=int,
        default=18,
    )
    parser.add_argument(
        "--min-qual",
        help="minimal quality score for quality trimming (default=20)",
        type=int,
        default=20,
    )
    parser.add_argument(
        "--phred-base",
        help="phred score base used for qual trimming (default=33)",
        type=int,
        default=33,
    )

    parser.add_argument(
        "--threads-read",
        help="number of threads for reading bam_in (default=2)",
        type=int,
        default=2,
    )
    parser.add_argument(
        "--threads-write",
        help="number of threads for writing bam_out (default=2)",
        type=int,
        default=2,
    )
    parser.add_argument(
        "--threads-work",
        help="number of worker threads for actual trimming (default=1)",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--chunk-size",
        help="number of bam records per chunk. Chunks are processes in parallel by workers (default=100000)",
        type=int,
        default=100000,
    )

    parser.add_argument(
        "--stats-out",
        help="write tab-separated table with trimming results here",
        default="",
    )
    return parser.parse_args()


def load_adapters(right, left, max_errors=0.1, min_overlap=3):
    import cutadapt.adapters
    from spacemake.util import fasta_chunks
    import re

    def parse_seq_id(seq_id):
        err = max_errors
        ov = min_overlap
        name = seq_id.split()[0]
        m = re.search(r"max_errors=(\S+)", seq_id)
        if m:
            # overload max_errors default value
            err = float(m.groups()[0])

        m = re.search(r"min_overlap=(\S+)", seq_id)
        if m:
            # overload min_overlap default value
            ov = int(m.groups()[0])

        return name, err, ov

    adapters_right = []
    if args.adapters_right:
        for seq_id, seq in fasta_chunks(open(right)):
            name, err, ov = parse_seq_id(seq_id)
            adapters_right.append(
                (
                    name,
                    seq,
                    cutadapt.adapters.BackAdapter(
                        seq, name=name, max_errors=err, min_overlap=ov
                    ),
                )
            )

    adapters_left = []
    if args.adapters_left:
        for seq_id, seq in fasta_chunks(open(left)):
            name, err, ov = parse_seq_id(seq_id)
            adapters_left.append(
                (
                    name,
                    seq,
                    cutadapt.adapters.NonInternalFrontAdapter(
                        seq, name=name, max_errors=err, min_overlap=ov
                    ),
                )
            )

    return adapters_right, adapters_left


class SimpleRead:
    def __init__(self, name, seq, qual, tags):
        self.query_name = name
        self.query_sequence = seq
        self.query_qualities = qual
        self.tags = tags

    @classmethod
    def from_BAM(cls, read):
        return cls(
            read.query_name,
            read.query_sequence,
            read.query_qualities,
            dict(read.get_tags()),
        )

    def set_tag(self, tag, value):
        self.tags[tag] = value

    @staticmethod
    def iter_BAM(bam_src):
        for read in bam_src:
            yield SimpleRead.from_BAM(read)

    @staticmethod
    def iter_to_BAM(sr_src, header=None):
        for read in sr_src:
            aln = pysam.AlignedSegment(header)
            aln.query_name = read.query_name
            aln.query_sequence = read.query_sequence
            aln.query_qualities = read.query_qualities
            aln.flag = 4
            aln.tags = list(read.tags.items())

            yield aln


def process_reads(read_source, args, stats={}, total={}, lhist={}):
    adapters_right, adapters_left = load_adapters(
        args.adapters_right, args.adapters_left
    )

    def check_discard(start, end, reason):
        if (end - start) < args.min_length:
            stats[reason] += 1
            stats["N_discarded"] += 1
            total["bp_discarded"] += end - start

            return True
        else:
            return False

    for read in read_source:
        # read_seq = read.query_sequence
        # read_qual = read.query_qualities
        # we got a string
        cols = read.split("\t")
        read_seq = cols[9]
        qual_str = cols[10]
        read_qual = np.array(bytearray(qual_str.encode("ASCII"))) - args.phred_base

        # print(read_qual)
        start = 0
        end = len(read_seq)
        tags = []

        stats["N_input"] += 1
        total["bp_input"] += end
        trimmed_names_right = []
        trimmed_bases_right = []

        trimmed_names_left = []
        trimmed_bases_left = []

        # quality trimming
        qtrim = np.array(read_qual) >= args.min_qual
        # print(qtrim)
        # print(
        #     f"qtrim.argmax()={qtrim.argmax()} L={end} new_end=qtrim[::-1].argmax() ={qtrim[::-1].argmax()}"
        # )
        # the first position where q >= min_q (from the 3' end) should be the new end.
        n_trimmed = (qtrim[::-1]).argmax()
        if n_trimmed:
            end -= n_trimmed
            trimmed_bases_right.append(n_trimmed)
            trimmed_names_right.append("Q")

            stats["N_Qtrimmed"] += 1
            total["bp_Qtrimmed"] += n_trimmed
            total["bp_trimmed"] += n_trimmed

        if check_discard(start, end, "N_too_short_after_Q"):
            continue

        if adapters_left:
            # left end adapter trimming
            for adap_name, adap_seq, adap in adapters_left:
                match = adap.match_to(read_seq[start:end])
                if match:
                    # print(adap_name, adap, match)
                    new_start = max(start, match.rstop)
                    n_trimmed = new_start - start
                    start = new_start
                    trimmed_bases_left.append(n_trimmed)
                    trimmed_names_left.append(adap_name)

                    stats["N_" + adap_name] += 1
                    total["bp_" + adap_name] += n_trimmed
                    total["bp_trimmed"] += n_trimmed

            if check_discard(start, end, f"N_too_short_after_left_{adap_name}"):
                continue

        if adapters_right:
            # right end adapter trimming
            for adap_name, adap_seq, adap in adapters_right:
                match = adap.match_to(read_seq[start:end])
                if match:
                    new_end = min(end, match.rstart)
                    n_trimmed = end - new_end
                    end = new_end
                    trimmed_bases_right.append(n_trimmed)
                    trimmed_names_right.append(adap_name)

                    stats["N_" + adap_name] += 1
                    total["bp_" + adap_name] += n_trimmed
                    total["bp_trimmed"] += n_trimmed

            if check_discard(start, end, f"N_too_short_after_right_{adap_name}"):
                continue

        # we've made it to the end!
        stats["N_kept"] += 1
        total["bp_kept"] += end
        lhist[end] += 1
        # print(f"keeping read up to {end}")
        # read.query_sequence = read_seq[start:end]
        # read.query_qualities = read_qual[start:end]
        cols[9] = read_seq[start:end]
        cols[10] = qual_str[start:end]

        if trimmed_names_right:
            tags.append(f"A3:Z:{','.join(trimmed_names_right)}")
            tags.append(f"T3:Z:{','.join([str(s) for s in trimmed_bases_right])}")

        if trimmed_names_left:
            tags.append(f"A5:Z:{','.join(trimmed_names_left)}")
            tags.append(f"T5:Z:{','.join([str(s) for s in trimmed_bases_left])}")

        if tags:
            cols[-1] = cols[-1] + " " + " ".join(tags)

        yield "\t".join(cols)
        # bam_out.write(read)
        # print(read)


def skim_reads(read_source, skim):
    for i, read in enumerate(read_source):
        if skim and i % skim != 0:
            continue
        else:
            yield read


def BAM_to_string(src):
    for rec in src:
        yield rec.tostring()


def string_to_BAM(src, header=None):
    for s in src:
        yield pysam.AlignedSegment.fromstring(s, header)


def main_single(args):
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger("spacemake.cutadapt_bam.main_single")

    bam_in = pysam.AlignmentFile(
        args.bam_in, "rb", check_sq=False, threads=args.threads_read
    )
    bam_out = pysam.AlignmentFile(
        args.bam_out,
        f"w{args.bam_out_mode}",
        header=util.make_header(bam_in, progname=os.path.basename(__file__)),
        threads=args.threads_write,
    )

    stats = defaultdict(int)
    total = defaultdict(int)
    lhist = defaultdict(int)

    t0 = time()
    for read in process_reads(
        skim_reads(bam_in.fetch(until_eof=True), args.skim),
        args,
        stats=stats,
        total=total,
        lhist=lhist,
    ):
        bam_out.write(read)

    dt = time() - t0
    logger.info(
        f"processed {stats['N_input']} reads in {dt:.1f} seconds ({stats['N_input']/dt:.1f} reads/second)."
    )

    if args.stats_out:
        with open(util.ensure_path(args.stats_out), "wt") as f:
            f.write("key\tcount\tpercent\n")
            for k, v in sorted(stats.items(), key=lambda x: -x[1]):
                f.write(f"reads\t{k}\t{v}\t{100.0 * v/stats['N_input']:.2f}\n")

            for k, v in sorted(total.items(), key=lambda x: -x[1]):
                f.write(f"bases\t{k}\t{v}\t{100.0 * v/total['bp_input']:.2f}\n")

            for k, v in sorted(lhist.items()):
                f.write(f"L_final\t{k}\t{v}\t{100.0 * v/stats['N_kept']:.2f}\n")


## Parallel implementation
def parallel_read(Qsam, args, Qerr, abort_flag, stat_list):
    """
    reads from BAM file, converts to string, groups the records into chunks for
    faster parallel processing, and puts these on a mp.Queue()
    """
    with ExceptionLogging(
        "spacemake.cutadapt_bam.parallel_read", Qerr=Qerr, exc_flag=abort_flag
    ) as el:
        bam_in = pysam.AlignmentFile(
            args.bam_in, "rb", check_sq=False, threads=args.threads_read
        )
        bam_header = stat_list[-1]
        bam_header.append(util.make_header(bam_in, progname=os.path.basename(__file__)))

        # read_source = BAM_to_string(skim_reads(bam_in.fetch(until_eof=True), args.skim))
        if args.skim:
            read_source = BAM_to_string(
                skim_reads(bam_in.fetch(until_eof=True), args.skim)
            )
        else:
            read_source = BAM_to_string(bam_in.fetch(until_eof=True))

        for chunk in chunkify(read_source, n_chunk=args.chunk_size):
            el.logger.debug(
                f"placing {chunk[0]} {len(chunk[1])} in queue of depth {Qsam.qsize()}"
            )
            if put_or_abort(Qsam, chunk, abort_flag):
                el.logger.warning("shutdown flag was raised!")
                break


def parallel_trim(Qsam, Qres, args, Qerr, abort_flag, stat_list):
    with ExceptionLogging(
        "spacemake.cutadapt_bam.parallel_trim", Qerr=Qerr, exc_flag=abort_flag
    ) as el:
        el.logger.debug(f"parallel_trim() starting up with args={args}")

        _stats = defaultdict(int)
        _total = defaultdict(int)
        _lhist = defaultdict(int)

        for n_chunk, reads in queue_iter(Qsam, abort_flag):
            result = list(
                process_reads(
                    reads,
                    args,
                    stats=_stats,
                    total=_total,
                    lhist=_lhist,
                )
            )
            Qres.put((n_chunk, result))

        el.logger.debug("synchronizing cache and counts")
        stats, total, lhist = stat_list[:3]
        stats.append(_stats)
        total.append(_total)
        lhist.append(_lhist)


def process_ordered_results(res_queue, args, Qerr, abort_flag, stat_list, timeout=10):
    with ExceptionLogging(
        "spacemake.cutadapt_bam.process_ordered_results", Qerr=Qerr, exc_flag=abort_flag
    ) as el:
        import heapq
        import time

        heap = []
        n_chunk_needed = 0
        t0 = time.time()
        t1 = t0
        n_rec = 0

        logger = el.logger
        # out = Output(args)
        bam_header = stat_list[-1]
        while (not len(bam_header)) and timeout > 0:
            if abort_flag:
                logger.info("aborting during startup...")
                return

            logger.info("waiting for header to get ready")
            time.sleep(1)
            timeout -= 1

        if timeout < 0:
            raise ValueError("header was not made available. Did the dispatcher die?")

        header = bam_header.pop()
        bam_out = pysam.AlignmentFile(
            args.bam_out,
            f"w{args.bam_out_mode}",
            header=header,
            threads=args.threads_write,
        )

        for n_chunk, results in queue_iter(res_queue, abort_flag):
            heapq.heappush(heap, (n_chunk, results))

            # as long as the root of the heap is the next needed chunk
            # pass results on to storage
            while heap and (heap[0][0] == n_chunk_needed):
                n_chunk, results = heapq.heappop(heap)  # retrieves heap[0]
                # for aln in SimpleRead.iter_to_BAM(results, header=bam_out.header):
                for aln in string_to_BAM(results, header=bam_out.header):
                    #     # print("record in process_ordered_results", record)
                    bam_out.write(aln)
                    n_rec += 1

                n_chunk_needed += 1

            # debug output on average throughput
            t2 = time.time()
            if t2 - t1 > 30:
                dT = t2 - t0
                rate = n_rec / dT
                logger.info(
                    "processed {0} reads in {1:.0f} seconds (average {2:.0f} reads/second).".format(
                        n_rec, dT, rate
                    )
                )
                t1 = t2

        # out.close()
        # by the time None pops from the queue, all chunks
        # should have been processed!
        if not abort_flag.value:
            assert len(heap) == 0
        else:
            logger.warning(
                f"{len(heap)} chunks remained on the heap due to missing data upon abort."
            )

        dT = time.time() - t0
        logger.info(
            "finished processing {0} reads in {1:.0f} seconds (average {2:.0f} reads/second)".format(
                n_rec, dT, n_rec / dT
            )
        )


## TODO: push these two into parallel/util modules


def count_dict_sum(sources):
    dst = defaultdict(float)
    for src in sources:
        for k, v in src.items():
            dst[k] += v

    return dst


def main_parallel(args):
    logging.basicConfig(level=logging.DEBUG)

    # queues for communication between processes
    Qsam = mp.Queue(args.threads_work * 10)  # reads from parallel_read->parallel_trim
    Qres = mp.Queue()  # extracted BCs from process_combinatorial->collector
    Qerr = mp.Queue()  # child-processes can report errors back to the main process here

    # Proxy objects to allow workers to report statistics about the run
    manager = mp.Manager()
    abort_flag = mp.Value("b")
    abort_flag.value = False

    stats = manager.list()
    total = manager.list()
    lhist = manager.list()
    bam_header = manager.list()
    stat_lists = [stats, total, lhist, bam_header]
    print(stat_lists[-1])
    with ExceptionLogging(
        "spacemake.cutadapt_bam.main_parallel", exc_flag=abort_flag
    ) as el:

        # read BAM in chunks and put them in Qsam
        dispatcher = mp.Process(
            target=parallel_read,
            name="dispatcher",
            args=(Qsam, args, Qerr, abort_flag, stat_lists),
        )

        dispatcher.start()
        el.logger.info("Started dispatch")

        # workers consume chunks of BAM from Qsam
        # process them, and put the results in Qres
        workers = []
        for i in range(args.threads_work):
            w = mp.Process(
                target=parallel_trim,
                name=f"worker_{i}",
                args=(Qsam, Qres, args, Qerr, abort_flag, stat_lists),
            )
            w.start()
            workers.append(w)

        el.logger.info("Started workers")

        collector = mp.Process(
            target=process_ordered_results,
            name="output",
            args=(Qres, args, Qerr, abort_flag, stat_lists),
        )
        collector.start()
        el.logger.info("Started collector")
        # wait until all sequences have been thrown onto Qfq
        qfq, qerr = join_with_empty_queues(dispatcher, [Qsam, Qerr], abort_flag)
        el.logger.info("The dispatcher exited")
        if qfq or qerr:
            el.logger.info(f"{len(qfq)} chunks were drained from Qfq upon abort.")
            log_qerr(qerr)

        # signal all workers to finish
        el.logger.info("Signalling all workers to finish")
        for n in range(args.threads_work):
            Qsam.put(None)  # each worker consumes exactly one None

        for w in workers:
            # make sure all results are on Qres by waiting for
            # workers to exit. Or, empty queues if aborting.
            qres, qerr = join_with_empty_queues(w, [Qres, Qerr], abort_flag)
            if qres or qerr:
                el.logger.info(f"{len(qres)} chunks were drained from Qres upon abort.")
                log_qerr(qerr)

        el.logger.info(
            "All worker processes have joined. Signalling collector to finish."
        )
        # signal the collector to stop
        Qres.put(None)

        # and wait until all output has been generated
        collector.join()
        el.logger.info("Collector has joined. Merging worker statistics.")

    if args.stats_out:
        stats = count_dict_sum(stats)
        total = count_dict_sum(total)
        lhist = count_dict_sum(lhist)

        with open(util.ensure_path(args.stats_out), "wt") as f:
            f.write("key\tcount\tpercent\n")
            for k, v in sorted(stats.items(), key=lambda x: -x[1]):
                f.write(f"reads\t{k}\t{v}\t{100.0 * v/stats['N_input']:.2f}\n")

            for k, v in sorted(total.items(), key=lambda x: -x[1]):
                f.write(f"bases\t{k}\t{v}\t{100.0 * v/total['bp_input']:.2f}\n")

            for k, v in sorted(lhist.items()):
                f.write(f"L_final\t{k}\t{v}\t{100.0 * v/stats['N_kept']:.2f}\n")

    if el.exception:
        return -1


if __name__ == "__main__":
    args = parse_cmdline()
    if args.threads_work == 1:
        ret_code = main_single(args)
    else:
        ret_code = main_parallel(args)

    import sys

    sys.exit(ret_code)
