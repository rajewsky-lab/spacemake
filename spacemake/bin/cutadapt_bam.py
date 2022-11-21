from spacemake.contrib import __version__, __license__, __author__, __email__
import numpy as np

from collections import defaultdict

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
    parser = util.make_minimal_parser(
        prog="cutadapt_bam.py",
        description="trim adapters from a BAM file using cutadapt",
    )
    parser.add_argument("--config", default="config.yaml", help="path to config-file")

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
        "--adapter-flavor",
        help="name of the adapter flavor used to retrieve sequences and parameters from the config.yaml",
        default="default",
    )

    # parser.add_argument(
    #     "--adapters-right",
    #     help="FASTA file with adapter sequences to trim from the right (3') end of the reads",
    #     default="",
    # )
    # parser.add_argument(
    #     "--adapters-left",
    #     help="FASTA file with adapter sequences to trim from the left (5') end of the reads",
    #     default="",
    # )
    parser.add_argument(
        "--skim",
        help="skim through the BAM by investigating only every <skim>-th record (default=1 off)",
        default=1,
        type=int,
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
    args = parser.parse_args()
    return util.load_config_with_fallbacks(args)


class QualityTrim:
    def __init__(self, min_base_qual=20):
        self.min_base_qual = min_base_qual
        self.name = "Q"

    def match_to(self, seq, qual):
        qtrim = np.array(qual) >= self.min_base_qual
        n_trimmed = (qtrim[::-1]).argmax()
        return n_trimmed


class AdapterTrim:
    def __init__(self, name="na", where="right", seq=None, **kw):
        import cutadapt.adapters

        self.name = name
        self.where = where
        if where == "left":
            self.adapter = cutadapt.adapters.NonInternalFrontAdapter(
                seq, name=name, **kw
            )
        else:
            self.adapter = cutadapt.adapters.BackAdapter(seq, name=name, **kw)

    def match_to(self, seq, qual):
        match = self.adapter.match_to(seq)
        n_trimmed = 0
        if match:
            if self.where == "left":
                n_trimmed = match.rstop
            else:
                n_trimmed = len(seq) - match.rstart

        return n_trimmed


class AdapterFlavor:
    def __init__(
        self,
        args,
        default_max_errors=0.1,
        default_min_overlap=3,
        default_min_base_qual=20,
        default_min_read_length=18,
        default_paired_end="single_end",
        stats={},
        total={},
    ):
        # import cutadapt.adapters

        self.stats = stats
        self.total = total
        self.flavor = args.adapter_flavor
        if not args.adapter_flavor in args.config["adapter_flavors"]:
            raise KeyError(
                f"adapter_flavor '{args.adapter_flavor}' not found in config.yaml! Need valid --adapter-flavor=... "
            )

        flavor_d = args.config["adapter_flavors"][self.flavor]
        # from pprint import pprint

        # pprint(flavor_d)
        self.adapter_sequences = args.config["adapters"]

        self.trimmers_right = []
        for adap_d in flavor_d.get("cut_right"):
            for name, param_d in adap_d.items():
                if name == "Q":
                    adapter = QualityTrim(
                        min_base_qual=param_d.get(
                            "min_base_qual", default_min_base_qual
                        )
                    )
                else:
                    adapter = AdapterTrim(
                        name=name,
                        seq=self.adapter_sequences[name],
                        where="right",
                        max_errors=param_d.get("max_errors", default_max_errors),
                        min_overlap=param_d.get("min_overlap", default_min_overlap),
                    )

                self.trimmers_right.append(adapter)

        self.trimmers_left = []
        for adap_d in flavor_d.get("cut_left"):
            for name, param_d in adap_d.items():
                adapter = AdapterTrim(
                    name=name,
                    seq=self.adapter_sequences[name],
                    where="left",
                    max_errors=param_d.get("max_errors", default_max_errors),
                    min_overlap=param_d.get("min_overlap", default_min_overlap),
                )
                self.trimmers_left.append(adapter)

        self.min_read_length = flavor_d.get("min_read_length", default_min_read_length)
        self.paired_end = flavor_d.get("paired_end", default_paired_end)

    def process_read(self, read_seq, read_qual):
        # print(read_qual)
        start = 0
        end = len(read_seq)

        self.stats["N_input"] += 1
        self.total["bp_input"] += end
        trimmed_names_right = []
        trimmed_bases_right = []

        trimmed_names_left = []
        trimmed_bases_left = []

        def check_discard(start, end, reason):
            if (end - start) < self.min_read_length:
                self.stats[reason] += 1
                self.stats["N_discarded"] += 1
                self.total["bp_discarded"] += end - start

                return True
            else:
                return False

        for trimmer in self.trimmers_right:
            n_trimmed = trimmer.match_to(read_seq[start:end], read_qual[start:end])
            if n_trimmed:
                end -= n_trimmed
                trimmed_bases_right.append(n_trimmed)
                trimmed_names_right.append(trimmer.name)

                self.stats[f"N_{trimmer.name}_trimmed"] += 1
                self.total[f"bp_{trimmer.name}_trimmed"] += n_trimmed
                self.total[f"bp_trimmed"] += n_trimmed

            if check_discard(start, end, f"N_too_short_after_{trimmer.name}"):
                return

        for trimmer in self.trimmers_left:
            n_trimmed = trimmer.match_to(read_seq[start:end], read_qual[start:end])
            if n_trimmed:
                start += n_trimmed
                trimmed_bases_left.append(n_trimmed)
                trimmed_names_left.append(trimmer.name)

                self.stats[f"N_{trimmer.name}_trimmed"] += 1
                self.total[f"bp_{trimmer.name}_trimmed"] += n_trimmed
                self.total[f"bp_trimmed"] += n_trimmed

            if check_discard(start, end, f"N_too_short_after_{trimmer.name}"):
                return

        # we've made it to the end!
        self.stats["N_kept"] += 1
        self.total["bp_kept"] += end

        tags = []
        if trimmed_names_right:
            tags.append(f"A3:Z:{','.join(trimmed_names_right)}")
            tags.append(f"T3:Z:{','.join([str(s) for s in trimmed_bases_right])}")

        if trimmed_names_left:
            tags.append(f"A5:Z:{','.join(trimmed_names_left)}")
            tags.append(f"T5:Z:{','.join([str(s) for s in trimmed_bases_left])}")

        return start, end, "\t".join(tags)


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
    flavor = AdapterFlavor(args, stats=stats, total=total)
    # TODO: paired-end processing
    for read in read_source:
        # read_seq = read.query_sequence
        # read_qual = read.query_qualities
        # we got a string
        cols = read.split("\t")
        read_seq = cols[9]
        qual_str = cols[10]
        read_qual = np.array(bytearray(qual_str.encode("ASCII"))) - args.phred_base

        result = flavor.process_read(read_seq, read_qual)
        if result:
            start, end, tags = result
            lhist[end - start] += 1
            cols[9] = read_seq[start:end]
            cols[10] = qual_str[start:end]

            if tags:
                cols[-1] = f"{cols[-1]} {tags}"

            yield "\t".join(cols)


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
    logger = util.setup_logging(args, name="spacemake.cutadapt_bam.main_single")

    bam_in = util.quiet_bam_open(
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
    for read in string_to_BAM(
        process_reads(
            BAM_to_string(skim_reads(bam_in.fetch(until_eof=True), args.skim)),
            args,
            stats=stats,
            total=total,
            lhist=lhist,
        )
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
        bam_in = util.quiet_bam_open(
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
    logger = util.setup_logging(args, name="spacemake.cutadapt_bam.main_parallel")

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
