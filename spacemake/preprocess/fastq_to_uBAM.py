#!/usr/bin/env python3
from spacemake.contrib import __version__, __author__, __license__, __email__
import logging
import os
import sys
import numpy as np
import pysam
import multiprocessing as mp
from collections import defaultdict

from spacemake.parallel import (
    put_or_abort,
    chunkify,
    queue_iter,
    order_results,
    join_with_empty_queues,
    ExceptionLogging,
    log_qerr,
)
import spacemake.util as util

# TODO:
# * turn Output into a context manager
# * load params from YAML (code from opseq)
# * save params to YAML in run-folder to document
# * refactor into multiple modules?


def dispatch_fastq_to_queue(Qfq, args, Qerr, abort_flag):
    """
    reads from two fastq files, groups the input into chunks for
    faster parallel processing, and puts these on a mp.Queue()
    """

    def read_source(args):
        for (id1, seq1, qual1), (id2, seq2, qual2) in zip(
            util.read_fq(args.read1), util.read_fq(args.read2)
        ):
            # assert id1.split()[0] == id2.split()[0]
            yield id1, seq1, id2, seq2, qual2

    with ExceptionLogging(
        "spacemake.preprocess.fastq.dispatcher", Qerr=Qerr, exc_flag=abort_flag
    ) as el:
        for chunk in chunkify(read_source(args), n_chunk=args.chunk_size):
            logging.debug(f"placing {chunk[0]} {len(chunk[1])} in queue")
            if put_or_abort(Qfq, chunk, abort_flag):
                el.logger.warning("shutdown flag was raised!")
                break


format_func_template = """
def format_func(qname=None, r2_qname=None, r2_qual=None, r1=None, r2=None):
    qparts = r2_qname.split(':N:0:')
    if len(qparts) > 1:
        i5i7 = [1].replace('+', '')
    else:
        i5i7 = None

    cell = {args.cell}
    raw = cell
    UMI = {args.UMI}
    seq = {args.seq}
    qual = {args.qual}

    return dict(cell=cell, raw=raw, UMI=UMI, seq=seq, qual=qual)
"""


class Output:
    logger = logging.getLogger("spacemake.preprocess.fastq.Output")

    def __init__(self, args, open_files=True):
        # assert Output.safety_check_eval(args.cell_raw)
        assert Output.safety_check_eval(args.cell)
        assert Output.safety_check_eval(args.UMI)
        assert Output.safety_check_eval(args.seq)
        assert Output.safety_check_eval(args.qual)

        # This trickery here compiles a python function that
        # embeds the user-specified code for cell, UMI
        # we replace globals() with an empty dictionary d,
        d = {}
        exec(format_func_template.format(args=args), d)
        # from said dict we can now retrieve the function object
        # and store it as self.format
        self.format = d["format_func"]

        self.fq_qual = args.fq_qual
        self.raw_cb_counts = defaultdict(int)
        self.count_cb = bool(args.save_cell_barcodes)

        self.tags = []
        for tag in args.bam_tags.split(","):
            self.tags.append(tag.split(":"))

        prog = os.path.basename(__file__)
        header = {
            "HD": {"VN": "1.6"},
            "PG": [
                {
                    "ID": "fastq_to_uBAM",
                    "PN": prog,
                    "CL": " ".join(sys.argv[1:]),
                    "VN": __version__,
                },
            ],
            "RG": [
                {"ID": "A", "SM": args.sample},
                # {"ID": "U", "SM": f"unassigned_{args.sample}"},
            ],
        }
        self.bam_header = pysam.AlignmentHeader.from_dict(header)
        fopen = lambda x: pysam.AlignmentFile(x, "wbu", header=header, threads=4)

        if open_files:
            self.out_bam = fopen(args.out_bam)

    @staticmethod
    def safety_check_eval(s, danger="();."):
        chars = set(list(s))
        if chars & set(list(danger)):
            return False
        else:
            return True

    def make_bam_record(self, flag=4, **kw):

        # pass what we already know about the read into format()
        # which contains the cmdline-specified code for cell and UMI.
        # since it returns a dictionary, lets just update kw and
        # keep using that to construct the BAM record.
        try:
            kw.update(self.format(**kw))
        except Exception as E:
            self.logger.error(
                "Unhandled exception inside user-defined read handler.\n"
                f"Input: {kw}\n"
                f"Error: {E}"
            )
            raise E

        a = pysam.AlignedSegment(self.bam_header)

        # STAR does not like spaces in read names so we have to split
        a.query_name = kw["r2_qname"].split()[0]
        a.query_sequence = kw["seq"]
        a.flag = flag
        a.query_qualities = pysam.qualitystring_to_array(kw["qual"])
        a.tags = [(name, templ.format(**kw)) for name, templ in self.tags]
        if self.count_cb:
            self.raw_cb_counts[kw["cell"]] += 1

        return a.to_string()

    def write_bam(self, rec):
        self.out_bam.write(pysam.AlignedSegment.fromstring(rec, self.bam_header))

    def close(self):
        self.out_bam.close()


def collect_from_queue(res_queue, args, Qerr, abort_flag):
    with ExceptionLogging(
        "spacemake.preprocess.fastq.collector", Qerr=Qerr, exc_flag=abort_flag
    ) as el:
        out = Output(args)

        for record in order_results(res_queue, abort_flag, el.logger):
            # el.logger.debug(f"got record: {record}")
            out.write_bam(record)

        out.close()


def parallel_worker(Qfq, Qres, args, Qerr, abort_flag, stat_lists):
    def quality_trim_read2(reads, min_qual=20, phred_base=33, min_len=18):
        for (name1, seq1, name2, seq2, qual2) in reads:
            end = len(seq2)
            q2 = np.array(bytearray(qual2.encode("ASCII"))) - phred_base
            qtrim = q2 >= min_qual
            new_end = end - (qtrim[::-1]).argmax()

            # TODO: yield A3,T3 adapter-trimming tags
            if new_end != end:
                qual2 = qual2[:new_end]
                seq2 = seq2[:new_end]

            if len(seq2) >= min_len:
                yield (name1, seq1, name2, seq2, qual2)

    with ExceptionLogging(
        "spacemake.preprocess.fastq.worker", Qerr=Qerr, exc_flag=abort_flag
    ) as el:
        el.logger.debug(
            f"process_dropseq starting up with Qfq={Qfq}, Qres={Qres} and args={args}"
        )
        out = Output(args, open_files=False)
        N = defaultdict(int)
        for n_chunk, reads in queue_iter(Qfq, abort_flag):
            el.logger.debug(f"received chunk {n_chunk} of {len(reads)} reads")
            results = []
            if args.min_qual_trim:
                reads = quality_trim_read2(
                    reads, min_qual=args.min_qual_trim, phred_base=args.phred_base
                )
            if args.paired_end:
                for fqid, r1, fqid2, r2, qual2 in reads:
                    N["total"] += 1
                    rec1 = out.make_bam_record(
                        qname=fqid,
                        r1=r1,
                        r2=r1,
                        r2_qual=qual2,
                        r2_qname=fqid2,
                        flag=69,  # unmapped, paired, first in pair
                    )
                    rec2 = out.make_bam_record(
                        qname=fqid,
                        r1=r1,
                        r2=r2,
                        r2_qual=qual2,
                        r2_qname=fqid2,
                        flag=133,  # unmapped, paired, second in pair
                    )
                    results.append(rec1)
                    results.append(rec2)

            else:
                for fqid, r1, fqid2, r2, qual2 in reads:
                    N["total"] += 1
                    rec = out.make_bam_record(
                        qname=fqid,
                        r1=r1,
                        r2=r2,
                        r2_qual=qual2,
                        r2_qname=fqid2,
                    )
                    results.append(rec)

            Qres.put((n_chunk, results))

        # return our counts and observations
        # via these list proxies
        el.logger.debug("synchronizing cache and counts")
        Ns = stat_lists[0]
        Ns.append(N)

        cb_counts = stat_lists[1]
        cb_counts.append(out.raw_cb_counts)


def count_dict_sum(sources):
    dst = defaultdict(float)
    for src in sources:
        for k, v in src.items():
            dst[k] += v

    return dst


def dict_merge(sources):
    dst = {}
    for src in sources:
        dst.update(src)
    return dst


def main_parallel(args):
    # queues for communication between processes
    Qfq = mp.Queue(args.parallel * 5)
    Qres = mp.Queue(args.parallel * 5)  # BAM records as string
    Qerr = mp.Queue()  # child-processes can report errors back to the main process here

    # Proxy objects to allow workers to report statistics about the run
    manager = mp.Manager()
    abort_flag = mp.Value("b")
    abort_flag.value = False
    Ns = manager.list()
    cb_counts = manager.list()
    stat_lists = [Ns, cb_counts]
    workers = []

    res = 0
    with ExceptionLogging(
        "spacemake.preprocess.fastq.main_parallel", exc_flag=abort_flag
    ) as el:

        # read FASTQ in chunks and put them in Qfq
        dispatcher = mp.Process(
            target=dispatch_fastq_to_queue,
            name="dispatcher",
            args=(Qfq, args, Qerr, abort_flag),
        )

        dispatcher.start()
        el.logger.info("Started dispatch")

        # workers consume chunks of FASTQ from Qfq,
        # process them, and put the results in Qres
        for i in range(args.parallel):
            w = mp.Process(
                target=parallel_worker,
                name=f"worker_{i}",
                args=(Qfq, Qres, args, Qerr, abort_flag, stat_lists),
            )
            w.start()
            workers.append(w)

        el.logger.info("Started workers")
        collector = mp.Process(
            target=collect_from_queue,
            name="output",
            args=(Qres, args, Qerr, abort_flag),
        )
        collector.start()
        el.logger.info("Started collector")
        # wait until all sequences have been thrown onto Qfq
        qfq, qerr = join_with_empty_queues(dispatcher, [Qfq, Qerr], abort_flag)
        el.logger.info("The dispatcher exited")
        if qfq or qerr:
            el.logger.error("Some exception occurred! Tearing it all down")
            res = 1
            el.logger.info(f"{len(qfq)} chunks were drained from Qfq upon abort.")
            log_qerr(qerr)
            dispatcher.terminate()
            collector.terminate()
            for w in workers:
                w.terminate()
        else:
            # signal all workers to finish
            el.logger.info("Signalling all workers to finish")
            for n in range(args.parallel):
                Qfq.put(None)  # each worker consumes exactly one None

        for w in workers:
            # make sure all results are on Qres by waiting for
            # workers to exit. Or, empty queues if aborting.
            qres, qerr = join_with_empty_queues(w, [Qres, Qerr], abort_flag)
            if qres or qerr:
                res = 1
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

        N = count_dict_sum(Ns)
        if N["total"]:
            el.logger.info(f"Run completed, {N['total']} reads processed.")
        else:
            el.logger.error("No reads were processed!")
        # el.logger.debug(print(str(N.keys())[:200]))
        cb_count = count_dict_sum(cb_counts)
        if args.save_cell_barcodes:
            el.logger.info(
                f"writing {len(cb_count)} barcode counts to '{args.save_cell_barcodes}'"
            )
            import gzip

            with gzip.open(util.ensure_path(args.save_cell_barcodes), "wt") as f:
                f.write("cell_bc\traw_read_count\n")
                for bc, count in sorted(cb_count.items()):
                    f.write(f"{bc}\t{count}\n")

        if args.save_stats:
            with open(util.ensure_path(args.save_stats), "w") as f:
                for k, v in sorted(N.items()):
                    f.write(f"freq\t{k}\t{v}\t{100.0 * v/max(N['total'], 1):.2f}\n")

    if el.exception:
        res = 1
        for w in workers:
            w.terminate()
            w.join()

    return res


def parse_args():
    import spacemake.util as util

    parser = util.make_minimal_parser(
        "fastq.py",
        description="Convert raw reads1 and reads2 FASTQ into a single BAM file with cell barcode and UMI as BAM-tags",
    )

    parser.add_argument(
        "--read1",
        default=None,
        help="source from where to get read1 (FASTQ format)",
    )
    parser.add_argument(
        "--read2",
        default="/dev/stdin",
        help="source from where to get read2 (FASTQ format)",
        required=True,
    )
    parser.add_argument("--cell", default="r1[8:20][::-1]")
    parser.add_argument("--UMI", default="r1[0:8]")
    parser.add_argument("--seq", default="r2")
    parser.add_argument("--qual", default="r2_qual")
    parser.add_argument(
        "--paired-end",
        default=False,
        action="store_true",
        help="read1 and read2 are paired end mates and store both in the BAM",
    )
    parser.add_argument(
        "--out-bam",
        default="/dev/stdout",
        help="output for unaligned BAM records (default=/dev/stdout) ",
    )
    parser.add_argument(
        "--save-stats",
        default="preprocessing_stats.txt",
        help="store statistics in this file",
    )
    parser.add_argument(
        "--save-cell-barcodes",
        default="",
        help="store (raw) cell barcode counts in this file. Numbers add up to number of total raw reads.",
    )
    parser.add_argument(
        "--fq-qual",
        default="E",
        help="phred qual for assigned barcode bases in FASTQ output (default='E')",
    )
    parser.add_argument(
        "--min-qual-trim",
        default=0,
        type=int,
        help="clip low quality-bases from a read2 3' end (e.g. pre-detected adapter)",
    )
    parser.add_argument(
        "--min-len",
        default=18,
        type=int,
        help="minimum read2 length to keep after quality trimming (default=18)",
    )
    parser.add_argument(
        "--phred-base",
        default=33,
        type=int,
        help="phred quality base (default=33)",
    )
    parser.add_argument(
        "--parallel", default=1, type=int, help="how many processes to spawn"
    )
    parser.add_argument(
        "--chunk-size",
        default=10000,
        type=int,
        help="how many reads (mate pairs) are grouped together for parallel processing (default=10000)",
    )
    # SAM standard now supports CB=corrected cell barcode, CR=original cell barcode, and MI=molecule identifier/UMI
    parser.add_argument(
        "--bam-tags",
        default="CR:{cell},CB:{cell},MI:{UMI},RG:A",
        help="a template of comma-separated BAM tags to generate. Variables are replaced with extracted cell barcode, UMI etc.",
    )
    args = parser.parse_args()

    return args


def cmdline():
    args = parse_args()
    util.setup_logging(args, name="spacemake.preprocess")

    with ExceptionLogging("spacemake.preprocess.main") as el:
        if not args.read2:
            raise ValueError("bam output format requires --read2 parameter")

        res = main_parallel(args)
    if el.exception:
        res = 1

    el.logger.info(f"exit code={res}")
    return res


if __name__ == "__main__":
    import os, sys

    ## We catch errors in the workers, but somehow
    # multiprocessing masks our non-zero exit code.
    # I think this is a bug...
    # sys.stderr.write("running that shit\n")
    res = cmdline()
    # sys.stderr.write(f"calling sys.exit({res})\n")
    # os._exit(1)
    # sys.exit(1)
    sys.exit(res)
