from spacemake.contrib import __version__, __author__, __license__, __email__
import logging
import time


def put_or_abort(Q, item, abort_flag, timeout=1):
    """
    Small wrapper around queue.put() to prevent
    dead-locks in the event of (detectable) errors
    that might cause put() to block forever.
    Expects a shared mp.Value instance as abort_flag

    Returns: False if put() was succesful, True if execution
    should be aborted.
    """
    import queue

    sent = False
    # logging.warning(f"sent={sent} abort_flag={abort_flag}")
    while not (sent or abort_flag.value):
        try:
            Q.put(item, timeout=timeout)
        except queue.Full:
            time.sleep(0.1)
        else:
            sent = True

    return abort_flag.value


def queue_iter(Q, abort_flag, stop_item=None, timeout=1, logger=logging):
    """
    Small generator/wrapper around multiprocessing.Queue allowing simple
    for-loop semantics:

        for item in queue_iter(queue, abort_flag):
            ...
    The abort_flag is handled analogous to put_or_abort, only
    that it ends the iteration instead
    """
    import queue

    # logging.debug(f"queue_iter({queue})")
    while True:
        if abort_flag.value:
            break
        try:
            item = Q.get(timeout=timeout)
        except queue.Empty:
            pass
        else:
            if item is stop_item:
                logger.debug("received stop_item")
                # signals end->exit
                try:
                    queue.task_done()
                except AttributeError:
                    # We do not have a JoinableQueue. Fine, no problem.
                    pass

                break
            else:
                # logging.debug(f"queue_iter->item {item}, not {stop_item}")
                # print(f"queue_iter->item {type(item)}, not {type(stop_item)}")
                yield item


def iter_queues_round_robin(Qs, abort_flag, stop_item=None, timeout=1, logger=logging):
    """
    Small generator/wrapper around multiprocessing.Queue allowing simple
    for-loop semantics over a list of Queues that contain data in a round-robin
    fashion:

        for item in iter_queues_round_robin(queues, abort_flag):
            ...
    The abort_flag is handled analogous to put_or_abort, only
    that it ends the iteration instead
    """
    import queue

    # logging.debug(f"queue_iter({queue})")
    n = 0
    N = len(Qs)
    while True:
        # which queue do we need to poll
        i = n % N
        if abort_flag.value:
            break
        try:
            item = Qs[i].get(timeout=timeout)
        except queue.Empty:
            pass
        else:
            if item == stop_item:
                logger.debug("received stop_item")
                # signals end->exit
                try:
                    queue.task_done()
                except AttributeError:
                    # We do not have a JoinableQueue. Fine, no problem.
                    pass

                break
            else:
                # logging.debug(f"queue_iter->item {type(item)}")
                yield item
                n += 1


def join_with_empty_queues(proc, Qs, abort_flag, timeout=1, logger=logging):
    """
    joins() a process that writes data to queues Qs w/o deadlock.
    In case of an abort, the subprocess normally would not join
    until the Qs are emptied. join_with_empty_queues() monitors a global
    abort flag and empties the queues if needed, allowing the sub-process
    to terminate properly.
    """

    def drain(Q):
        content = []
        while not Q.empty():
            try:
                item = Q.get(timeout=timeout)
            except queue.Empty:
                pass
            else:
                content.append(item)

        return content

    contents = [list() for i in range(len(Qs))]
    while proc.exitcode is None:
        proc.join(timeout)
        logger.debug(f"Q sizes: {[Q.qsize() for Q in Qs]}")
        if abort_flag.value:
            for Q, content in zip(Qs, contents):
                content.extend(drain(Q))

    return contents


def order_results(res_queue, abort_flag, logger=logging):
    import heapq
    import time

    heap = []
    n_chunk_needed = 0
    t0 = time.time()
    t1 = t0
    n_rec = 0

    for n_chunk, results in queue_iter(res_queue, abort_flag):
        logger.debug(f"received chunk {n_chunk}")
        heapq.heappush(heap, (n_chunk, results))

        # as long as the root of the heap is the next needed chunk
        # pass results on to storage
        while heap and (heap[0][0] == n_chunk_needed):
            n_chunk, results = heapq.heappop(heap)  # retrieves heap[0]
            logger.debug(f"picked {n_chunk} from heap")
            for record in results:
                # print("record in process_ordered_results", record)
                yield record
                n_rec += 1

            n_chunk_needed = n_chunk + 1

        # debug output on average throughput
        t2 = time.time()
        if t2 - t1 > 30:
            dT = t2 - t0
            rate = n_rec / dT
            logger.debug(
                "processed {0} records in {1:.0f} seconds (average {2:.0f} reads/second).".format(
                    n_rec, dT, rate
                )
            )
            t1 = t2

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
        "finished processing {0} records in {1:.0f} seconds (average {2:.0f} reads/second)".format(
            n_rec, dT, n_rec / dT
        )
    )


def chunkify(src, n_chunk=1000):
    """
    Iterator which collects up to n_chunk items from iterable <src> and yields them
    as a list.
    """
    chunk = []
    n = 0
    for x in src:
        chunk.append(x)
        if len(chunk) >= n_chunk:
            yield n, chunk
            n += 1
            chunk = []

    if chunk:
        yield n, chunk


def log_qerr(qerr):
    "helper function for reporting errors in sub processes"
    for name, lines in qerr:
        for line in lines:
            logging.error(f"subprocess {name} exception {line}")


class ExceptionLogging:
    """
    A context manager that handles otherwise uncaught exceptions by logging
    the event and traceback info, optinally raises a flag.
    Very handy for wrapping the main function in a sub-process!
    """

    def __init__(self, name, Qerr=None, exc_flag=None, set_proc_title=True):
        # print('__init__ called')
        self.Qerr = Qerr
        self.exc_flag = exc_flag
        self.name = name
        self.logger = logging.getLogger(name.split()[0])
        self.exception = None

        if set_proc_title:
            import setproctitle

            setproctitle.setproctitle(name)

    def __enter__(self):
        self.t0 = time.time()
        # print('__enter__ called')
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        # print('__exit__ called')
        self.t1 = time.time()
        self.logger.debug(f"CPU time: {self.t1 - self.t0:.3f} seconds.")
        if exc_type and (exc_type != SystemExit):
            import traceback

            lines = "\n".join(
                traceback.format_exception(exc_type, exc_value, exc_traceback)
            ).split("\n")
            self.exception = lines
            self.logger.error(f"an unhandled exception occurred")
            for l in lines:
                self.logger.error(l)

            if self.Qerr is not None:
                self.Qerr.put((self.name, lines))

            if self.exc_flag is not None:
                self.logger.error(f"raising exception flag {self.exc_flag}")
                self.exc_flag.value = True


class SplitBAMbyUMI:
    def __init__(self, Qin, buf_size=10000, k=4):
        import re
        self.logger = logging.getLogger("spacemake.quant.SplitBAMbyUMI")
        self.Qin = Qin
        self.n = len(Qin)
        self.k = k
        self.regexp = re.compile(f'MI:Z:({"." * self.k})')

        self.buffers = [list() for n in range(self.n)]
        self.buffer_map = {}

        from spacemake.util import generate_kmers
        for i, kmer in enumerate(generate_kmers(k, nts='ACGTN')):
            j = i % self.n
            self.buffer_map[kmer] = self.buffers[j]
            self.logger.debug(f"assigning {kmer} -> worker {j}")

        self.logger.debug(f"splitting reads by first 4 bases of the UMI to assign to {self.n} queues")
        self.max_buf_size = buf_size
        self.n_chunks = 0

    def push_to_queues(self, bam_name):
        import numpy as np
        sizes = []
        self.logger.debug(f"push_to_queues, {bam_name}, buffer_sizes={self.get_buffer_sizes()}")
        for buf, Q in zip(self.buffers, self.Qin):
            sizes.append(len(buf))
            if buf:
                # print("pushing buffer", buf)
                # MJ: just spent the better part of an hour debugging 
                # to realize there is some async black magic happening 
                # in Q.put() or whatever.
                # cause if you send the actual buf and then clear() it
                # right after put() the receiver gets truncated or no data.
                # Yeah. For. real.
                # so buf.copy() it is... :-/
                Q.put((self.n_chunks, bam_name, buf.copy()))
                buf.clear()
                self.n_chunks += 1

        # self.logger.debug(f"push_to_queues done, {bam_name}, buffer_sizes={self.get_buffer_sizes()}")
        return np.array(sizes)
  
    def broadcast(self, rec):
        # print(f"broadcast, buffers={self.buffers}")
        for buf in self.buffers:
            buf.append(rec)

    def shutdown(self):
        self.logger.info("sending None to all queues to signal end of input")
        for Q in self.Qin:
            Q.put(None)

    def iter_bam(self, bam_path):
        import subprocess
        import re
        import os

        bam_name = os.path.basename(bam_path)
        proc = subprocess.Popen(['samtools', 'view', '-h', '--no-PG', '--threads=4', bam_path], stdout=subprocess.PIPE, text=True)
        while line := proc.stdout.readline():
            if line.startswith('@'):
                self.broadcast(line)
            else:
                M = re.search(self.regexp, line)
                buf = self.buffer_map[M.groups()[0]]
                # print(f"appending to buffer of L={len(buf)}")
                buf.append(line)

                if len(buf) > self.max_buf_size:
                    yield self.push_to_queues(bam_name)
        
        self.logger.debug(f"Completed ingesting from {bam_name}. Flushing buffers and moving on")
        yield self.push_to_queues(bam_name)

    def get_buffer_sizes(self):
        return [len(buf) for buf in self.buffers]

class CountingStatistics:
    def __init__(self):
        from collections import defaultdict
        # lambda functions can not be pickled! Hence the copy
        self.stats_by_ref = defaultdict(defaultdict(int).copy)

    def count(self, ref, name):
        self.stats_by_ref[ref][name] += 1

    def add_other_stats(self, other):
        for ref, stats in other.items():
            for k, v in stats.items():
                self.stats_by_ref[ref][k] += v

    def get_stats_df(self):
        import pandas as pd
        dfs = []
        for ref in sorted(self.stats_by_ref.keys()):
            data = {}
            data['ref'] = ref
            for k, v in self.stats_by_ref[ref].items():
                data[k] = [v]

            df = pd.DataFrame(data)
            dfs.append(df[sorted(df.columns)])

        if not dfs:
            dfs = [pd.DataFrame({})]

        return pd.concat(dfs)

    def save_stats(self, path, sep='\t', **kw):
        self.get_stats_df().to_csv(path, sep=sep, **kw)



def BAM_reader(Qin, bam_in, Qerr, abort_flag, stat_list, buffer_size=10000, log_domain="spacemake.parallel", counter_class=CountingStatistics):
    import os
    with ExceptionLogging(
        f"{log_domain}.BAM_reader", Qerr=Qerr, exc_flag=abort_flag
    ) as el:
        from time import time
        stats = counter_class() # statistics on disambiguation and counting
        ref_stats = None
        last_ref = None
        dispatch = SplitBAMbyUMI(Qin, buf_size=buffer_size)

        # iterate over all BAMs with countable alignments
        N = 0
        T0 = time()
        for bam_name in bam_in:
            reference_name = os.path.basename(bam_name).split(".")[0]
            el.logger.info(f"reading alignments from {bam_name} to reference '{reference_name}'.")
            t0 = time()
            t1 = time()
            N_ref = 0

            # split-by-UMI dispatch
            for n_pushed in dispatch.iter_bam(bam_name):
                n = n_pushed.sum()
                N += n
                N_ref += n
                # ref_stats['N_records'] = N_ref
                if time() - t1 > 2:
                    t = time() - t0
                    el.logger.debug(f"worker load distribution: {n_pushed / float(n)}")
                    el.logger.info(f"ingested {N_ref} BAM records in {t:.1f} seconds ({0.001 * N_ref/t:.2f} k/sec).")
                    t1 = time()

        t = time() - T0
        el.logger.debug(f"worker load distribution: {n_pushed / float(n)}")
        el.logger.info(f"ingested {N} BAM records in {t:.1f} seconds ({0.001 * N/t:.2f} k/sec).")
        t1 = time()

        el.logger.debug("done. closing queues.")
        dispatch.shutdown()
        
        el.logger.debug("syncing stats...")
        stat_list.append(stats.stats_by_ref)
        
        el.logger.debug("shutting down...")


class AlignedSegmentsFromQueue:
    def __init__(self, Qin, abort_flag):
        self.Qin = Qin
        self.abort_flag = abort_flag
        self.header = None
        self.header_lines = []
        self.last_ref = None
        self.logger = logging.getLogger("spacemake.parallel.AlignedSegmentsFromQueue")

    def __iter__(self):
        import pysam
        for n_chunk, ref, sam_lines in queue_iter(self.Qin, self.abort_flag):
            if ref != self.last_ref:
                self.logger.debug(f"switched BAM {self.last_ref} -> {ref}, expecting new header.")
                self.header = None
                self.header_lines = []
                self.last_ref = ref

            self.logger.debug(f"received chunk {n_chunk} ref={ref} with {len(sam_lines)} line 0: {sam_lines[0]}")

            for line in sam_lines:
                # print(line)
                if line.startswith('@'):
                    self.header_lines.append(line)
                else:
                    if not self.header:
                        # print("about to construct header:")
                        # print(self.header_lines)
                        self.header = pysam.AlignmentHeader.from_text("".join(self.header_lines))
                        # print(self.header)
                    else:
                        sam = pysam.AlignedSegment.fromstring(line, self.header)
                        # print(sam)
                        yield ref, sam
            
            self.Qin.task_done()


def parallel_BAM_workflow(bam_inputs, worker_func, collector_func, log_domain="spacemake.parallel", n_workers=8, counter_class=CountingStatistics, buffer_size=10000, **kw):
    import spacemake.util as util
    import multiprocessing as mp

    # queues for communication between processes
    Qin = [mp.JoinableQueue(n_workers * 10) for i in range(n_workers)]
    Qout = mp.JoinableQueue()
    Qerr = mp.Queue()  # child-processes can report errors back to the main process here

    def get_params(keys, scope={}, from_scratch=True):
        "used to gather kwargs for sub-processes"
        if from_scratch:
            d = {}
        else:
            d = kw.copy()

        for k in keys:
            d[k] = scope.get(k, None)

        return d

    # Proxy objects to allow workers to report statistics about the run
    with mp.Manager() as manager:
        abort_flag = mp.Value("b")
        abort_flag.value = False

        stat_list = manager.list()
        with ExceptionLogging(
            f"{log_domain}.main_parallel", exc_flag=abort_flag
        ) as el:

            # read BAM in chunks and put them in Qsam
            dispatcher = mp.Process(
                target=BAM_reader,
                name="BAM_reader",
                args=(Qin, bam_inputs, Qerr, abort_flag, stat_list),
                kwargs=get_params(['log_domain', 'buffer_size', 'counter_class'], scope=locals()),
            )

            dispatcher.start()
            el.logger.debug("Started dispatch")

            # workers consume chunks of BAM from Qsam
            # process them, and put the results in Qres
            workers = []
            for i in range(n_workers):
                w = mp.Process(
                    target=worker_func,
                    name=f"worker_{i}",
                    args=(Qin[i], Qout, Qerr, abort_flag, stat_list),
                    kwargs=get_params(['log_domain', 'buffer_size', 'counter_class'], scope=locals(), from_scratch=False)
                )
                w.start()
                workers.append(w)

            el.logger.debug("Started workers")

            collector = mp.Process(
                target=collector_func,
                name="output",
                args=(Qout, Qerr, abort_flag, stat_list),
                kwargs=get_params(['log_domain', 'buffer_size', 'counter_class'], scope=locals(), from_scratch=False)
            )
            collector.start()
            el.logger.debug("Started collector")
            # wait until all sequences have been thrown onto Qfq
            qerr = join_with_empty_queues(dispatcher, Qin + [Qerr], abort_flag, logger=el.logger)
            el.logger.debug("The dispatcher exited")

            if qerr[-1]:
                el.logger.info(f"{len(qerr)} chunks were drained from Qfq upon abort.")
                log_qerr(qerr[-1])

            # wait until the workers have completed every task and placed
            # every output chunk onto Qres.

            import time
            el.logger.debug("Waiting for workers to process all data")
            for Q, w in zip(Qin, workers):
                while Q.qsize() > 0:
                    time.sleep(.1)
                # el.logger.debug(f"Q-size: {Q.qsize()}")
                # w.join()
                # el.logger.debug(f"Q-size: {Q.qsize()}")
                # el.logger.debug(".")

            el.logger.debug("Waiting for writer to accumulate all data")
            # Qout.join()
            while Qout.qsize():
                time.sleep(.1)
                # print(Qout.qsize())

            el.logger.debug("Telling writer to stop.")
            Qout.put(None)


            collector.join()
            el.logger.debug("Collector has joined.")

            # el.logger.debug("Joining input queues")
            # [Q.join() for Q in Qin]

            el.logger.debug("Waiting for workers to exit")
            for i, w in enumerate(workers):
                w.join()

            el.logger.debug(
                "All worker processes have joined. Gathering statistics."
            )
            
            stats = counter_class()
            for s in stat_list:
                stats.add_other_stats(s)
                
        if el.exception:
            return None
    
        return stats
    


def batched(iterable, n):
    "Batch data into lists of length n. The last batch may be shorter."
    from itertools import islice
    # batched('ABCDEFG', 3) --> ABC DEF G
    it = iter(iterable)
    while True:
        batch = list(islice(it, n))
        if not batch:
            return
        yield batch


from contextlib import contextmanager

@contextmanager
def create_named_pipes(names):
    import os
    import tempfile
    with tempfile.TemporaryDirectory() as base:
        paths = [os.path.join(base, name) for name in names]
        #print(paths)
        # create new fifos (named pipes)
        [os.mkfifo(fname) for fname in paths]

        try:
            yield paths
        finally:
            # Clean up the named pipes
            for fname in paths:
                os.remove(fname)

def open_named_pipe(path, mode='rt+', buffer_size=1000000):
    import fcntl
    F_SETPIPE_SZ = 1031  # Linux 2.6.35+
    # F_GETPIPE_SZ = 1032  # Linux 2.6.35+

    fifo_fd = open(path, mode)
    fcntl.fcntl(fifo_fd, F_SETPIPE_SZ, buffer_size)
    return fifo_fd


def igzip_reader(input_files, pipe):
    with ExceptionLogging("spacemake.parallel.igzip_reader") as el:
        el.logger.info(f"writing to {pipe}")
        from isal import igzip
        out_file = open_named_pipe(pipe, mode='wb')

        try:
            for fname in input_files:
                el.logger.info(f"reading from {fname}")
                in_file = igzip.IGzipFile(fname, 'r')
                while True:
                    block = in_file.read(igzip.READ_BUFFER_SIZE)
                    if block == b"":
                        break

                    out_file.write(block)
                
                in_file.close()
        finally:
            el.logger.info(f"closing down {pipe}")
            out_file.close()

import pyximport; pyximport.install()
import spacemake.cython.fast_loop as fast_loop

def fastq_distributor(in_pipe, worker_in_pipes):
    with ExceptionLogging("spacemake.parallel.fastq_distributor") as el:
        # from .fast_loop import distribute
        el.logger.info(f"reading from {in_pipe}, writing to {worker_in_pipes} chunk_size=4")
        fast_loop.distribute(in_pipe, worker_in_pipes, chunk_size=80)


def fastq_to_sam_worker(w_in1, w_in2, w_out, args):
    with ExceptionLogging("spacemake.parallel.fastq_to_SAM_worker") as el:

        from spacemake.bin.fastq_to_uBAM import Formatter
        formatter = Formatter(args)

        fin1 = open_named_pipe(w_in1, mode='rt', buffer_size=2**19)
        fin2 = open_named_pipe(w_in2, mode='rt', buffer_size=2**19)
        fout = open_named_pipe(w_out, mode='wt', buffer_size=2**19)

        from spacemake.util import FASTQ_src
        for (qname, r1, r1_qual), (r2_qname, r2, r2_qual) in zip(FASTQ_src(fin1), FASTQ_src(fin2)):
            # el.logger.warning(f"{w_in1} {w_in2} -> got read. writing to {w_out}")
            #fout.write(f"{name1.rstrip()}\t{seq1.rstrip()}\t{seq2.rstrip()}\n")
            kw = formatter.format(flags=4, qname=qname, r1=r1, r1_qual=r1_qual, r2_qname=r2_qname, r2=r2, r2_qual=r2_qual)
            #fout.write(formatter.make_bam_record(**locals()))
            # QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [ALIGNMENT SECTION]
            fout.write("{qname}\t{flags}\t*\t0\t255\t*\t*\t0\t0\t{seq}\t{qual}\tRG:Z:{rg}\tCB:Z:{cell}\tUMI:Z:{UMI}\n".format(qname=qname, flags=4, rg='A', **kw))


def collector(worker_out_pipes, out_pipe, args):
    from time import time
    with ExceptionLogging("spacemake.parallel.collector") as el:
        # from spacemake.experiment.fast_loop import collect
        el.logger.warning(f"collecting from {worker_out_pipes} into {out_pipe}")
        
        # make a BAM header and write it to a temporary file 
        from spacemake.bin.fastq_to_uBAM import Formatter
        header=str(Formatter(args).bam_header).rstrip()

        print(f"header='{header}'")
        import os
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as fhead:
            fhead.write(header + '\n')
            fhead.close()
            t0 = time()
            n = fast_loop.collect(worker_out_pipes, out_pipe, chunk_size=20, header_fifo=fhead.name)
            dt = time() - t0
            el.logger.info(f"collected {n:,.0f} records in {dt:.2f} seconds. Rate of {n/dt/1000:.2f}k rec/sec")
            os.unlink(fhead.name)


def fastq_to_ubam_main(input_files1, input_files2, args, n_workers=32):
    import multiprocessing as mp
    
    pipe_names = (['igzip1', 'igzip2'] + 
        [f'in_1_{n}' for n in range(n_workers)] + 
        [f'in_2_{n}' for n in range(n_workers)] + 
        [f'out_{n}' for n in range(n_workers)] +
        ['out_sam'])


    with create_named_pipes(pipe_names) as pipe_paths:
        igzip_pipe1, igzip_pipe2 = pipe_paths[:2]
        worker_inputs1 = pipe_paths[2:2+n_workers]
        worker_inputs2 = pipe_paths[2+n_workers:2+2*n_workers]
        worker_outputs = pipe_paths[2+2*n_workers:-1]

        p_gzip1 = mp.Process(target=igzip_reader, args=(input_files1, igzip_pipe1,))
        p_gzip2 = mp.Process(target=igzip_reader, args=(input_files2, igzip_pipe2,))
                                   
        p_dist1 = mp.Process(target=fastq_distributor, args=(igzip_pipe1, worker_inputs1))
        p_dist2 = mp.Process(target=fastq_distributor, args=(igzip_pipe2, worker_inputs2))

        workers = []
        for w_in1, w_in2, w_out in zip(worker_inputs1, worker_inputs2, worker_outputs):
            p = mp.Process(target=fastq_to_sam_worker, args=(w_in1, w_in2, w_out, args))
            workers.append(p)

        p_collect = mp.Process(target=collector, args=(worker_outputs, '/tmp/out', args))
        p_collect.start()

        for w in workers:
            w.start()
        
        p_dist1.start()
        p_dist2.start()
        p_gzip1.start()
        p_gzip2.start()

        p_gzip1.join()
        p_gzip2.join()
        p_dist1.join()
        p_dist2.join()

        for w in workers:
            w.join()
        
        p_collect.join()

    
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    import spacemake.util as util
    from spacemake.bin.fastq_to_uBAM import parse_args
    args = parse_args()
    util.setup_logging(args, name="spacemake.bin.fastq_to_uBAM")

    fastq_to_ubam_main(
        ['/data/rajewsky/sequencing/human/sc_smRNA_dss_69/dss_69-1-mRNA_S1_R1_001.fastq.gz'], 
        ['/data/rajewsky/sequencing/human/sc_smRNA_dss_69/dss_69-1-mRNA_S1_R2_001.fastq.gz'],
        args
    )
