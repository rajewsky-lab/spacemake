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


def queue_iter(Q, abort_flag, stop_item=None, timeout=1):
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
            if item == stop_item:
                # signals end->exit
                break
            else:
                # logging.debug(f"queue_iter->item {item}")
                yield item


def join_with_empty_queues(proc, Qs, abort_flag, timeout=1):
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
        if abort_flag.value:
            for Q, content in zip(Qs, contents):
                content.extend(drain(Q))

    return contents


def order_results(res_queue, abort_flag, logger):
    import heapq
    import time

    heap = []
    n_chunk_needed = 0
    t0 = time.time()
    t1 = t0
    n_rec = 0

    for n_chunk, results in queue_iter(res_queue, abort_flag):
        heapq.heappush(heap, (n_chunk, results))

        # as long as the root of the heap is the next needed chunk
        # pass results on to storage
        while heap and (heap[0][0] == n_chunk_needed):
            n_chunk, results = heapq.heappop(heap)  # retrieves heap[0]
            for record in results:
                # print("record in process_ordered_results", record)
                yield record
                n_rec += 1

            n_chunk_needed += 1

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
