__version__ = "0.9"
__author__ = ["Marvin Jens"]
__license__ = "GPL"
__email__ = ["marvin.jens@mdc-berlin.de"]

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
            pass
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


class ExceptionLogging:
    """
    A context manager that handles otherwise uncaught exceptions by logging
    the event and traceback info, optinally raises a flag.
    Very handy for wrapping the main function in a sub-process!
    """

    def __init__(self, name, Qerr=None, exc_flag=None):
        # print('__init__ called')
        self.Qerr = Qerr
        self.exc_flag = exc_flag
        self.name = name
        self.logger = logging.getLogger(name)
        self.exception = None

    def __enter__(self):
        self.t0 = time.time()
        # print('__enter__ called')
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        # print('__exit__ called')
        self.t1 = time.time()
        self.logger.info(f"CPU time: {self.t1 - self.t0:.3f} seconds.")
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

            if self.exc_flag:
                self.logger.error(f"raising exception flag {self.exc_flag}")
                self.exc_flag.value = True
