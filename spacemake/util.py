import errno
import os
import logging

from contextlib import ContextDecorator, contextmanager
from spacemake.errors import SpacemakeError, FileWrongExtensionError

LINE_SEPARATOR = '-'*50+'\n'

bool_in_str = ['True', 'true', 'False', 'false']

def assert_file(file_path, default_value='none', extension='all'):
    if file_path == default_value:
        # file doesn't exist but has the default value,
        # so we do not need to assert anything
        return False

    # check if file exists, raise error if not
    if not os.path.isfile(file_path):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), file_path)


    if not file_path.endswith(extension) and extension != 'all':
        raise FileWrongExtensionError(file_path, extension)

    # return true if file exists and every test was good
    return True

def str2bool(var):
    if isinstance(var, bool):
        return var

    if var in ['True', 'true']:
        return True
    elif var in ['False', 'false']:
        return False
    else:
        raise ValueError(f'variable should be boolean, or one of: {bool_in_str}')

def ensure_path(path):
    import os

    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path


def read_fq(fname):
    import gzip
    from more_itertools import grouper

    if fname.endswith(".gz"):
        src = gzip.open(fname, mode="rt")
    elif type(fname) is str:
        src = open(fname)
    else:
        src = fname  # assume its a stream or file-like object already

    for name, seq, _, qual in grouper(src, 4):
        yield name.rstrip()[1:], seq.rstrip(), qual.rstrip()

@contextmanager
def message_aggregation(log_listen='spacemake', print_logger=False):
    message_buffer = []

    log = logging.getLogger(log_listen)
    log.setLevel(logging.INFO)

    class MessageHandler(logging.NullHandler):
        def handle(this, record):
            if record.name == log_listen:
                message_buffer.append(record.msg)

    log.addHandler(MessageHandler())

    try:
        yield True
        if print_logger:
            msg = f'{log_listen}: '.join([m + '\n' for m in message_buffer])
        else:
            msg = '\n'.join(message_buffer)

        msg = f'{msg}\n{LINE_SEPARATOR}SUCCESS!'

        print(msg)

    except SpacemakeError as e:
        print(e)
