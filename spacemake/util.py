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
