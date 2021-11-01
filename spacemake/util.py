import errno
import os
import logging

from contextlib import ContextDecorator, contextmanager
from spacemake.errors import SpacemakeError, FileWrongExtensionError

LINE_SEPARATOR = "-" * 50 + "\n"

bool_in_str = ["True", "true", "False", "false"]

def assert_file(file_path, default_value=None, extension=["all"]):
    if file_path == default_value:
        # file doesn't exist but has the default value,
        # so we do not need to assert anything
        return False

    if isinstance(extension, str):
        extension = [extension]

    if not isinstance(file_path, list):
        file_path = [file_path]

    for fp in file_path:
        # check if file exists, raise error if not
        if not os.path.isfile(fp):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), fp)

        # check for all extensions
        if extension != ["all"]:
            extension_check = [fp.endswith(ex) for ex in extension]
            if not any(extension_check):
                raise FileWrongExtensionError(fp, extension)

    # return true if file exists and every test was good
    return True


def str2bool(var):
    if isinstance(var, bool):
        return var

    if var in ["True", "true"]:
        return True
    elif var in ["False", "false"]:
        return False
    else:
        raise ValueError(f"variable should be boolean, or one of: {bool_in_str}")


def ensure_path(path):
    import os


def ensure_path(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path


def FASTQ_src(src):
    from more_itertools import grouper

    for name, seq, _, qual in grouper(src, 4):
        yield name.rstrip()[1:], seq.rstrip(), qual.rstrip()


def BAM_src(src):
    import pysam

    bam = pysam.AlignmentFile(src, "rb", check_sq=False)
    for read in bam.fetch(until_eof=True):
        yield read.query_name, read.query_sequence, read.query_qualities


def read_fq(fname):
    import gzip

    if fname.endswith(".gz"):
        src = FASTQ_src(gzip.open(fname, mode="rt"))
    elif fname.endswith(".bam"):
        src = BAM_src(fname)
    elif type(fname) is str:
        src = FASTQ_src(open(fname))
    else:
        src = FASTQ_src(fname)  # assume its a stream or file-like object already

    for name, seq, qual in src:
        yield name, seq, qual


def dge_to_sparse(dge_path):
    import anndata
    import numpy as np
    import gzip
    import pandas as pd
    from scipy.sparse import csr_matrix, vstack

    ix = 0
    mix = 0
    matrices = []
    gene_names = []

    with gzip.open(dge_path, "rt") as dge:
        first_line = dge.readline().strip().split("\t")

        barcodes = first_line[1:]
        ncol = len(barcodes)
        M = np.zeros((1000, ncol))

        for line in dge:
            vals = line.strip().split("\t")
            # first element contains the gene name
            gene_names.append(vals[0])
            vals = vals[1:]
            vals = np.array(vals, dtype=np.int64)

            M[ix] = vals

            if ix % 1000 == 999:
                print(mix)
                mix = mix + 1
                matrices.append(csr_matrix(M))
                ix = 0
                M = np.zeros((1000, ncol))
            else:
                ix = ix + 1

        # get the leftovers
        M = M[:ix]
        matrices.append(csr_matrix(M))

        # sparse expression matrix
        X = vstack(matrices, format="csr")
        print(len(gene_names))

        adata = anndata.AnnData(
            X.T, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=gene_names)
        )

        return adata


def compute_neighbors(adata, min_dist=None, max_dist=None):
    """Compute all direct neighbors of all spots in the adata. Currently
        tailored for 10X Visium.
    Args:
        adata: an AnnData object
        min_dist: int, minimum distance to consider neighbors
        max_dist: int, maximum distance to consider neighbors
    Returns:
        neighbors: a dictionary holding the spot IDs for every spot
    """
    import numpy as np
    from scipy.spatial.distance import cdist

    # Calculate all Euclidean distances on the adata
    dist_mtrx = cdist(adata.obsm["spatial"], adata.obsm["spatial"])
    neighbors = dict()

    for i in range(adata.obs.shape[0]):
        neighbors[i] = np.where(
            (dist_mtrx[i, :] > min_dist) & (dist_mtrx[i, :] < max_dist)
        )[0]

    return neighbors


def compute_islands(adata, min_umi):
    """Find contiguity islands

    Args:
        adata: an AnnData object
        cluster: a cell type label to compute the islands for
    Returns:
        islands: A list of lists of spots forming contiguity islands
    """
    import numpy as np
    import itertools

    # this is hard coded for now for visium, to have 6 neighbors per spot
    # TODO: define an iterative approach where the key is to have around 6
    # neighbors per spot on average
    neighbors = compute_neighbors(adata, min_dist=0, max_dist=3)
    spots_cluster = np.where(np.array(adata.obs["total_counts"]) < min_umi)[0]

    # create a new dict holding only neighbors of the cluster
    islands = []

    for spot in neighbors:
        if spot not in spots_cluster:
            continue
        islands.append({spot}.union({x for x in neighbors[spot] if x in spots_cluster}))

    # merge islands with common spots
    island_spots = set(itertools.chain.from_iterable(islands))

    for each in island_spots:
        components = [x for x in islands if each in x]
        for i in components:
            islands.remove(i)
        islands += [list(set(itertools.chain.from_iterable(components)))]

    return islands


def detect_tissue(adata, min_umi):
    """Detect tissue: first find beads with at least min_umi UMIs, then detect island in the rest

    Args
        adata: an AnnData object, with spatial coordinates
        min_umi: integer, the min umi to be assigned as tissue bead by default
    Returns:
        tissue_indices: a list of indices which should be kept for this AnnData object
    """
    import numpy as np

    islands = compute_islands(adata, min_umi)

    # find the sizes of the islands. remove the biggest, as either the tissue has a big hole in it
    # or there are not so many big islands in which case removal is OK.
    # to be evaluated later...
    island_sizes = [len(island) for island in islands]

    tissue_islands = np.delete(islands, np.argmax(island_sizes))

    # get the indices of the islands
    tissue_indices = np.where(np.array(adata.obs["total_counts"]) >= min_umi)[0]

    tissue_indices = np.append(tissue_indices, np.hstack(tissue_islands))

    return tissue_indices


COMPLEMENT = {
    "a": "t",
    "t": "a",
    "c": "g",
    "g": "c",
    "k": "m",
    "m": "k",
    "r": "y",
    "y": "r",
    "s": "s",
    "w": "w",
    "b": "v",
    "v": "b",
    "h": "d",
    "d": "h",
    "n": "n",
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "K": "M",
    "M": "K",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "B": "V",
    "V": "B",
    "H": "D",
    "D": "H",
    "N": "N",
    "-": "-",
    "=": "=",
    "+": "+",
}


def complement(s):
    return "".join([COMPLEMENT[x] for x in s])


def rev_comp(seq):
    return complement(seq)[::-1]


def fasta_chunks(lines, strip=True, fuse=True):
    chunk = ""
    data = []

    for l in lines:
        if l.startswith("#"):
            continue
        if l.startswith(">"):
            if data and chunk:
                # print chunk
                yield chunk, "".join(data)

                if strip:
                    data = []
                else:
                    data = [l]

            chunk = l[1:].strip()
        else:
            if fuse:
                data.append(l.rstrip())
            else:
                data.append(l)

    if data and chunk:
        yield chunk, "".join(data)


@contextmanager
def message_aggregation(
    log_listen="spacemake",
    print_logger=False,
    print_success=True
):
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
            msg = f"{log_listen}: ".join([m + "\n" for m in message_buffer])
        else:
            msg = "\n".join(message_buffer)

        if print_success:
            msg = f"{msg}\n{LINE_SEPARATOR}SUCCESS!"

        print(msg)

    except SpacemakeError as e:
        print(e)


def str_to_list(value):
    # if list in string representation, return the list
    if value is None:
        return None

    if type(value) is str and value.startswith("["):
        return eval(value)
    # else create a list
    else:
        return [value]

def check_star_index_compatibility(star_index_dir):
    import os
    
    star_version = os.popen('STAR --version').read().strip()

    with open(os.path.join(star_index_dir, 'Log.out'), 'r') as f:
        first_line = f.readline().strip()
        index_version = first_line.split('=')[-1].split('_')[-1]

        if index_version != star_version:
            raise SpacemakeError(f'STAR index version ({index_version}) is' +
                f' incompatible with your STAR version ({star_version})')
