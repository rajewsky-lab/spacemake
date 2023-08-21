import errno
import os
import logging

from contextlib import ContextDecorator, contextmanager
from spacemake.errors import SpacemakeError, FileWrongExtensionError
from spacemake.contrib import __version__, __license__, __author__, __email__

LINE_SEPARATOR = "-" * 50 + "\n"

bool_in_str = ["True", "true", "False", "false"]

def generate_kmers(k, nts='ACGT'):
    if k == 0:
        yield ''
    elif k > 0:
        for x in nts:
            for mer in generate_kmers(k-1, nts=nts):
                yield x + mer

class dotdict(dict):
    """dot.notation access to dictionary attributes"""

    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __str__(self):
        buf = ["dotdict"]
        for k, v in self.items():
            if not k.startswith("__"):
                buf.append(f"  {k} = {v}")

        return "\n".join(buf)


def wc_fill(x, wc):
    """
    Versatile snakemake helper function that can render any string template used in the mapping and reports modules,
    either filling in from a Wildcards object, or from a dotdict.
    """
    return x.format(
        sample_id=getattr(wc, "sample_id", "NO_sample_id"),
        project_id=getattr(wc, "project_id", "NO_project_id"),
        species=getattr(wc, "species", "NO_species"),
        ref_name=getattr(wc, "ref_name", "NO_ref_name"),
        mapper=getattr(wc, "mapper", "NO_mapper"),
        link_name=getattr(wc, "link_name", "NO_link_name"),
        polyA_adapter_trimmed=getattr(wc, "polyA_adapter_trimmed", ""),
        data_root_type=getattr(wc, "data_root_type", "complete_data"),
    )

def quiet_bam_open(*argc, **kw):
    """_summary_

    This wrapper around pysam.AlignmentFile() simply silences warnings about missing BAM index etc.
    We don't care about the index and therefore these error messages are just spam.

    Returns:
        pysam.AlignmentFile: the sam/bam object as returned by pysam.
    """
    import pysam

    save = pysam.set_verbosity(0)
    bam = pysam.AlignmentFile(*argc, **kw)
    pysam.set_verbosity(save)
    return bam


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
    dirname = os.path.dirname(path)
    if dirname:
        os.makedirs(dirname, exist_ok=True)
    return path

def timed_loop(
    src,
    logger,
    T=5,
    chunk_size=10000,
    template="processed {i} BAM records in {dT:.1f}sec. ({rate:.3f} k rec/sec)",
    skim=0,
):
    from time import time

    t0 = time()
    t_last = t0
    i = 0
    for i, x in enumerate(src):
        if skim:
            if i % skim == 0:
                yield x
        else:
            yield x

        if i % chunk_size == 0:
            t = time()
            if t - t_last > T:
                dT = t - t0
                rate = i / dT / 1000.0
                logger.info(template.format(**locals()))
                t_last = t

    t = time()
    dT = t - t0
    rate = i / dT / 1000.0
    logger.info("Finished! " + template.format(**locals()))


def FASTQ_src(src):
    from more_itertools import grouper

    for name, seq, _, qual in grouper(src, 4):
        yield name.rstrip()[1:], seq.rstrip(), qual.rstrip()


def BAM_src(src):
    import pysam

    bam = quiet_bam_open(src, "rb", check_sq=False)
    for read in bam.fetch(until_eof=True):
        yield read.query_name, read.query_sequence, read.query_qualities


def read_fq(fname, skim=0):
    import gzip

    logger = logging.getLogger("spacemake.util.read_fq")
    if str(fname) == "None":
        logger.warning("yielding empty data forever")
        while True:
            yield ("no_qname", "no_seq", "no_qual")

    if "*" in fname:
        logger.warning("EXPERIMENTAL: fname contains wildcards")
        from glob import glob

        for match in sorted(glob(fname)):
            for rec in read_fq(match, skim=skim):
                yield rec

    logger.info(f"iterating over reads from '{fname}'")
    if fname.endswith(".gz"):
        src = FASTQ_src(gzip.open(fname, mode="rt"))
    elif fname.endswith(".bam"):
        src = BAM_src(fname)
    elif type(fname) is str:
        src = FASTQ_src(open(fname))
    else:
        src = FASTQ_src(fname)  # assume its a stream or file-like object already

    n = 0
    for record in timed_loop(src, logger, T=15, template="processed {i} reads in {dT:.1f}sec. ({rate:.3f} k rec/sec)", skim=skim):
        yield record
        n += 1

    logger.info(f"processed {n} FASTQ records from '{fname}'")


def make_header(bam, progname):
    import os
    import sys

    header = bam.header.to_dict()
    # if "PG" in header:
    # for pg in header['PG']:
    #     if pg["ID"] == progname:
    #         progname = progname + ".1"

    pg_list = header.get("PG", [])
    pg = {
        "ID": progname,
        "PN": progname,
        "CL": " ".join(sys.argv),
        "VN": __version__,
    }
    if len(pg_list):
        pg["PP"] = pg_list[-1]["ID"]

    header["PG"] = pg_list + [pg]
    return header


def header_dict_to_text(header):
    buf = []
    for group, data in header.items():
        if type(data) is dict:
            data = [data]

        for d in data:
            valstr = "\t".join([f"{key}:{value}" for key, value in d.items()])
            buf.append(f"@{group}\t{valstr}")

    return "\n".join(buf)



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


# @contextmanager
# def message_aggregation(log_listen="spacemake", print_logger=False, print_success=True):
#     message_buffer = []

#     log = logging.getLogger(log_listen)
#     log.setLevel(logging.INFO)

#     class MessageHandler(logging.NullHandler):
#         def handle(this, record):
#             if record.name == log_listen:
#                 if print_logger:
#                     print(f"{log_listen}: {record.msg}")
#                 else:
#                     print(record.msg)

#     log.addHandler(MessageHandler())

#     try:
#         yield True

#         if print_success:
#             print(f"{LINE_SEPARATOR}SUCCESS!")

#     except SpacemakeError as e:
#         print(e)

def message_aggregation(log_listen="spacemake", print_logger=False, print_success=True):
    from functools import wraps

    def the_decorator(func):
        @wraps(func)
        def wrapper(*argc, **kw):
            message_buffer = []

            log = logging.getLogger(log_listen)
            log.setLevel(logging.INFO)

            class MessageHandler(logging.NullHandler):
                def handle(this, record):
                    if record.name == log_listen:
                        if print_logger:
                            print(f"{log_listen}: {record.msg}")
                        else:
                            print(record.msg)

            log.addHandler(MessageHandler())

            try:
                res = func(*argc, **kw)
            except SpacemakeError as err:
                # print(f"CAUGHT A SPACEMAKE ERROR {err}")
                return err
            else:
                if print_success:
                    print(f"{LINE_SEPARATOR}SUCCESS!")
                return res

        return wrapper        

    return the_decorator

def str_to_list(value):
    # if list in string representation, return the list
    if value is None:
        return None

    if type(value) is str and value.startswith("["):
        if value == "[nan]":
            return []
        else:
            return eval(value)
    # else create a list
    else:
        return [value]


def check_star_index_compatibility(star_index_dir):
    import os

    star_version = os.popen("STAR --version").read().strip()

    with open(os.path.join(star_index_dir, "Log.out"), "r") as f:
        first_line = f.readline().strip()
        index_version = first_line.split("=")[-1].split("_")[-1]

        if index_version != star_version:
            raise SpacemakeError(
                f"STAR index version ({index_version}) is"
                + f" incompatible with your STAR version ({star_version})"
            )


def load_yaml(path, mode="rt"):
    if not path:
        return {}
    else:
        import yaml

        return yaml.load(path, mode)


def setup_logging(
    args,
    name="spacemake.main",
    log_file="",
    FORMAT="%(asctime)-20s\t{sample:30s}\t%(name)-50s\t%(levelname)s\t%(message)s",
):
    sample = getattr(args, "sample", "na")
    import setproctitle
    if name != "spacemake.main":
        setproctitle.setproctitle(f"{name} {sample}")

    FORMAT = FORMAT.format(sample=sample)

    log_level = getattr(args, "log_level", "INFO")
    lvl = getattr(logging, log_level)
    logging.basicConfig(level=lvl, format=FORMAT)
    root = logging.getLogger("spacemake")
    root.setLevel(lvl)

    log_file = getattr(args, "log_file", log_file)
    if log_file:
        fh = logging.FileHandler(filename=ensure_path(log_file), mode="a")
        fh.setFormatter(logging.Formatter(FORMAT))
        root.debug(f"adding log-file handler '{log_file}'")
        root.addHandler(fh)

    if hasattr(args, "debug"):
        # cmdline requested debug output for specific domains (comma-separated)
        for logger_name in args.debug.split(","):
            if logger_name:
                root.info(f"setting domain {logger_name} to DEBUG")
                logging.getLogger(logger_name.replace("root", "")).setLevel(
                    logging.DEBUG
                )

    logger = logging.getLogger(name)
    logger.debug("started logging")
    for k, v in sorted(vars(args).items()):
        logger.debug(f"cmdline arg\t{k}={v}")

    return logger

def setup_smk_logging(name="spacemake.smk", **kw):
    import argparse
    args = argparse.Namespace(**kw)
    return setup_logging(args, name=name)

default_log_level = "INFO"


def make_minimal_parser(prog="", usage="", **kw):
    import argparse

    parser = argparse.ArgumentParser(prog=prog, usage=usage, **kw)
    parser.add_argument(
        "--log-file",
        default=f"{prog}.log",
        help=f"place log entries in this file (default={prog}.log)",
    )
    parser.add_argument(
        "--log-level",
        default=default_log_level,
        help=f"change threshold of python logging facility (default={default_log_level})",
    )
    parser.add_argument(
        "--debug",
        default="",
        help=f"comma-separated list of logging-domains for which you want DEBUG output",
    )
    parser.add_argument(
        "--sample", default="sample_NA", help="sample_id (where applicable)"
    )
    return parser


def load_config_with_fallbacks(args, try_yaml="config.yaml"):
    """
    Tries to load spacemake configuration from
        1) args.config
        2) try_yaml
        3) builtin default from spacamake package

    The entire configuration is attached to the args namespace such that
    args.config["barcode_flavors"]["dropseq] can work.
    """
    from spacemake.config import ConfigFile

    config = None
    if hasattr(args, "config") and os.access(args.config, os.R_OK):
        config = ConfigFile.from_yaml(args.config)
    elif os.access(try_yaml, os.R_OK):
        config = ConfigFile.from_yaml(try_yaml)
    else:
        builtin = os.path.join(os.path.dirname(__file__), "data/config/config.yaml")
        config = ConfigFile.from_yaml(builtin)

    args_kw = vars(args)
    # try:
    #     args_kw["log_level"] = config.variables["logging"]["level"]
    # except KeyError:
    #     pass

    args_kw["config"] = config.variables
    import argparse

    return argparse.Namespace(**args_kw)
