import logging
import datetime
import argparse
import numpy as np

counted_regions = ["UTR", "CODING"]


def get_tag_on_split(row, tag="XF:Z:", default="INTERGENIC"):
    for col in row[11:]:
        if col.startswith(tag):
            return col.split(tag)[1]

    return default


def select_alignment(alignments):
    pieces = [aln.split("\t") for aln in alignments]
    read_names = [p[0] for p in pieces]
    if read_names.count(read_names[0]) != len(read_names):
        # print(read_names)
        raise Exception(f"input alignments do not come from the same read")

    def is_exonic(row):
        xf = get_tag_on_split(row)
        return xf in counted_regions

    alignments_are_exonic = np.array([is_exonic(row) for row in pieces])
    exonic_ix = np.where(alignments_are_exonic == True)[0]

    num_exonic = exonic_ix.shape[0]

    if num_exonic == 1:
        # if only one exonic reads from the group
        # return the exonic indices
        return alignments[exonic_ix[0]]
    else:
        return None


def load_barcodes(fname):
    logger = logging.getLogger("spacemake.scripts.filter_mm_reads.load_barcodes")
    if fname.endswith(".gz"):
        import isal.igzip

        _open = isal.igzip.open
    else:
        _open = open

    bcs = set()
    for line in _open(fname, "rt"):
        if line.startswith("cell"):
            continue

        bcs.add(line.split("\t")[0])

    logger.info(f"loaded {len(bcs)} barcodes from {fname}")
    return bcs


def filter_mm(input, _out, bcs=set(), **kw):
    logger = logging.getLogger("spacemake.scripts.filter_mm_reads")
    from collections import defaultdict

    counter = defaultdict(int)

    import re

    CB_pattern = re.compile(f"CB:Z:(\S+)")

    from time import time

    T0 = time()

    output = open(_out, "wt")
    multi_mappers = []
    qname = None
    for aln in open(input, "rt"):
        # header line. Just pass through
        if aln.startswith("@"):
            output.write(aln)

        # counting and rate info
        counter["N_alignments"] += 1
        if counter["N_alignments"] % 10000000 == 0:
            dT = time() - T0
            rate = counter["N_alignments"] / 1000.0 / dT
            logger.info(
                f"processed {counter['N_alignments']} alignments in {dT:.1f} seconds ({rate:.1f}k/sec)"
            )

        # restrict output to only desired barcode subset
        if bcs:
            m_CB = re.search(CB_pattern, aln)
            if not m_CB:
                # This should never happen at this stage
                counter["N_no_CB"] += 1
                continue

            CB = m_CB.groups(0)[0]
            if CB not in bcs:
                counter["N_CB_not_selected"] += 1
                continue
            else:
                counter["N_CB_selected"] += 1

        query_name = aln.split("\t", maxsplit=1)[0]
        if query_name != qname:
            # new read
            if len(multi_mappers) == 1:
                # fast path
                counter["N_unique"] += 1
                output.write(aln)

            elif len(multi_mappers) > 1:
                counter["N_multi"] += 1
                # decide which, if any, to keep
                aln_to_keep = select_alignment(multi_mappers)
                if aln_to_keep is not None:
                    counter["N_salvaged"] += 1
                    # set aln secondary flag to 0, so that it is flagged as primary
                    # secondary flag is at 0x100, so 8th bit (starting from 0)
                    cols = aln_to_keep.split("\t")
                    flag = int(cols[1])
                    flag = flag & ~(1 << 8)
                    cols[1] = str(flag)
                    output.write("\t".join(cols))
                else:
                    counter["N_not_salvaged"] += 1

            # reset multimapper list
            multi_mappers = []
            qname = query_name

        # add the last alignment
        multi_mappers.append(aln)

    # final iteration:
    if len(multi_mappers) == 1:
        counter["N_unique"] += 1
        output.write(aln)

    elif len(multi_mappers) > 1:
        counter["N_multi"] += 1
        # decide which, if any, to keep
        aln_to_keep = select_alignment(multi_mappers)
        if aln_to_keep is not None:
            counter["N_salvaged"] += 1
            # set aln secondary flag to 0, so that it is flagged as primary
            # secondary flag is at 0x100, so 8th bit (starting from 0)
            cols = aln_to_keep.split("\t")
            flag = int(cols[1])
            flag = flag & ~(1 << 8)
            cols[1] = str(flag)
            output.write("\t".join(cols))
        else:
            counter["N_not_salvaged"] += 1

    output.flush()
    return counter


if __name__ == "__main__":
    import spacemake.util as util

    parser = util.make_minimal_parser(
        description="Filter out ambiguous multi-mapper reads"
    )

    parser.add_argument("--in-sam", help="input sam", default="/dev/stdin")
    parser.add_argument(
        "--barcode-list",
        help="[optional] only pass reads with CB:Z:<barcode> from this (compressed) table's first column",
    )
    parser.add_argument("--out-sam", help="output sam", default="/dev/stdout")
    # parser.add_argument("--sample", help="sample_id", default="NA")

    args = parser.parse_args()
    logger = util.setup_logging(args, name="spacemake.scripts.filter_mm_reads")
    logger.info("starting up")
    if args.barcode_list:
        bcs = load_barcodes(args.barcode_list)
    else:
        bcs = set()

    start_time = datetime.datetime.now()
    counter = filter_mm(args.in_sam, args.out_sam, bcs=bcs)
    finish_time = datetime.datetime.now()

    formatted_time = finish_time.strftime("%Y-%m-%d %H:%M:%S")
    total_elapsed_seconds = (finish_time - start_time).total_seconds()

    logger.info(
        f"Finished processing {counter['N_alignments']:,} records in {total_elapsed_seconds:,.0f} seconds. Current time: {formatted_time}"
    )
