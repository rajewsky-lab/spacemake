import mrfifo as mf
import logging


def parse_args():
    from spacemake.util import make_minimal_parser

    parser = make_minimal_parser("BamTagHistogram")

    parser.add_argument("--parallel", type=int, default=8)
    parser.add_argument("--input", default="/dev/stdin")
    parser.add_argument("--output", default="/dev/stdout")
    parser.add_argument(
        "--prefix-size",
        default=4,
        type=int,
        help=(
            "how many letters of the tag value are used to split the stream. "
            "default=4 allows for up to (alphabet_size)^4 distinct parallel workers. "
            "will be spread across workers by mod <args.parallel>"
        ),
    )
    parser.add_argument("--prefix-alphabet", default="ACGTN")
    parser.add_argument("--min-count", default=10, type=int)
    parser.add_argument(
        "--sort-mem",
        default=8,
        type=int,
        help="how many GB are allowed to be used for sorting (default=8)",
    )
    parser.add_argument(
        "--tag", default="CB", help="which BAM tag to count (default='CB')"
    )

    return parser.parse_args()


def CB_distributor(
    input, outputs, tag="CB", prefix_size=3, prefix_alphabet="ACGTN", n=8, **kw
):
    "ensure that the FIFOs are not managed"
    assert type(input) is str
    logger = logging.getLogger("mrfifo.parts.CB_distributor")
    logger.info(
        f"reading from {input}, writing to {outputs} "
        f"tag={tag} prefix_size={prefix_size} prefix_alphabet={prefix_alphabet} "
        f"kw={kw}"
    )

    lkup = {}
    from itertools import product

    i = 0
    for letters in product(*([prefix_alphabet] * prefix_size)):
        prefix = "".join(letters).encode("ascii")
        lkup[prefix] = i % n
        i += 1

    # for k, v in sorted(lkup.items()):
    #     print(f"{k}\t{v}")

    from mrfifo.fast_loops import distribute_by_substr

    tag_lead = b"\t" + tag.encode("ascii") + b":Z:"
    logger.debug(
        f"scanning for tag-lead {tag_lead} and using next {prefix_size} bytes as prefix"
    )
    res = distribute_by_substr(
        fin_name=input,
        fifo_names=outputs,
        sub_lookup=lkup,
        sub_size=prefix_size,
        sub_lead=tag_lead,
        # **kw,
    )
    logger.info("distribution complete")
    return res


def tag_counter(input, output, tag="CB", min_count=10):
    from collections import defaultdict

    counter = defaultdict(int)
    import re

    n_lines = 0
    n_tags = 0
    pattern = re.compile(f"{tag}:Z:(\S+)")
    for sam_line in input:
        n_lines += 1
        if m := re.search(pattern, sam_line):
            n_tags += 1
            tag_val = m.groups(0)[0]
            counter[tag_val] += 1

    n_above_cut = 0
    for value, count in counter.items():
        if count >= min_count:
            n_above_cut += 1
            output.write(f"{value}\t{count}\n")

    return n_lines, n_tags, len(counter), n_above_cut


def sort_function(input, output, n=8, sort_mem_gigs=8):
    import os

    if output.endswith(".gz"):
        os.system(
            f"sort -rnk 2 -S {sort_mem_gigs}G --parallel={n} {input} "
            f"| python -m isal.igzip -c > {output}"
        )
    else:
        os.system(f"sort -rnk 2 -S {sort_mem_gigs}G --parallel={n} {input} > {output}")


def main(args):
    w = (
        mf.Workflow("BamTagHistogram")
        .BAM_reader(
            input=args.input,
            mode="S",
        )
        .distribute(
            input=mf.FIFO("input_sam", "rt"),
            outputs=mf.FIFO("dist_{n}", "wt", n=args.parallel),
            func=CB_distributor,
            tag=args.tag,
            prefix_size=args.prefix_size,
            prefix_alphabet=args.prefix_alphabet,
            n=args.parallel,
        )
        .workers(
            func=tag_counter,
            tag=args.tag,
            input=mf.FIFO("dist_{n}", "rt"),
            output=mf.FIFO("counts_{n}", "wt"),
            n=args.parallel,
            min_count=args.min_count,
        )
        .collect(
            inputs=mf.FIFO("counts_{n}", "rt", n=args.parallel),
            output=mf.FIFO("unsorted", "wt"),
            chunk_size=1,
        )
        .funnel(
            input=mf.FIFO("unsorted", "rt"),
            output=args.output,
            func=sort_function,
            _manage_fifos=False,
        )
        .run()
    )
    print(w.result_dict)
    return w


if __name__ == "__main__":
    args = parse_args()
    import spacemake.util as util

    util.setup_logging(args)
    main(args)
