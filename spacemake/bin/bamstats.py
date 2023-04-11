import spacemake.util as util
from spacemake.parallel import CountingStatistics, AlignedSegmentsFromQueue, ExceptionLogging

def parse_args():
    parser = util.make_minimal_parser(prog="bamstats.py")
    parser.add_argument("bamfile", help="input bam file")
    parser.add_argument("--umi-cutoff", type=int, default=0, help="keep CB only if UMI >= cutoff (default=0)")
    parser.add_argument("--out-stats", default="{bampath}/{bamname}.bc_stats.tsv", help="barcode and UMI statistics output")
    parser.add_argument("--out-length", default="{bampath}/{bamname}.readlen.tsv", help="read-length distribution output file")
    args = parser.parse_args()

    return args


def scan(Qin, Qout, Qerr, abort_flag, stat_list, args={}, count_class=CountingStatistics, skim=0, **kw):
    from collections import defaultdict

    with ExceptionLogging(
        "spacemake.bamstats.scan", Qerr=Qerr, exc_flag=abort_flag
    ) as el:

        UMI_counter = defaultdict(set)
        read_counter = defaultdict(int)
        read_len_hist = defaultdict(int)

        # print(fname)
        bam_src = AlignedSegmentsFromQueue(Qin, abort_flag)

        # sam = util.quiet_bam_open(fname, "rb", check_sq=False, threads=2)
        last_qname = ""
        for ref, read in util.timed_loop(bam_src, el.logger, skim=skim):
            # count each read only once regardless of how many times it aligns
            if read.query_name == last_qname:
                continue

            last_qname = read.query_name

            CB = read.get_tag("CB")
            MI = read.get_tag("MI")

            UMI_counter[CB].add(hash(MI))
            read_counter[CB] += 1

            seq = read.query_sequence
            read_len_hist[len(seq)] += 1

        # print("done processing")
        import numpy as np

        stats = count_class()
        for cb, umis in UMI_counter.items():
            stats.stats_by_ref['UMI_counts'][cb] = len(umis)

        stats.stats_by_ref['read_counts'] = read_counter
        stats.stats_by_ref['readlen'] = read_len_hist

        stat_list.append(stats.stats_by_ref)


def main(args):
    import os

    logger = util.setup_logging(args, "spacemake.bamstats")
    bamname = os.path.basename(args.bamfile).split('.bam')[0]
    bampath = os.path.dirname(args.bamfile)
    sample = args.sample.format(bamname=bamname)

    from spacemake.parallel import parallel_BAM_workflow
    f_stats_name = util.ensure_path(args.out_stats.format(bamname=bamname, bampath=bampath))
    f_len_name = util.ensure_path(args.out_length.format(bamname=bamname, bampath=bampath))

    logger = util.setup_logging(args, f"spacemake.bamstats.main_parallel")
    logger.info("startup")

    stats = parallel_BAM_workflow(
        [args.bamfile],
        scan,
        None,
        n_workers=4,
        buffer_size=100000,
        log_domain="spacemake.bamstats",
        args=args
    )

    if not stats:
        return -1

    # format and save output
    import numpy as np
    import pandas as pd

    CBs = []
    UMI_counts = []
    read_counts = []
    d = stats.stats_by_ref
    for cb, n_umis in sorted(d['UMI_counts'].items(), key=lambda x: -x[1]):
        if n_umis >= args.umi_cutoff:
            CBs.append(cb)
            UMI_counts.append(n_umis)
            read_counts.append(d['read_counts'][cb])
            # print(f"{cb}\t{n_umis}")

    if len(d['readlen']):
        Lmax = np.array(list(d['readlen'].keys())).max()
    else:
        Lmax = 0

    Lhist = np.array([d['readlen'][x] for x in np.arange(Lmax + 1)])

    df = pd.DataFrame(dict(sample=sample, bam=bamname, CB=CBs, counts=UMI_counts, reads=read_counts))
    df.to_csv(f_stats_name, sep="\t", index=None)

    df = pd.DataFrame(dict(sample=sample, bam=bamname, readlen=np.arange(len(Lhist)), count=Lhist))
    df.to_csv(f_len_name, sep="\t", index=None)


if __name__ == "__main__":
    args = parse_args()
    main(args)
