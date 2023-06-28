from spacemake.util import read_fq
import os
import numpy as np
from collections import defaultdict
import logging
import matplotlib
import spacemake.util as util
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


class NTCounter:
    def __init__(self, name="name"):
        self.name = name
        self.n = 0
        self.L = 0
        self.min_x = 0
        self.max_x = 0
        self.counter = defaultdict(lambda: defaultdict(int))

    def count_seq(self, seq, ofs=0):
        self.n += 1
        self.L = max(self.L, len(seq))
        self.min_x = min(self.min_x, -ofs)
        self.max_x = max(self.max_x, len(seq) - ofs)

        for i, nt in enumerate(seq):
            self.counter[nt][i - ofs] += 1

    def get_percentages(self, nt):
        pos = np.arange(self.min_x, self.max_x)
        freq = np.array([self.counter[nt][i] for i in pos])
        if self.n <= 0:
            return np.ones(self.L) * np.nan
        else:
            return 100.0 * freq / float(self.n)

    def get_df(self, bases="ACGTN", **kw):
        data = dict(kw)
        data['pos']= np.arange(self.min_x, self.max_x)
        for nt in bases:
            data[nt] = self.get_percentages(nt)

        return pd.DataFrame(data).set_index('pos')

    def plot_percentages(self, ax, legend=True, bases = "ACGTN"):
        ax.set_axisbelow(True)
        ax.grid(True, axis="y")
        ax.set_title(f"{self.name}")

        pos = np.arange(self.min_x, self.max_x)
        bottom = np.zeros_like(pos, dtype=float)

        colors = {
            "A": "#f9c734",
            "C": "#F9575A",
            "G": "#6E74F0",
            "T": "#37A967",
            "N": "#999999",
        }

        for nt in bases:
            # print(pos)
            freq = self.get_percentages(nt)
            # print(freq)
            ax.bar(
                pos,
                height=freq,
                bottom=bottom,
                label=nt,
                width=0.8,
                align="center",
                color=colors[nt],
            )
            bottom += freq

        ax.set_xlim(-0.6, self.L - 0.6)
        ax.set_xticks(pos[0::4])
        ax.set_ylim(0, 100)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlabel("position [nt]")
        ax.set_ylabel("frequency [%]")
        ax.set_yticks([0, 25, 50, 75, 100])
        # ax.set_xticks([1, 4, 8, 12, 16, 20])
        if legend:
            ax.legend(ncol=1, frameon=False, loc="center left", bbox_to_anchor=(1, 0.5))

        ax.set_xlim(self.min_x - 1, self.max_x + 1)

def parse_args():
    parser = util.make_minimal_parser()
    parser.add_argument("fname")
    parser.add_argument(
        "--skim",
        default=10,
        type=int,
        help="skim through the file by counting only every n-th read. set to 0 to disable (default=10)",
    )
    parser.add_argument("--out-pdf", default="nt_freq_{sample}.pdf")
    parser.add_argument("--out-csv", default="nt_freq_{sample}.csv")
    parser.add_argument("--anchor", default="", help="OPTIONAL: specify an anchoring regexp. Only matching reads will be scanned and the 3' end of the match which will define position 0")
    return parser.parse_args()

def filter_seq(seq, regexp):
    if not regexp:
        return 0

    import re
    M = re.search(regexp, seq)
    if M:
        return M.end()
    else:
        return False

def main(args):
    sample = os.path.basename(args.fname)
    logger = logging.getLogger("spacemake.bin.plot_base_freqs")

    if sample.endswith("fastq.gz"):
        logger.debug("reading from FASTQ")
        counter = NTCounter(sample)

        for name, seq, qual in read_fq(args.fname, skim=args.skim):
            x = filter_seq(seq, args.anchor)
            if x or (not args.anchor):
                counter.count_seq(seq, ofs=x)

        if args.out_csv:
            df = counter.get_df(sample=sample, anchor=args.anchor)
            df.to_csv(args.out_csv, sep='\t')

        fig, ax = plt.subplots(figsize=(10, 5))
        counter.plot_percentages(ax)
        fig.tight_layout()
        fig.savefig(args.out_pdf.format(sample=sample))

    elif sample.endswith("bam"):
        logger.debug("reading from BAM")
        import pysam

        bam = util.quiet_bam_open(args.fname, check_sq=False, threads=4)
        R2_counter = NTCounter("read2")
        CB_counter = NTCounter("CB")
        UMI_counter = NTCounter("UMI")

        n = 0
        for read in bam.fetch(until_eof=True):
            n += 1
            if args.skim and n % args.skim != 0:
                continue

            x = filter_seq(read.query_sequence, args.anchor)
            if x or (not args.anchor):
                R2_counter.count_seq(read.query_sequence, ofs=x)
                CB_counter.count_seq(read.get_tag("CB"))
                UMI_counter.count_seq(read.get_tag("MI"))

        if args.out_csv:
            df = pd.concat([
                R2_counter.get_df(name="read2"),
                CB_counter.get_df(name="BC"),
                UMI_counter.get_df(name="UMI"),
            ])
            df['sample'] = sample
            df['anchor'] = args.anchor
            df.to_csv(args.out_csv, sep='\t')
        
        fig, (ax_r2, ax_cb, ax_umi) = plt.subplots(
            nrows=1,
            ncols=3,
            figsize=(12, 5),
            gridspec_kw=dict(
                width_ratios=[R2_counter.L, CB_counter.L, UMI_counter.L], wspace=0.1
            ),
            sharey=True,
        )

        fig.suptitle(sample)
        R2_counter.plot_percentages(ax_r2, legend=False)
        CB_counter.plot_percentages(ax_cb, legend=False)
        UMI_counter.plot_percentages(ax_umi)
        
        pdfout = args.out_pdf.format(sample=sample)
        logger.info(f"storing barcode & UMI nt-freq profile in '{pdfout}'")
        fig.tight_layout()
        fig.savefig(pdfout)
    else:
        raise ValueError(f"unknown file format {args.fname}")

def cmdline():
    args = parse_args()
    util.setup_logging(args, "spacemake.bin.plot_base_freqs")
    main(args)

if __name__ == "__main__":
    cmdline()