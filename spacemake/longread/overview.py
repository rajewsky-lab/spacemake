import os
import pandas as pd
import numpy as np
from glob import glob
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


def main(args):
    dfs = []
    for fname in list(args.fnames) + list(glob(args.glob_pattern)):
        print(f"loading {fname}")
        df = pd.read_csv(fname, sep='\t')
        df['stats_file'] = fname
        dfs.append(df)

    df = pd.concat(dfs)

    repriming = ['TSO,TSO_RC', 'dN-SMRT,dN-SMRT_RC', ]
    concatenation = [c for c in df.columns if c.endswith('+') and ',' not in c]
    bead = ["bead_complete", "bead_only_handle", "bead_no_dT", "bead_no_opseq"][::-1]

    # avoid crash if columns are missing
    for r in repriming + concatenation + bead:
        if r not in df.columns:
            df[r] = 0

    # print(df)
    # print(f"concat columns {concatenation}")
    # print(f"bead columns {bead}")
    df['reprimed'] = df[repriming].sum(axis=1)
    df['bead_complete'] = np.nan_to_num(df['bead_complete'], nan=0.0)
    df['concat'] = df[concatenation].sum(axis=1)
    df['bead_related'] = np.nan_to_num(df[bead].sum(axis=1), nan=0.0)
    df['bead_dropseq'] = np.nan_to_num(df['bead_no_opseq'], nan=0.0)
    df['bead_incomplete'] = df['bead_related'] - df['bead_complete'] - df['bead_dropseq']
    df['non_bead'] = 100 - df['bead_related']
    df['bead_fidelity'] = 100 * df['bead_complete'] / df['bead_related']
    df = df.fillna(0)
    # print(df)
    if args.csv_out:
        df.to_csv(args.csv_out, float_format='%.2f', sep='\t', index=False)

    def clean(txt):
        txt = os.path.basename(txt)
        t = txt\
            .replace('source/','') \
            .replace('sts_', '') \
            .replace('pb_', '') \
            .replace('ds_', '') \
            .replace('.fq', '') \
            .replace('.bam', '') \
            .replace('lima.', '')
        
        if t.count('_') > 1:
            t = "_".join(t.split('_')[:2])
        
        return t

    df['name'] = df['qfa'].apply(clean)
    # df = df.sort_values('bead_related')
    df = df.sort_values('name')

    def guess_rRNA_file(path):
        # print("guessrRNA raw path", path)
        name = os.path.basename(path).replace('.summary', '.rRNA')
        
        if args.rRNA_same_place:
            place = os.path.dirname(path)
        else:
            place = args.rRNA

        return [
            os.path.join(place, name.replace(".fq", ".txt")),
            os.path.join(place, name.replace(".fq", ".txt")).replace('.rRNA.tsv', '.txt'),
            os.path.join(place, name.replace(".fq", ".txt")).replace('.rRNA.tsv', '.rRNA.txt'),
            os.path.join(place, name.replace(".bam", ".txt").replace("lima.", "")),
            os.path.join(place, name.replace(".bam", ".txt").replace("lima.", "")).replace('.rRNA.tsv', '.txt'),
            os.path.join(place, name.replace(".bam", ".txt").replace("lima.", "")).replace('.rRNA.tsv', '.rRNA.txt'),
        ]

    rRNA_fracs = []
    for row in df[['stats_file', 'N_reads']].itertuples():
        rcount = np.nan
        for fname in guess_rRNA_file(row.stats_file):
            print(fname)
            try:
                rcount = int(open(fname).read())
            except (FileNotFoundError, ValueError):
                pass
            else:
                break
        if rcount == np.nan:
            raise ValueError

        rRNA_fracs.append(100. * rcount / row.N_reads)

    df['rRNA'] = rRNA_fracs
    # print(df[['qfa', 'rRNA']])

    def make_bars(ax, df, kinds, labels, cmap=plt.get_cmap('tab10'), w=0.9, colors=None):
        n = len(kinds)
        if colors is None:
            colors = cmap(np.linspace(0, 1, n))

        x = np.arange(len(df)) - w/2.0
        y0 = np.zeros(len(x), dtype=float)
        for kind, label, color in zip(kinds, labels, colors):
            y = np.nan_to_num(df[kind], nan=0.0)
            # print(kind)
            # print(y)
            ax.bar(x, y, bottom=y0, label=label, width=w, color=color)
            y0 += y

        ax.set_ylabel('fraction of library')
        ax.set_xticks(x)
        labels = df['name']  # [clean(fq) for fq in df['qfa']]
        ax.set_xticklabels(labels, rotation=90)
        ax.set_ylim(0, 100)

    marie = ["non_bead", "bead_incomplete", "bead_dropseq", "bead_complete", ]
    marie_colors = ["gray", "royalblue", "green", "gold"]
    
    w = max(8 / 25. * len(df), 3)
    if args.multi_page:
        pdf = PdfPages(args.breakdown)
        fig, ax1 = plt.subplots(1, figsize=(w, 4))
    else:
        fig, (ax1, ax2) = plt.subplots(2, figsize=(w, 6), sharex=True)

    make_bars(ax1, df, marie, labels=[b.replace('bead_', '') for b in marie], colors=marie_colors)
    ax1.legend(title='Marie-stats', ncol=len(marie))
    if args.multi_page:
        fig.tight_layout()
        pdf.savefig()
        plt.close()
        fig, ax2 = plt.subplots(1, figsize=(w, 4))

    make_bars(ax2, df, ["bead_fidelity"], labels=["bead fidelity"])
    ax2.set_ylabel("bead fidelity")
    if args.multi_page:
        fig.tight_layout()
        pdf.savefig()
        pdf.close()
    else:
        fig.tight_layout()
        plt.savefig(args.breakdown)

    plt.close()

    if args.multi_page:
        pdf = PdfPages(args.output)
        fig, ax1 = plt.subplots(1, figsize=(w, 4))
    else:
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize=(w, 12), sharex=True)

    # print("bead related", bead)
    make_bars(ax1, df, bead, labels=[b.replace('bead_', '') for b in bead])
    ax1.legend(title='bead-related', ncol=len(bead))
    if args.multi_page:
        fig.tight_layout()
        pdf.savefig()
        plt.close()
        fig, ax2 = plt.subplots(1, figsize=(w, 4))

    # print("repriming events", repriming)
    make_bars(ax2, df, repriming, labels=[r.split(',')[0] for r in repriming], cmap=plt.get_cmap('tab20c'))
    ax2.legend(title='repriming', ncol=len(repriming))
    if args.multi_page:
        fig.tight_layout()
        pdf.savefig()
        plt.close()
        fig, ax3 = plt.subplots(1, figsize=(w, 4))

    # print("concat events", concatenation)
    make_bars(ax3, df, concatenation, labels=concatenation, cmap=plt.get_cmap('tab20b'))
    ax3.legend(title='concatamers', ncol=len(concatenation))
    if args.multi_page:
        fig.tight_layout()
        pdf.savefig()
        plt.close()
        fig, ax4 = plt.subplots(1, figsize=(w, 4))

    make_bars(ax4, df, ["rRNA",], labels = ["rRNA"], cmap=plt.get_cmap('tab20c'))
    ax4.legend(title='human rRNA', ncol=1)
    if args.multi_page:
        fig.tight_layout()
        pdf.savefig()
        pdf.close()
    else:
        fig.tight_layout()
        plt.savefig(args.output)

    plt.close()


def setup_parser(parser):
    parser.add_argument("fnames", nargs='*')
    parser.add_argument("--output", default="pb_overview.pdf",
                        help="path/name of detailed report PDF")
    parser.add_argument("--csv-out", default="all_pb_stats.csv",
                        help="path/name of detailed report PDF")
    parser.add_argument("--breakdown", default="bead_overview.pdf",
                        help="path/name of bead report (Marie style) PDF")
    parser.add_argument("--glob-pattern", default="stats/*summary.tsv",
                        help="search pattern to gather summary files generated by the scan command")
    parser.add_argument("--rRNA", default="rRNA/", 
                        help="path to search for rRNA counts corresponding to samples")
    parser.add_argument("--rRNA-same-place", default=False, action='store_true',
                        help="If set, look for rRNA txt file with same sample name in same directory")
    parser.add_argument("--multi-page", default=False, action="store_true",
                        help="If set, generate multiple PDF pages instead of subplots")


if __name__ == "__main__":
    # setup own parser
    import argparse
    parser = argparse.ArgumentParser(prog='pb_overview')
    setup_parser(parser)
    main(parser.parse_args())
