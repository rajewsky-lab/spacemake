import numpy as np
import pandas as pd
import matplotlib as mpl
import spacemake.util as util

mpl.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.tab20.colors)


def reflines(ax, N_total_x=1, N_total_y=1):
    ax.plot([0, N_total_x], [0, N_total_y], "-k", lw=0.5, label="100%")
    ax.plot([0, N_total_x], [0, N_total_y/10], "--k", lw=0.5, label="10%")
    ax.plot([0, N_total_x], [0, N_total_y/100], "-.k", lw=0.5, label="1%")
    ax.plot([0, N_total_x], [0, N_total_y/1000], ":k", lw=0.5, label="0.1%")


def parse_args():
    parser = util.make_minimal_parser(
        prog="plot_complexity.py",
        description="plot library complexity from sub-sampling analysis"
    )

    parser.add_argument(
        "subsample_results",
        help="list of subsampling result files to plot",
        default="",
        nargs="+",
    )

    parser.add_argument(
        "--out-plot",
        help="name of PDF (default=complexity.pdf)'",
        default="complexity.pdf",
    )

    parser.add_argument(
        "--cutadapt",
        help="path to adapter-trimming stats",
        default="",
    )

    args = parser.parse_args()
    return args

def get_N_kept(fname):
    N = np.NaN
    if fname:
        df = pd.read_csv(fname, sep='\t')
        N = df.loc['reads'].query('key == "N_kept"')['count'].iloc[0]
    
    return float(N)

def main(args):
    logger = util.setup_logging(args, "spacemake.bin.plot_complexity")
    fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(10, 10))
    fig.suptitle(args.sample)
    N_kept = get_N_kept(args.cutadapt)

    maxy = 0
    for fname in args.subsample_results:
        df = pd.read_csv(fname, sep="\t", comment="#", index_col=False)
        name = df['bamname'].iloc[0]
        N_reads = df["n_reads"].values.max()

        frac = N_reads/N_kept
        logger.debug(f"N_kept={N_kept} {name} frac={frac:.3e}")
        df["kept_est"] = df["n_reads"] / frac
        df["percent"] = 100.0 * df["n_reads"] / N_reads

        label = f"{name} ({N_reads/1e6:.1f} M)"
        
        ax1.plot(df["kept_est"]/1e6, df["n_umis"]/1e6, ".-", label=label)
        ax2.plot(df["percent"], df["n_umis"]/1e6, ".-")
        ax3.plot(df["percent"], df["n_reads"] / df["n_umis"], ".-")

        maxy = max(maxy, df["n_umis"].max()/1e6)


    reflines(ax1, N_kept/1e6, N_kept/1e6)
    reflines(ax2, 100.0, N_kept/1e6)

    ax1.set_ylim(0, maxy * 1.1)
    ax2.set_ylim(0, maxy * 1.1)
    ax1.legend(bbox_to_anchor=(0.5, 1.02), loc="lower center", ncol=3)
    ax1.set_ylabel("UMIs [M]")
    ax2.set_ylabel("UMIs [M]")
    ax2.set_xlabel("subsample %")
    ax1.set_xlabel("trimmed reads (est.) [M]")
    ax3.set_xlabel("subsample %")
    ax3.set_ylabel("reads-to-UMI ratio")
    fig.tight_layout()
    logger.info(f"saving '{args.out_plot}' ...")
    plt.savefig(args.out_plot)

if __name__ == "__main__":
    args = parse_args()
    main(args)