import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np

def plot_cell_hist(df, key="n_reads", bins=100, xlog=True, ylog=False, ax=None, **kw):
    if ax is None:
        _, ax = plt.subplots(figsize=(4, 3))

    values = df[key].values
    med = np.median(values)
    # print(values)
    if xlog:
        lv = np.log10(values + 1)
    else:
        lv = values

    # exclude extreme outliers
    if len(lv) > 1000:
        lmin, lmax = np.percentile(lv, [0.1, 99.9])
        lv = lv[(lv >= lmin) & (lv <= lmax)]

    hist, bin_edges = np.histogram(lv, bins=bins)
    if xlog:
        bin_edges = 10 ** bin_edges

    w = bin_edges[1:] - bin_edges[:-1]
    ax.bar(bin_edges[:-1], hist, width=w, **kw)

    ax.axvline(med, color='r', label=f'median={med:.1f}')
    if xlog:
        ax.set_xscale("log")

    if ylog:
        ax.set_yscale("log")

    ax.set_xlabel(f"{key}")
    ax.set_ylabel("no. cell barcodes")

    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)

    # Only show ticks on the left and bottom spines
    # ax.yaxis.set_ticks_position("left")
    # ax.xaxis.set_ticks_position("bottom")
    ax.legend(frameon=False)
    for axis in [ax.xaxis, ax.yaxis]:
        formatter = ScalarFormatter()
        formatter.set_scientific(False)
        axis.set_major_formatter(formatter)

        # formatter = ScalarFormatter()
        # formatter.set_scientific(False)
        # axis.set_minor_formatter(formatter)


def plot_loglog_knee(adata, key="n_counts", ax=None):
    df = adata.obs
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))

    # UMI = np.sort(np.array(df[key].values, dtype=int))[::-1]
    UMI = np.array(df[key].values, dtype=int)
    if not len(UMI):
        return

    UMIsum = UMI.sum()
    # ax.plot(UMI.values.cumsum()[:n_rank], label=f"{key} (total={UMIsum/1e6:.1f} M)")
    # upper_q = np.percentile(UMI, 95)
    # if idx == "miRNA":
    # axh.hist(UMI, bins=100, label=idx, histtype="step", range=(0, upper_q))
    ax.loglog(
        len(UMI) - np.bincount(UMI).cumsum(),
        label=f"{key} (total={UMIsum/1e6:.1f} M)",
        alpha=0.75,
        lw=3,
    )

    # ax.axvline(100, color="k", ls="dashed", lw=0.1)
    ax.set_ylabel("cumulative no. cells")
    ax.set_xlabel(f"{key} > x")

    # ax.set_ylabel(f"cumulative {key}")
    # ax.set_xlabel("cell BC rank")
    # ax.set_xlim(0, n_rank)
    ax.legend()
    # fig.tight_layout()
    ax.grid(axis="y")
    if fig:
        sample_name = adata.uns.get("sample_name", "no-sample-name")
        fig.suptitle(sample_name)
        fig.tight_layout()

    return fig


def plot_dge_stats(adata, hist_kw={}, fig_kw={'figsize': (8,5)}):
    import scanpy as sc
    sc.pp.filter_cells(adata, min_counts=1)
    # re-generating the metrics
    adata.obs["n_exonic_reads"] = adata.layers["exonic_reads"].sum(axis=1)[:, 0]
    adata.obs["n_reads"] = adata.layers["reads"].sum(axis=1)[:, 0]
    adata.obs["n_exonic_counts"] = adata.X.sum(axis=1)[:, 0]
    adata.obs["n_genes"] = (adata.X > 0).sum(axis=1)[:, 0]
    adata.obs["reads_per_counts"] = (
        adata.obs["n_exonic_reads"] / adata.obs["n_exonic_counts"]
    )

    # print(adata)
    # print(adata.obs)
    # print(adata.var)

    # cell metrics plot
    fig, axes = plt.subplots(2, 2, **fig_kw)
    axes = np.array(axes).ravel()
    plot_cell_hist(adata.obs, key="n_exonic_reads", ax=axes[0], **hist_kw)
    plot_cell_hist(adata.obs, key="n_counts", ax=axes[1], **hist_kw)
    plot_cell_hist(adata.obs, key="n_genes", ax=axes[2], **hist_kw)
    plot_cell_hist(adata.obs, key="reads_per_counts", xlog=False, ax=axes[3], **hist_kw)
    fig.subplots_adjust(hspace=0.4)
    # fig.tight_layout()

    return fig



def scatter_plot(df, x, y, hilight=[], cutoff=10, xlabel=None, ylabel=None, cmap='tab10', title='counts', ax=None, default_figsize=(5,5)):
    """
    Produce a square, log-scaled scatter plot of the values in df, plotting column 'y' against column 'x'.
    If column names '{x}_lo' and '{x}_hi' are detected, they define errorbars for the plot.
    axes can be labelled and specific genes (df rows found by index) can be highlighted.
    """
    import matplotlib.pyplot as plt
    # restrict to genes that are properly expressed
    df = df.loc[(df > cutoff).any(axis=1)]

    if ax is None:
        fig, ax = plt.subplots(figsize=default_figsize)

    ax.set_prop_cycle(color=plt.get_cmap(cmap).colors)
    ax.set_aspect(1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    xerr = []
    yerr = []
    if x+'_lo' in df.columns:
        xerr = np.array([df[x] - df[x+'_lo'], df[x+'_hi'] - df[x]])

    if y+'_lo' in df.columns:
        yerr = np.array([df[y] - df[y+'_lo'], df[y+'_hi'] - df[y]])

    from scipy.stats import pearsonr
    R, pval = pearsonr(np.log10(df[x] + 1), np.log10(df[y] + 1))
    print(f"R={R:.2f} pval={pval:.2e}")
    #if xerr or yerr:
    #    ax.errorbar(x, y, data=df+1, xerr=xerr, yerr=yerr, color="gray", label=f'all (R={R:.2f})', alpha=0.4, fmt='.')
    #else:
    ax.plot(x, y, '.', data=df+1, label=f'all (R={R:.2f})', color="gray",
            ms=5, mew=0,
            alpha=0.4,)

    if y+'_lo' in df.columns and x+'_lo' in df.columns:
        up = df[y+'_lo'] > df[x+'_hi']
        do = df[y+'_hi'] < df[x+'_lo']
    
    #hilight += df.index[up].to_list()
    #hilight += df.index[do].to_list()
    for h in hilight:
        select = df.index.str.contains(h)
        if not select.sum():
            continue
    
        if len(xerr) or len(yerr):
            ax.errorbar(x, y, data=df.loc[select] + 1, xerr=xerr[:, select], yerr=yerr[:, select], alpha=0.4, marker='o', ms=6, label=h, mfc='none', mew=2, zorder=len(df) + 1, ls='')
        else:
            ax.plot(x, y, "o", data=df.loc[select] + 1, ms=6, label=h, mfc='none', mew=2, zorder=len(df) + 1)            
            
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), frameon=False, handletextpad=0.4)

    if xlabel is None: xlabel = x
    if ylabel is None: ylabel = y
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
      ax.title.set_text(title)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    xmin = min(xmin, ymin)
    xmax = max(xmax, ymax)

    ax.plot([xmin, xmax], [xmin, xmax], "--k", lw=0.5)

    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")

    return R, pval    


def mRNA_miRNA_overview_plot(adata, gridsize=30, bins_mrna=100, bins_mirna=100, title="overview"):
    import matplotlib.pyplot as plt
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=False, figsize=(12, 3))
    fig.suptitle(title)
    ax1.set_xscale('log')
    ax2.set_xscale('log')
    y_mrna = adata.obs['n_coding_counts']
    ly_mrna = np.log10(y_mrna)
    bins = 10**np.linspace(ly_mrna.min(), ly_mrna.max(), bins_mrna)
    ax1.hist(y_mrna, bins=bins, alpha=1, label="coding genes")
    ax1.set_xlabel("mRNA UMIs/cell")
    ax1.set_ylabel("number of cells")

    y_mirna = adata.obs['n_miRNA_counts']
    ly_mirna = np.log10(y_mirna)
    bins = 10**np.linspace(ly_mirna.min(), ly_mirna.max(), bins_mirna)
    ax2.hist(y_mirna, bins=bins, alpha=1, label="mature miRNA")
    ax2.set_xlabel("miRNA UMIs/cell")
    ax2.set_ylabel("number of cells")

    ax1.spines.right.set_visible(False)
    ax1.spines.top.set_visible(False)
    ax2.spines.right.set_visible(False)
    ax2.spines.top.set_visible(False)
    ax3.spines.right.set_visible(False)
    ax3.spines.top.set_visible(False)

    from matplotlib.colors import LogNorm
    ax3.set_aspect(1)
    cm = ax3.hexbin(y_mrna,y_mirna, xscale='log', yscale='log', linewidths=0, gridsize=gridsize, mincnt=1, cmap="viridis") #, norm=LogNorm()
    plt.colorbar(cm, shrink=0.5, label="number of cells")
    ax3.set_xlabel("mRNA UMIs/cell")
    ax3.set_ylabel("miRNA UMIs/cell")

    from matplotlib.ticker import ScalarFormatter
    for axis in [ax1.xaxis, ax2.xaxis, ax3.xaxis, ax3.yaxis]:
        formatter = ScalarFormatter()
        formatter.set_scientific(False)
        axis.set_major_formatter(formatter)

    fig.tight_layout()
    print(f"Thus, we keep {len(adata)} cells with a median protein coding mRNA count of {int(y_mrna.median())} and miRNA count of {int(y_mirna.median())} for further analysis.")


def plot_cell_metrics(adata, x="coding", y="mt", gridsize=30, bins_x=100, bins_y=100, title="per-cell metrics"):
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    fig.suptitle(title)
    ax1.hist(adata.obs[f'pct_{x}'], bins=bins_x, label=f"{x} genes")
    ax1.hist(adata.obs[f'pct_{y}'], bins=bins_y, label=f"{y} genes")

    ax1.legend(ncol=2, frameon=False, loc='lower right', bbox_to_anchor=(1.0, 1.0))
    ax1.set_xlabel("% of counts")
    ax1.set_xlim(0, 100)
    ax1.set_ylabel("number of cells")

    from matplotlib.colors import LogNorm
    ax2.set_aspect(1)
    cm = ax2.hexbin(adata.obs[f'n_{x}_counts'] + 1, adata.obs[f'n_{y}_counts'] + 1, xscale='log', yscale='log', linewidths=0, gridsize=gridsize, mincnt=1, cmap="viridis") # , norm=LogNorm()
    plt.colorbar(cm, shrink=0.5, label="number of cells")
    ax2.set_xlabel(f"{x} gene counts")
    ax2.set_ylabel(f"{y} gene counts")

    ax1.spines.right.set_visible(False)
    ax1.spines.top.set_visible(False)
    ax2.spines.right.set_visible(False)
    ax2.spines.top.set_visible(False)

    fig.tight_layout()
    return fig


def plot_reads_per_UMI(adata, ax=None):
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots()

    n_gene_reads = np.array(adata.layers['exonic_reads'].sum(axis=1))[:,0]
    n_gene_UMIs = np.array(adata.layers['exonic_counts'].sum(axis=1))[:,0]

    ratio = n_gene_reads / n_gene_UMIs

    ax.hist(
        ratio,
        bins=50, #np.linspace(pcr_coding.min(), np.percentile(pcr_coding, 99), 100),
        alpha=1.0, label="reads per UMI")

    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    
    return ratio

