import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np

def cell_hist(adata, key="n_reads", bins=100, xlog=True, ylog=False, ax=None, **kw):
    df = adata.obs
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

    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)


def loglog_knee(adata, key="n_counts", ax=None, min_cells=500, title=None, debug=False):
    """
    we will plot the number of cell barcodes (CBs) on the y-axis which have a UMI 
    count larger than a sliding threshold, which is plotted on the x-axis
    due to the dynamic ranges of both dimensions we plot this double-logarithmic
    
      - for low cutoffs (left side of plot) we will first see a large number of cell barcodes 
          which are background (encapsulation of free floating RNA)
      - eventually, the number of CBs begins to drop, that's often the regime of
          sequencing or synthesis errors in the CBs, think of it as a "shadow" of the real
          cells shifted to the left by ~2 orders of magnitude (error rate ~1%).
      - then we get to a region where the slope is less negative, or even approximates a 
          plateau (if cells are very homogeneous, e.g. cell culture). These are the captured
          transcriptomes of real cells
      - at very high cutoffs, we start seeing a rapid drop in cells that have that many UMIs.
          This is to be expected.
    
    This function also tries to autodetect a UMI cutoff value that would separate 
    properly captured cellular transcriptomes by finding the peak of the second derivate at 
    the inflection point between "plateau" of cells and background further on the left.
    If you want to see the details on this, pass `debug=True`. Note that when cell numbers get 
    low on the high end of UMI cutoffs, the derivative estimation from the spline fit gets noisy.
    To mask this noise, we employ a min_cells cutoff.
    """
    df = adata.obs
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))

    UMI = np.array(df[key].values, dtype=int)
    if not len(UMI):
        return

    UMIsum = UMI.sum()
    # bincount is like a histogram of bin-size 1 and always starting at 0
    n_CBs = len(UMI)
    y = np.array(n_CBs - np.bincount(UMI).cumsum(), dtype=np.float32)

    import scipy.interpolate
    x = np.arange(len(y)) + 1
    ly = np.log10(y + 1)
    lx = np.log10(x)
    # linear interpolation on the log-scale
    f = scipy.interpolate.interp1d(lx, ly)
    # bicubic spline approximation of the linear interpolation,
    # with constant extensions beyond the boundary. 
    # We use the lin-approx first to get an evenly spaced sampling 
    # of the data on log(x), which makes everything better...
    x = np.arange(0, lx.max(), 0.01)
    spl = scipy.interpolate.UnivariateSpline(x, f(x), ext=3, k=5, s=0.01)
    # mask noisy estimates of second derivative when we reach low cell count
    valid = spl(x) > np.log10(min_cells)
    d2_spl = spl.derivative(2)
    d3_spl = spl.derivative(3)
    dd = np.where(valid, d2_spl(x), 0)
    ddd = np.where(valid, d3_spl(x), 0)

    # we want to find roots of the third derivative as possible inflection points
    from scipy.interpolate import PPoly, splrep
    ppoly = PPoly.from_spline(splrep(x[valid], ddd[valid]))
    roots = ppoly.roots(extrapolate=False)
    x_cuts = []
    for r in roots:
        # we have a local maximum in second derivative. 
        # That's what we're looking for.
        if d2_spl(r) > 0: 
            x_cuts.append(r)
    
    # keep the right-most inflection point with 
    # correct curvature: -> \_ NOT =\
    if len(x_cuts):
        x_cut = int(10**x_cuts[-1])
    else:
        x_cut = np.nan

    # reasonable UMI cutoff value from peak of second derivative
    # x_cut = int(10**x[dd.argmax()])
    n_cells = (UMI >= x_cut).sum()
    n_UMI = ((UMI >= x_cut) * UMI).sum()
    if debug:
        ax.plot(x, f(x), label='linear interpolation of ')
        ax.plot(x, spl(x), label='spline interpolation of linear interpolation')
        ax.plot(x, spl.derivative(1)(x), label='first derivative of spline')
        # ax.plot(x, spl.derivative(2)(x), label="2nd")
        ax.plot(x, dd, label='masked, second derivative of spline')
        ax.plot(x, ddd, label='masked, third derivative of spline')
        ax.axvline(roots[-1])
    else:
        ax.loglog(
            y,
            label=f"{key} (total={UMIsum/1e6:.1f}M)",
            alpha=0.75,
            lw=3,
        )
        
        ax.axvline(x_cut, color='r', lw=0.5, label=f'suggested UMI cutoff:\n  UMI >= {x_cut}\n  n_cells={n_cells}\n  n_UMI={n_UMI/1e6:.1f}M ({100.0 * n_UMI/UMIsum:.1f} %)')
    
    ax.set_ylabel("no. cells above cutoff")
    ax.set_xlabel(f"{key} > x")

    ax.legend(loc='lower left', fancybox=False, markerscale=0.5)
    ax.grid(axis="y")

    if fig:
        
        if not title:
            sample_name = adata.uns.get("sample_name", "no-sample-name")
            fig.suptitle(sample_name)
        else:
            fig.suptitle(title)

        fig.tight_layout()

    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)

    return fig, x_cut


def dge_stats(adata, hist_kw={}, fig_kw={'figsize': (8,5)}):
    # cell metrics plot
    fig, axes = plt.subplots(2, 2, **fig_kw)
    axes = np.array(axes).ravel()
    cell_hist(adata, key="n_reads", ax=axes[0], **hist_kw)
    cell_hist(adata, key="n_counts", ax=axes[1], **hist_kw)
    cell_hist(adata, key="n_genes", ax=axes[2], **hist_kw)
    cell_hist(adata, key="reads_per_counts", xlog=False, ax=axes[3], **hist_kw)
    fig.subplots_adjust(hspace=0.4)
    # fig.tight_layout()

    return fig

def reads_per_UMI(adata, ax=None):
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots()

    n_gene_reads = np.array(adata.layers['exonic_reads'].sum(axis=1))[:,0]
    n_gene_UMIs = np.array(adata.layers['exonic_counts'].sum(axis=1))[:,0]

    ratio = n_gene_reads / n_gene_UMIs

    ax.hist(
        ratio,
        bins=50, #np.linspace(pcr_coding.min(), np.percentile(pcr_coding, 99), 100),
        alpha=1.0, label="exonic reads per UMI")

    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    
    return ratio

def scatter(df, x, y, hilight=[], cutoff=10, xlabel=None, ylabel=None, cmap='tab10', title='counts', ax=None, default_figsize=(5,5)):
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


def mRNA_miRNA_overview(adata, gridsize=30, bins_mrna=100, bins_mirna=100, title="overview"):
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


def cell_metrics(adata, x="coding", y="mt", gridsize=30, bins_x=100, bins_y=100, title="per-cell metrics"):
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


