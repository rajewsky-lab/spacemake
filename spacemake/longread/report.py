import numpy as np
import logging
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def multi_row_barplots(dfs, row_labels, samples, attr, color_dict={}):
    import spacemake.longread.util as util

    fig, axes = plt.subplots(len(row_labels), sharex=True, sharey=True, figsize=(10, 8))
    for df, part, ax in zip(dfs, row_labels, axes):
        x, y = util.gather_data_from_overview(df, samples, attr)
        c = [color_dict[s] for s in samples]
        ax.bar(x, y, color=c)
        ax.set_ylabel(part)
        ax.set_xticks(x)
        ax.set_xticklabels(samples, rotation=45, horizontalalignment="right")

        for _x, _y in zip(x, y):
            ax.text(
                _x,
                0.5,
                f"{_y:.2f}",
                horizontalalignment="center",
                verticalalignment="center",
            )
    return fig, axes


def donut_plot(
    ax,
    labels,
    counts,
    sa=10,
    explode=None,
    colors=None,
    w=0.5,
    radius=1,
    label_inside=False,
    **kw,
):
    wedges, texts = ax.pie(
        counts,
        wedgeprops=dict(width=w * radius),
        startangle=sa,
        explode=explode,
        colors=colors,
        radius=radius,
        **kw,
    )
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.5)
    kw = dict(arrowprops=dict(arrowstyle="-"), bbox=bbox_props, zorder=0, va="center")
    c = np.array(counts)
    pcts = 100.0 * c / float(c.sum())
    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1) / 2.0 + p.theta1
        y = np.sin(np.deg2rad(ang)) * radius
        x = np.cos(np.deg2rad(ang)) * radius
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        # print(i, labels[i], pcts[i], x, y)
        if pcts[i] > 1.0:
            if label_inside:
                label = f"{labels[i]}\n{pcts[i]:.1f}"
            else:
                label = f"{pcts[i]:.1f}"

            # draw text inside the wedge
            ax.text(
                x * (radius - w / 2),  # * radius,
                y * (radius - w / 2),  # * radius,
                label,
                horizontalalignment="center",
                verticalalignment="center",
            )
            # draw annotation label outside, linked with solid line
            if not label_inside:
                ax.annotate(
                    labels[i],
                    xy=(x, y),
                    xytext=((radius + 0.15) * np.sign(x), y * 1.1),
                    horizontalalignment=horizontalalignment,
                    **kw,
                )


def make_colors_explode(
    labels, cmap="Blues", hilight="bead-related", hicolor="red", color_order=[]
):
    ex = np.zeros(len(labels))
    colors = list(plt.get_cmap(cmap)(np.linspace(0.2, 0.8, len(labels))))
    # assigned = ['gray'] * len(labels)
    # for i, l in enumerate(labels):
    #     if l in color_order:
    #         assigned[i] = colors[color_order.index(l)]
    # for j, name in enumerate(color_order):
    #     try:
    #         z = labels.index(name)
    #     except ValueError:
    #         pass
    #     else:
    #         colors[z], colors[j] = colors[j], colors[z]
    try:
        i = labels.index(hilight)
    except ValueError:
        pass
    else:
        ex[i] = 0.1
        colors[i] = hicolor
    return ex, colors


def plot_results(
    sig_counts,
    donut_labels,
    donut_counts,
    bead_normed_labels,
    bead_normed_counts,
    syn_rates,
    all_parts=["bead_start", "polyT"],
    fname="donuts.pdf",
    suptitle="",
):

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    [(ax1, ax2), (ax0, ax3)] = axes
    n_total = sig_counts["n_total"]

    if suptitle:
        fig.suptitle(f"{suptitle} (n_reads={n_total})")

    ax0.set_title("top 10 read signatures")
    ax0.set_xscale("log")

    # Hide the right and top spines
    ax0.spines["right"].set_visible(False)
    ax0.spines["top"].set_visible(False)

    # Only show ticks on the left and bottom spines
    sig_items = sorted(sig_counts.items(), key=lambda x: x[1])[-10:]
    sig_labels = [f"{x[0]} ({x[1]})" for x in sig_items]
    sig_cnt = [x[1] for x in sig_items]
    y_pos = np.arange(len(sig_labels))
    ax0.barh(y_pos, sig_cnt, height=0.2)
    ax0.set_yticks(y_pos)
    ax0.set_yticklabels(sig_labels)

    ax1.set_title("library overview")
    ax1.set_aspect("equal")
    ex, colors = make_colors_explode(
        donut_labels, cmap="GnBu", hilight="bead-related", hicolor="lawngreen"
    )
    donut_plot(ax1, donut_labels, donut_counts, sa=40, explode=ex, colors=colors)

    ax2.set_title("bead completeness")
    ax2.set_aspect("equal")
    ex, colors = make_colors_explode(
        bead_normed_labels,
        cmap="YlOrRd",
        hilight="complete",
        hicolor="seagreen",
        color_order=["missing_OP1", "missing_pT", "only_bead_start"],
    )
    donut_plot(
        ax2, bead_normed_labels, bead_normed_counts, sa=30, explode=ex, colors=colors
    )
    # estimate synthesis completion rates based on co-linear detection
    # of the constitutive parts of the capture oligos
    ax3.set_title("synthesis rates")
    ax3.spines["right"].set_visible(False)
    ax3.spines["top"].set_visible(False)
    x = range(1, len(syn_rates) + 1)
    ax3.barh(x, syn_rates[::-1])
    ax3.set_yticks(x)
    ax3.set_yticklabels([all_parts[i - 1] for i in x[::-1]])
    ax3.set_xlabel("synthesis completion")
    ax3.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax3.set_xlim(0, 1)

    fig.subplots_adjust(
        top=0.942, bottom=0.16, left=0.125, right=0.813, hspace=0.167, wspace=0.2
    )

    # fname = os.path.basename(sys.argv[1]).replace('.fq', '.report')
    logger = logging.getLogger("pb_annotate.plot_results")
    logger.info(f"saving donut-charts with annotation breakdown to '{fname}'")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(fname)

    return fig, axes


def obs_to_arrays(df):
    x = df["value"].values
    y = df["freq"].values
    return x, y


def plot_histograms(df, fname, parts=["bead_start", "OP1", "pT"], n_total=1):
    import seaborn as sns

    all_ends = df.query("attr == 'end'")["value"]
    max_x = int((np.percentile(all_ends, 95) + 1) / 100) * 100

    print(
        f"position 95th percentile would be {np.percentile(all_ends, 95)} -> max_x={max_x}"
    )

    n_parts = len(parts)
    fig, axes = plt.subplots(4, n_parts, figsize=(2 + n_parts * 3, 8))
    # axes = np.array(axes)
    # print(axes, n_parts)
    # first row: start coordinates, second row end coordinates, third row score distributions
    colors = plt.get_cmap("tab10")(np.arange(len(parts)))
    for ax_row, attr in zip(axes, ["start", "end", "score", "len"]):
        for part, ax, color in zip(parts, ax_row, colors):
            # ax.set_yscale("log")
            ctrl_x, ctrl_y = obs_to_arrays(
                df.query(
                    f"signature == 'anywhere' and oligo == '{part}' and attr == '{attr}'"
                )
            )
            df_sig = df.query(
                f"signature != 'anywhere' and oligo == '{part}' and attr == '{attr}'"
            )
            # print(df_sig)
            sig_x, sig_y = obs_to_arrays(df_sig)

            if attr == "len":
                ys = sig_y / n_total
                yc = ctrl_y / n_total
            else:
                ys = sig_y.cumsum() / n_total
                yc = ctrl_y.cumsum() / n_total

            ax.plot(
                sig_x,
                ys,
                "--",
                color=color,
                label=f"intact",
            )
            ax.plot(
                ctrl_x,
                yc,
                "-",
                color=color,
                label=f"all",
            )
            if attr == "start":
                ax.set_title(part)
            ax.legend(frameon=False)
            ax.set_xlabel(f"match {attr}")
            if attr != "len":
                ax.set_ylabel("cumulative rel. frequency")
                ax.set_ylim(0, 1)
            else:
                ax.set_ylabel("rel. frequency")
            if attr in ["start", "end"]:
                ax.set_xlim(0, max_x)

    fig.tight_layout()
    fig.savefig(fname)


def plot_edits_heatmap(ax, oname, oseq, edits, nmatch=1):
    N = len(oseq)
    counts = np.zeros((N, 6), dtype=float)
    x = np.arange(N)
    for i in x:
        for j, nn in enumerate("ACGT-"):
            if nn == oseq[i]:
                counts[i, j] = np.nan
            else:
                counts[i, j] = edits.iloc[i][oseq[i] + nn]

        is_ins = np.array([edits.iloc[i]["-" + ins] for ins in "ACGT"]).sum()
        counts[i, 5] = is_ins

    # print(counts, "non-finite entries", (~np.isfinite(counts)).sum())
    # print(counts / nmatch)
    im = ax.imshow(counts.T / nmatch, interpolation="none", cmap="viridis")
    plt.colorbar(im, ax=ax, shrink=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(list(oseq))
    ax.set_yticks(range(6))
    ax.set_yticklabels(list("ACGT-+"))


def plot_edits(df, fname, parts=["bead_start", "OP1", "pT"]):
    n_parts = len(parts)
    fig, axes = plt.subplots(n_parts, 2, figsize=(10, 1 + n_parts * 2))
    colors = plt.get_cmap("tab10")(np.arange(len(parts)))
    for part, (ax1, ax2), color in zip(parts, axes, colors):
        dfp = df.query(f"oligo == '{part}'")
        # print(dfp["pos"], dfp["fmatch"])
        ax1.step(dfp["pos"].values, dfp["fmatch"].values, where="mid", lw=2)
        ax1.set_xlabel(f"{part} pos. [nt]")
        ax1.set_ylabel("match frequency")
        #
        if len(dfp["seq"]):
            plot_edits_heatmap(
                ax2,
                part,
                dfp["seq"].iloc[0],
                dfp["ed_dict"],
                nmatch=df["nmatch"].iloc[0],
            )
            # ax2.title.set_text(f"{region} n={len(qnames)}")
            ax2.set_xlabel(f"{part} sequence")
        ax1.set_ylim(0, 1)

    fig.tight_layout()
    fig.savefig(fname)
