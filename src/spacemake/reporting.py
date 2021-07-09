import os
import sys
import logging
import pandas as pd
import numpy as np


def count_dict_collapse_misc(
    counts, misc_thresh=0.01, total=1, add_up=None, sig_intact=None
):
    out_counts = {}
    out_frac = {}

    misc = 0
    sum = 0
    if sig_intact is not None:
        complete = ",".join(sig_intact)
        everything = set(sig_intact)
    else:
        complete = None
        everything = set()

    def relkey(key):
        if sig_intact is None:
            return key

        if key == complete:
            return "complete"

        obs = set(key.split(","))
        there = obs & everything
        extra = obs - everything
        missing = everything - obs

        if len(missing) <= len(there):
            res = "missing_" + ",".join(sorted(missing))
        else:
            res = "only_" + ",".join(sorted(there))
        if extra:
            res += "_extra_" + ",".join(sorted(extra))

        return res

    for key, n in sorted(counts.items()):
        key = relkey(key)
        sum += n
        f = n / float(total)
        if f < misc_thresh:
            misc += n
        else:
            out_counts[key] = n
            out_frac[key] = f

        if misc > 0:
            out_counts["misc"] = misc
            out_frac["misc"] = misc / float(total)

    if add_up is None:
        other = total - sum
    else:
        other = total - counts[add_up]

    if other > 0:
        out_counts["NA"] = other
        out_frac["NA"] = other / float(total)
    return out_counts, out_frac


def count_dict_out(counts, title, misc_thresh=0.01, total=1, **kw):
    print(f"### {title}")
    out_counts, out_frac = count_dict_collapse_misc(counts, misc_thresh, total, **kw)
    for key in sorted(out_counts.keys()):
        print(f"{key}\t{out_counts[key]}\t{out_frac[key]:.3f}")


def to_hist(d, normed=True):
    x = np.array(list(d.keys()))
    x0 = x.min()
    x1 = x.max() + 1
    counts = np.zeros(x1, dtype=np.float32)

    for i in x:
        counts[i] = d[i]

    n = counts.sum()
    if normed:
        counts /= n

    return counts, n


def donut_plot(
    ax, data, sa=10, explode=None, colors=None, labels=None, title="", cmap="tab20"
):
    import matplotlib.pyplot as plt

    if labels is None:
        labels = sorted(data.keys())

    counts = [data.get(k, 0) for k in labels]

    if colors is None:
        colors = list(plt.cm.get_cmap(cmap)(np.linspace(0, 1, len(labels))))

    wedges, texts = ax.pie(
        counts,
        wedgeprops=dict(width=0.5),
        startangle=sa,
        explode=explode,
        colors=colors,
    )

    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.5)
    kw = dict(arrowprops=dict(arrowstyle="-"), bbox=bbox_props, zorder=0, va="center")
    c = np.array(counts)
    pcts = 100.0 * c / float(c.sum())
    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1) / 2.0 + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        if pcts[i] > 0:
            ax.text(x * 0.75, y * 0.75, f"{pcts[i]:.1f}", horizontalalignment="center")
            ax.annotate(
                labels[i],
                xy=(x, y),
                xytext=(1.4 * np.sign(x), 1.4 * y),
                horizontalalignment=horizontalalignment,
                **kw,
            )

    if title:
        ax.set_title(title)

    return labels, colors


def approximate(intvalue):
    suffixes = {9: "G", 6: "M", 3: "k", 0: ""}
    dec = int(np.floor(np.log10(intvalue) / 3)) * 3
    x = np.round(intvalue / 10 ** dec, decimals=2)
    return f"{x:.2f} {suffixes.get(dec, '?')}"


def len_plot(
    ax,
    data,
    labels=None,
    colors=None,
    xlabel="aligned bases",
    ylabel="fraction",
    title="type",
    cmap="tab20",
    min_count=10,
    cumulative=False,
    legend=True,
):
    import matplotlib.pyplot as plt

    if labels is None:
        labels = sorted(data.keys())

    if colors is None:
        colors = plt.cm.get_cmap(cmap)(np.linspace(0, 1, len(labels)))

    color_dict = {}
    for cig_type, color in zip(labels, colors):
        color_dict[cig_type] = color

        if not cig_type in data:
            continue
        ml, n = to_hist(data[cig_type], normed=True)
        if n < min_count:
            continue

        x = np.arange(len(ml))
        y = ml.cumsum() if cumulative else ml

        ax.step(
            x,
            y,
            where="mid",
            label=f"{cig_type} ({approximate(n)})",
            color=color,
            lw=2,
            solid_capstyle="round",
        )

    if cumulative:
        ax.axhline(0.5, lw=0.5, ls="dashed", color="k")

    if legend:
        ax.legend(title=title, bbox_to_anchor=(0.5, 1.05), loc="lower center", ncol=2)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return color_dict


# def make_colors_explode(labels, cmap="Blues", hilight="bead-related", hicolor="red"):
#     import matplotlib.pyplot as plt
#     ex = np.zeros(len(labels))
#     colors = list(plt.get_cmap(cmap)(np.linspace(0.2, 0.8, len(labels))))
#     try:
#         i = labels.index(hilight)
#     except ValueError:
#         pass
#     else:
#         ex[i] = 0.1
#         colors[i] = hicolor
#     return ex, colors
