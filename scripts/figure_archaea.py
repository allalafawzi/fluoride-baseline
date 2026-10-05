#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Draw the four-panel archaeal baseline figure from baseline_archaea.json:
prevalence per phylum, what the sampling design changes, prevalence against
genome size, and the size-matched comparison between phyla."""
import json, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

R = json.load(open(sys.argv[1]))
OUT = sys.argv[2] if len(sys.argv) > 2 else "figure_archaea"

# Validated palette (six controls, light mode, surface #fcfcfb):
#   teal #00998c / brick #b8391f - CVD gap 12.3 (deutan), 26.3 in normal vision.
TEAL, BRICK = "#00998c", "#b8391f"
INK, MUTED, FAINT, RULE = "#101d1b", "#5d6f6a", "#8b9a95", "#dfe6e3"
SURF = "#fcfcfb"
REDUCED = 1.5e6         # threshold for a "reduced genome", in base pairs

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10.5, "figure.facecolor": SURF,
    "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": RULE, "axes.labelcolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.family": ["DejaVu Sans"],
})
fig = plt.figure(figsize=(12.2, 9.4))
gs = fig.add_gridspec(2, 2, width_ratios=[1.32, 1], height_ratios=[1.18, 1],
                      hspace=0.34, wspace=0.30)

def dot_color(t):
    return BRICK if (t == t and t < REDUCED) else TEAL

def dots(ax, items, lab_key="phylum", n_key="N"):
    """Points + horizontal intervals, marks >= 8 px, thin bars."""
    y = np.arange(len(items))[::-1]
    for yy, d in zip(y, items):
        c = dot_color(d.get("median_size", float("nan")))
        ax.plot([100 * d["lo"], 100 * d["hi"]], [yy, yy], color=c, lw=2,
                solid_capstyle="round", alpha=.55, zorder=2)
        ax.plot([100 * d["prop"]], [yy], "o", ms=7.5, color=c,
                mec=SURF, mew=1.4, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{d[lab_key]}" for d in items], fontsize=8.6, color=INK)
    for yy, d in zip(y, items):
        ax.text(-1.4, yy, f"{d[n_key]}", ha="right", va="center",
                fontsize=7.4, color=FAINT, transform=ax.get_yaxis_transform(
                    which="grid") if False else ax.transData)
    ax.grid(axis="x", color=RULE, lw=.7, zorder=0)
    ax.set_axisbelow(True)

# ---------------------------------------------------------------- A -------
ax = fig.add_subplot(gs[0, 0])
pp = [d for d in R["by_phylum"] if d["N"] >= 10]
dots(ax, pp)
ax.axvspan(85, 101, color=FAINT, alpha=.13, zorder=0)
ax.text(93, len(pp) - 0.35, "bacteria\n> 85 %", ha="center", va="top",
        fontsize=7.6, color=MUTED, style="italic")
b = R["raw_proportion"]; c = R["mean_by_phylum"]
ax.axvline(100 * c["estimate"], color=INK, lw=1, ls="--", alpha=.55, zorder=1)
ax.text(100 * c["estimate"] + 1.6, len(pp) - 0.42,
        f"mean-by-phylum\n{100*c['estimate']:.0f} %", fontsize=7.6, color=INK,
        va="top")
ax.set_xlim(-6, 101); ax.set_xlabel("genomes carrying Fluc/CrcB (%)")
ax.set_title("A - Prevalence by phylum, balanced sample",
             loc="left", fontweight="bold", color=INK)
ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=7.5, color=TEAL,
                          mec=SURF, mew=1.2, label="median genome $\\geq$ 1.5 Mb"),
                   Line2D([], [], marker="o", ls="", ms=7.5, color=BRICK,
                          mec=SURF, mew=1.2, label="median genome < 1.5 Mb")],
          fontsize=7.6, frameon=False, loc="center right",
          bbox_to_anchor=(1.0, .40), handletextpad=.4)
s = R["threshold_control"]
if s:
    ax.text(0.985, 0.015,
            f"threshold control: without a calibrated threshold,\n"
            f"{s['gained_by_relaxing']} genome gained out of {s['n_common']}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.4,
            color=MUTED, bbox=dict(fc=SURF, ec=RULE, lw=.7, pad=3.5), zorder=5)

# ---------------------------------------------------------------- B -------
ax = fig.add_subplot(gs[0, 1])
ests = [("unbalanced sample\n(1,510 genomes, 246 families)", .499, .473, .524, BRICK),
        ("balanced sample\n(1,259 genomes, 591 families)", b["prop"], b["lo"], b["hi"], TEAL),
        ("balanced, mean-by-phylum", c["estimate"], c["lo"], c["hi"], TEAL)]
y = np.arange(len(ests))[::-1]
for yy, (lab, p, lo, hi, col) in zip(y, ests):
    ax.plot([100 * lo, 100 * hi], [yy, yy], color=col, lw=2.4, alpha=.55,
            solid_capstyle="round", zorder=2)
    ax.plot([100 * p], [yy], "o", ms=9, color=col, mec=SURF, mew=1.5, zorder=3)
    ax.text(100 * p, yy + .21, f"{100*p:.1f} %", ha="center", fontsize=9,
            color=INK, fontweight="bold")
ax.set_yticks(y); ax.set_yticklabels([e[0] for e in ests], fontsize=8.2, color=INK)
ax.axvspan(85, 101, color=FAINT, alpha=.13, zorder=0)
ax.text(93, y[0] + .55, "bacteria > 85 %", ha="center", fontsize=7.8,
        color=MUTED, style="italic")
ax.annotate("", xy=(14, y[1] + .42), xytext=(49, y[0] - .42),
            arrowprops=dict(arrowstyle="->", color=INK, lw=1.1,
                            connectionstyle="arc3,rad=.25"))
ax.text(31, (y[0] + y[1]) / 2 + .05, "the sampling design,\non its own",
        fontsize=7.8, color=INK, ha="center", style="italic")
ax.set_xlim(0, 101); ax.set_ylim(-.75, len(ests) - .1)
ax.set_xlabel("genomes carrying Fluc/CrcB (%)")
ax.grid(axis="x", color=RULE, lw=.7); ax.set_axisbelow(True)
ax.set_title("B - What the sampling design changes",
             loc="left", fontweight="bold", color=INK)

# ---------------------------------------------------------------- C -------
ax = fig.add_subplot(gs[1, 0])
pt = R["by_size"]
x = np.arange(len(pt))
lab = [f"< 1.0" if d["max_mb"] == 1.0 else
       ("$\\geq$ 3.5" if d["max_mb"] > 90 else f"{d['min_mb']:.1f}-{d['max_mb']:.1f}") for d in pt]
for xx, d in zip(x, pt):
    col = BRICK if d["max_mb"] <= 1.5 else TEAL
    ax.plot([xx, xx], [100 * d["lo"], 100 * d["hi"]], color=col, lw=2,
            alpha=.55, solid_capstyle="round", zorder=2)
    ax.plot([xx], [100 * d["prop"]], "o", ms=8, color=col, mec=SURF, mew=1.4, zorder=3)
    ax.text(xx, 100 * d["hi"] + 2.6, f"{100*d['prop']:.1f} %", ha="center",
            fontsize=8, color=INK)
    ax.text(xx, -5.2, f"n={d['N']}", ha="center", fontsize=7.2, color=FAINT)
ax.plot(x, [100 * d["prop"] for d in pt], color=MUTED, lw=1, alpha=.4, zorder=1)
ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=8.4)
ax.set_xlabel("genome size (Mb)"); ax.set_ylabel("Fluc/CrcB (%)")
ax.set_ylim(-9, 78); ax.set_xlim(-.6, len(pt) - .4)
ax.grid(axis="y", color=RULE, lw=.7); ax.set_axisbelow(True)
m = R["size_phylum_model"]
ax.set_title("C - Genome size weighs heavily", loc="left",
             fontweight="bold", color=INK)
ax.text(.03, .95, f"$\\chi^2$ = {m['chi2_size']:.0f} (1 df)", transform=ax.transAxes,
        fontsize=8.4, color=MUTED, va="top")

# ---------------------------------------------------------------- D -------
ax = fig.add_subplot(gs[1, 1])
ap = R["size_matched"]
tm = {d["phylum"]: d.get("median_size", float("nan")) for d in R["by_phylum"]}
for d in ap:
    d["median_size"] = tm.get(d["phylum"], float("nan"))
dots(ax, ap)
ax.set_xlim(-6, 78); ax.set_xlabel("genomes carrying Fluc/CrcB (%)")
ax.set_title("D - ...but it does not explain the clade", loc="left",
             fontweight="bold", color=INK)
# The size restriction describes the dots, that is the size-matched comparison
# plotted here. The chi-square describes the logistic model, which is fitted on
# all 1259 genomes. The two are labelled separately so that the restriction is
# not read as a condition of the test.
ax.text(.97, .06, f"dots: size-matched comparison, 1.0-2.5 Mb\n"
                  f"test: all {R['sample']['N']} genomes, phylum after size\n"
                  f"$\\chi^2$ = {m['chi2_phylum_after_size']:.0f} "
                  f"({m['df_phylum']} df reported, "
                  f"{m.get('df_phylum_effective', m['df_phylum'])} effective)",
        transform=ax.transAxes, ha="right", fontsize=7.2, color=MUTED,
        bbox=dict(fc=SURF, ec=RULE, lw=.7, pad=3.5))

fig.suptitle("Fluc/CrcB fluoride exporter in the archaea - "
             f"{R['sample']['N']} genomes, {R['sample']['n_families']} families, "
             f"{R['sample']['n_phyla']} phyla",
             fontsize=12.5, fontweight="bold", color=INK, y=.977)
icc = R["measured_icc"]
fig.text(.5, .012,
         f"One genome per family, families in reproducible pseudo-random order | "
         f"95 % Wilson intervals | figures on the left = counts | "
         f"measured ICC: family {icc['family']:.2f}, phylum {icc['phylum']:.2f}",
         ha="center", fontsize=7.6, color=FAINT)
fig.subplots_adjust(left=.135, right=.985, top=.915, bottom=.075)
fig.savefig(f"{OUT}.png", dpi=200)
fig.savefig(f"{OUT}.pdf")
print(f"[ok] {OUT}.png", file=sys.stderr)
