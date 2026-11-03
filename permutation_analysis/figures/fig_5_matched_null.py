#!/usr/bin/env python
"""Fig. 5: S-LDSC enrichment of coloc-anchored gene sets against the corrected
matched-anchor null.

  A  ASTCO Rest: null enrichment histograms (200 perms) for anchors,
     regulators, targets; real value as a line in the set colour.
  B  All 5 runs: real value vs null median and 2.5-97.5% range, for
     enrichment (log axis) and coefficient z. Filled dot = empirical p < 0.05.
  C  Fraction of the real regulator / target sets recovered by matched
     pseudo-anchors at Rest (median, 95% range over 500 draws; k-NN matching,
     from scripts/matching_check_knn.py).

  uv run python figures/fig_5_matched_null.py
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

from fig_common import (
    COND_NAMES, DISEASE_NAMES, ENR_FLOOR, ENR_TICKS, FIG_WIDTH, NULL_GREY, RUNS,
    SET_COLOURS, SET_LABELS, SETS, SOURCE, STATS, emp_p, fmt_p, load_perperm,
    null_and_real, panel_letter, plain_log_axis, run_label, run_stats, save,
    setup_style,
)

PANEL_A_RUN = ("ASTCO_Ferreira_30929738_1", "Rest")
ENR_XLIM = (ENR_FLOOR / 1.25, 150)


HIST_BINS = np.geomspace(ENR_FLOOR, 150, 61)


def log_hist(ax, null, real, colour):
    clipped = np.clip(null, ENR_FLOOR, None)
    ax.hist(clipped, bins=HIST_BINS, color=NULL_GREY, linewidth=0)
    ax.axvline(real, color=colour, lw=1.4)
    plain_log_axis(ax, "x", [0.1, 1, 10, 100], floor_label=ENR_FLOOR)
    ax.set_xlim(*ENR_XLIM)
    n_low = int(np.sum(null < ENR_FLOOR))
    if n_low:
        ax.annotate(f"{n_low} ≤ {ENR_FLOOR:g}", xy=(ENR_FLOOR * 1.06, 0.97),
                    xycoords=("data", "axes fraction"), xytext=(4, 0),
                    textcoords="offset points", va="top", fontsize=6)


def panel_a(axes):
    df = load_perperm(*PANEL_A_RUN)
    for ax, s in zip(axes, SETS):
        real, null = null_and_real(df, s, "enrichment")
        log_hist(ax, null, real, SET_COLOURS[s])
        ax.set_title(f"{SET_LABELS[s].split(' (')[0]}\n{fmt_p(emp_p(real, null))}, n = {null.size}",
                     fontsize=7.5, color=SET_COLOURS[s], linespacing=1.3)
        ax.set_xlabel("Enrichment")
    axes[0].set_ylabel("Permutations")
    axes[1].text(0.5, 1.33, run_label(*PANEL_A_RUN),
                 transform=axes[1].transAxes, ha="center", fontsize=8)


def panel_b(ax_enr, ax_z, stats: pl.DataFrame):
    ypos, labels, y = {}, [], 0.0
    group_y = {}
    for s in SETS:
        group_y[s] = y
        y += 0.9
        for run in RUNS:
            ypos[(s, run)] = y
            labels.append((y, run_label(*run)))
            y += 1
        y += 0.5
    for ax, stat in ((ax_enr, "enrichment"), (ax_z, "coefficient_z")):
        sub = stats.filter(pl.col("statistic") == stat)
        for r in sub.iter_rows(named=True):
            yy = ypos[(r["set"], (r["study"], r["condition"]))]
            lo, hi, med, real = r["null_q025"], r["null_q975"], r["null_median"], r["real"]
            if stat == "enrichment":
                clipped = lo < ENR_FLOOR
                lo, med = max(lo, ENR_FLOOR), max(med, ENR_FLOOR)
            else:
                clipped = False
            ax.plot([lo, hi], [yy, yy], color="#d0d0d0", lw=3.2,
                    solid_capstyle="butt", zorder=1)
            if clipped:
                ax.plot(lo, yy, marker="<", ms=3, color=NULL_GREY, zorder=2)
            ax.plot([med, med], [yy - 0.3, yy + 0.3], color="#6e6e6e", lw=0.9, zorder=2)
            c = SET_COLOURS[r["set"]]
            filled = r["empirical_p"] < 0.05
            ax.plot(real, yy, "o", ms=4, mec=c, mew=1.0,
                    mfc=c if filled else "white", zorder=3)
        ax.set_ylim(y - 0.3, -0.4)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.set_xlabel(STATS[stat])
        if stat == "enrichment":
            plain_log_axis(ax, "x", ENR_TICKS, floor_label=ENR_FLOOR)
            ax.set_xlim(*ENR_XLIM)
            ax.axvline(1, color="#bbbbbb", lw=0.5, ls=":", zorder=0)
        else:
            ax.axvline(0, color="#bbbbbb", lw=0.5, ls=":", zorder=0)
            ax.set_xlim(-2.5, 5)
    for s in SETS:
        ax_enr.text(-0.02, group_y[s] + 0.05, SET_LABELS[s], transform=ax_enr.get_yaxis_transform(),
                    ha="right", va="center", fontsize=7.5, fontweight="bold",
                    color=SET_COLOURS[s])
    for yy, lab in labels:
        ax_enr.text(-0.02, yy, lab, transform=ax_enr.get_yaxis_transform(),
                    ha="right", va="center", fontsize=7)
    handles = [Line2D([], [], color="#d0d0d0", lw=3.2, label="Null 2.5–97.5%"),
               Line2D([], [], color="#6e6e6e", lw=0.9, marker="|", ms=6, ls="",
                      label="Null median"),
               Line2D([], [], color="k", marker="o", ms=4, ls="", label="Real, p < 0.05"),
               Line2D([], [], color="k", marker="o", ms=4, mfc="white", ls="",
                      label="Real, p ≥ 0.05")]
    ax_z.figure.legend(handles=handles, loc="upper center", ncol=4, frameon=False,
                       fontsize=7, handlelength=1.6,
                       bbox_to_anchor=(0.5, ax_z.get_position().y0 - 0.055))


def panel_c(ax):
    d = pl.read_csv(SOURCE / "matching_check_knn_draws.csv").filter(pl.col("method") == "knn")
    markers = {"AITD": "o", "ASTCO": "s"}
    for i, (s, col) in enumerate((("regulators", "reg_overlap_frac_of_real"),
                                  ("targets", "tgt_overlap_frac_of_real"))):
        for j, tr in enumerate(("AITD", "ASTCO")):
            x = d.filter(pl.col("study") == tr)[col].to_numpy()
            med, lo, hi = np.median(x), *np.quantile(x, [0.025, 0.975])
            xx = i + (j - 0.5) * 0.3
            ax.errorbar(xx, med, yerr=[[med - lo], [hi - med]], fmt=markers[tr],
                        color=SET_COLOURS[s], ms=4,
                        elinewidth=1, capsize=1.5)
    ax.set_xticks([0, 1], ["Regulators", "Targets"])
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of real set recovered")
    ax.set_title(f"Null overlap with real sets,\n{COND_NAMES['Rest']}", fontsize=7.5,
                 linespacing=1.3)
    handles = [Line2D([], [], color="#444444", marker=markers[t], ls="", ms=4,
                      label=DISEASE_NAMES[t])
               for t in ("AITD", "ASTCO")]
    ax.legend(handles=handles, frameon=False, fontsize=6.5, loc="upper center",
              handletextpad=0.2, bbox_to_anchor=(0.45, -0.2))


def main():
    setup_style()
    stats = run_stats()
    fig = plt.figure(figsize=(FIG_WIDTH, 5.6))
    gs = GridSpec(2, 1, figure=fig, height_ratios=[1.0, 2.55], hspace=0.42)
    top = gs[0].subgridspec(1, 5, width_ratios=[1, 1, 1, 0.35, 1.15], wspace=0.3)
    axa = [fig.add_subplot(top[0, i]) for i in range(3)]
    axc = fig.add_subplot(top[0, 4])
    bot = gs[1].subgridspec(1, 3, width_ratios=[1.45, 1.35, 1.1], wspace=0.12)
    ax_enr = fig.add_subplot(bot[0, 1])
    ax_z = fig.add_subplot(bot[0, 2])
    fig.subplots_adjust(bottom=0.1)
    panel_a(axa)
    panel_b(ax_enr, ax_z, stats)
    panel_c(axc)
    panel_letter(axa[0], "A", x=-0.3, y=1.3)
    panel_letter(axc, "C", x=-0.45, y=1.3)
    fig.text(0.02, ax_enr.get_position().y1 + 0.01, "B", fontsize=9, fontweight="bold")
    save(fig, "fig_5_matched_null")


if __name__ == "__main__":
    main()
