#!/usr/bin/env python
"""Supplementary: matched-null histograms for all 5 runs (rows) x anchors /
regulators / targets (columns). Page 1 = enrichment (log x), page 2 =
coefficient z. Real value as a line in the set colour; empirical p and n
in each panel. Writes fig_s_null_histograms.pdf (2 pages) and one PNG per page.

  uv run python figures/fig_s_null_histograms.py
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

from fig_common import (
    ENR_FLOOR, FIG_WIDTH, NULL_GREY, OUT, RUNS, SET_COLOURS,
    SET_LABELS, SETS, STATS, emp_p, fmt_p, load_perperm, null_and_real,
    plain_log_axis, run_label, setup_style,
)


def page(stat: str):
    fig, axes = plt.subplots(len(RUNS), len(SETS), figsize=(FIG_WIDTH, 8.2),
                             squeeze=False)
    for i, run in enumerate(RUNS):
        df = load_perperm(*run)
        for j, s in enumerate(SETS):
            ax = axes[i][j]
            real, null = null_and_real(df, s, stat)
            p = emp_p(real, null)
            if stat == "enrichment":
                x = np.clip(null, ENR_FLOOR, None)
                bins = np.geomspace(ENR_FLOOR, max(x.max(), real) * 1.05, 36)
            else:
                x = null
                lo, hi = min(null.min(), real), max(null.max(), real)
                bins = np.linspace(lo - 0.1, hi + 0.1, 31)
            ax.hist(x, bins=bins, color=NULL_GREY, edgecolor="white", linewidth=0.3)
            ax.axvline(real, color=SET_COLOURS[s], lw=1.3)
            if stat == "enrichment":
                plain_log_axis(ax, "x", [0.1, 1, 10, 100], floor_label=ENR_FLOOR)
                ax.set_xlim(ENR_FLOOR / 1.25, 150)
            ax.set_title(f"real = {real:.2f}\n{fmt_p(p)}, n = {null.size}",
                         fontsize=6.5, linespacing=1.25)
            n_low = int(np.sum(null < ENR_FLOOR)) if stat == "enrichment" else 0
            if n_low:
                ax.annotate(f"{n_low} ≤ {ENR_FLOOR:g}", xy=(ENR_FLOOR * 1.06, 0.97),
                            xycoords=("data", "axes fraction"), xytext=(4, 0),
                            textcoords="offset points", va="top", fontsize=6)
            if i == 0:
                ax.text(0.5, 1.45, SET_LABELS[s], transform=ax.transAxes, ha="center",
                        fontsize=8, fontweight="bold", color=SET_COLOURS[s])
            if j == 0:
                ax.set_ylabel("Permutations")
                ax.text(-0.3, 0.5, run_label(*run).replace(", ", ",\n")
                        .replace("thyroid disease", "thyroid\ndisease")
                        .replace(" (childhood", "\n(childhood"),
                        transform=ax.transAxes, ha="right", va="center",
                        fontsize=7, fontweight="bold")
            if i == len(RUNS) - 1:
                ax.set_xlabel(STATS[stat])
    fig.tight_layout(h_pad=0.9, w_pad=0.6)
    return fig


def main():
    setup_style()
    OUT.mkdir(parents=True, exist_ok=True)
    with PdfPages(OUT / "fig_s_null_histograms.pdf") as pdf:
        for k, stat in enumerate(STATS, start=1):
            fig = page(stat)
            pdf.savefig(fig, bbox_inches="tight")
            png = OUT / f"fig_s_null_histograms_p{k}_{stat}.png"
            fig.savefig(png, bbox_inches="tight")
            plt.close(fig)
            print(f"Wrote {png}")
    print(f"Wrote {OUT / 'fig_s_null_histograms.pdf'}")


if __name__ == "__main__":
    main()
