#!/usr/bin/env python
"""Supplementary forest plot: enrichment (top) and tau* (bottom) with 95% CIs.

One column per annotation (anchors, regulators, targets), one row per trait
(Fig. 2.2D order), conditions as colour. Each panel's x range covers the bulk
of the point estimates (x_range); CIs running past it are clipped with a
small arrowhead at the panel edge, and an estimate that is itself off-panel
is drawn as a larger arrowhead there. Filled: significant (*_sig columns:
BH FDR < 0.05 within expanded sets / within anchors, one-sided tau*,
negative-Prop._h2 sets never significant); open: not. Cells without an
estimate are left blank.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from hcommon import (FDR, FITTED, H_ANNOT_LABEL, H_ANNOTS, add_args, flagged_name,
                     load_diseases, load_long, resolve)
from common import (COND_COLORS, COND_LABEL, CONDS, PRINT_W, PT, axes_in,
                    color_ticklabels, panel_label, print_style, save, shade_groups)

fig4 = __import__("02_fig_4_heritability")

ROW_H = 0.17
DODGE = {"Rest": 0.25, "Stim8hr": 0.0, "Stim48hr": -0.25}
METRICS = [("enrichment", "enrichment_se", "enrichment_sig", 1.0, "Enrichment"),
           ("tau_star", "tau_star_se", "tau_star_sig", 0.0, "τ*")]


def x_range(df, est, ref):
    """Show the bulk of the point estimates: their 5-95% range widened by
    half its width on each side, always including the reference line."""
    lo, hi = df[est].quantile([0.05, 0.95])
    lo, hi = min(lo, ref), max(hi, ref)
    pad = 0.5 * (hi - lo)
    return lo - pad, hi + pad


def forest(ax, df, order, est, se, sig_col, ref, diseases, labels):
    n = len(order)
    ypos = {s: n - 1 - i for i, s in enumerate(order)}
    shade_groups(ax, order[::-1], diseases, label=False)
    lo_lim, hi_lim = x_range(df, est, ref)
    ax.axvline(ref, color="#8a8985", lw=0.5, ls="--", zorder=1)
    for _, r in df.iterrows():
        y = ypos[r.study] + DODGE[r.condition]
        c = COND_COLORS[r.condition]
        lo, hi = r[est] - 1.96 * r[se], r[est] + 1.96 * r[se]
        ax.plot([max(lo, lo_lim), min(hi, hi_lim)], [y, y], color=c, lw=0.6,
                solid_capstyle="butt", zorder=2)
        for edge, beyond, m in [(lo_lim, lo < lo_lim, "<"), (hi_lim, hi > hi_lim, ">")]:
            if beyond:
                ax.plot(edge, y, marker=m, color=c, markersize=2.6, mew=0, zorder=2)
        sig = bool(r[sig_col])
        face = c if sig else "white"
        if lo_lim <= r[est] <= hi_lim:
            ax.plot(r[est], y, marker="o", markersize=2.6, mew=0.6, color=c,
                    markerfacecolor=face, zorder=3)
        else:   # estimate itself off-panel: larger arrowhead at the edge
            edge, m = (lo_lim, "<") if r[est] < lo_lim else (hi_lim, ">")
            ax.plot(edge, y, marker=m, markersize=3.6, mew=0.6, color=c,
                    markerfacecolor=face, zorder=3, clip_on=False)
    ax.set_xlim(lo_lim, hi_lim)
    ax.set_ylim(-0.5, n - 0.5)
    ax.set_yticks(range(n))
    ax.tick_params(axis="y", length=0)
    if labels:
        ax.set_yticklabels([flagged_name(diseases, s) for s in order[::-1]])
        color_ticklabels(ax.get_yticklabels(), order[::-1], diseases)
    else:
        ax.tick_params(axis="y", labelleft=False)
    ax.spines["left"].set_visible(False)
    ax.locator_params(axis="x", nbins=3)


def draw_figure(long, diseases):
    order = fig4.trait_order(long, diseases)
    n = len(order)
    side = n * ROW_H
    names_w, panel_w, gap = 1.62, 1.36, 0.23
    W = PRINT_W
    bottom, mid_gap, top = 0.62, 0.62, 0.42
    H = bottom + 2 * side + mid_gap + top
    fig = plt.figure(figsize=(W, H))
    fitted = long[long.status == FITTED]
    for k, (est, se, sig_col, ref, lab) in enumerate(METRICS):
        y0 = bottom + (1 - k) * (side + mid_gap)
        for j, a in enumerate(H_ANNOTS):
            ax = axes_in(fig, names_w + j * (panel_w + gap), y0, panel_w, side)
            forest(ax, fitted[fitted.annot == a], order, est, se, sig_col, ref, diseases,
                   labels=j == 0)
            ax.set_xlabel(lab)
            ax.set_title(H_ANNOT_LABEL[a], fontsize=PT["title"], pad=3)
        panel_label(fig, 0.02, y0 + side + 0.36, "AB"[k],
                    "Enrichment (95% CI)" if k == 0 else "τ* (95% CI)")
    handles = [Line2D([], [], color=COND_COLORS[c], marker="o", markersize=3, lw=0.6,
                      label=COND_LABEL[c]) for c in CONDS]
    handles += [Line2D([], [], ls="", marker="o", markersize=3, color="#52514e",
                       label=f"FDR < {FDR}"),
                Line2D([], [], ls="", marker="o", markersize=3, color="#52514e",
                       markerfacecolor="white", mew=0.6, label=f"FDR ≥ {FDR}")]
    fig.legend(handles=handles, loc="lower left", ncol=5, borderaxespad=0,
               bbox_to_anchor=(names_w / W, 0.04 / H), handlelength=1.4,
               columnspacing=1.2)
    return fig


def main():
    ap = argparse.ArgumentParser()
    add_args(ap)
    args = resolve(ap.parse_args())
    long = load_long(args.long)
    diseases = fig4.with_ld_flag(load_diseases(), long)
    with print_style():
        fig = draw_figure(long, diseases)
    save(fig, "fig_s_heritability_forest", args.out_dir, exact=True)


if __name__ == "__main__":
    main()
