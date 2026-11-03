#!/usr/bin/env python
"""Figure 4 - heritability enrichment of anchor, regulator and target sets.

  A. Dot heatmap of S-LDSC enrichment: colour = log2 enrichment (diverging,
     centred on 0), dot area = -log10 p, black outline = FDR < 0.05.
  B. Same layout for tau* (standardised effect conditional on baselineLD
     v2.2 and the other custom annotations in the joint model); dot area =
     -log10 one-sided p (H1: tau* > 0).

Outlines use the *_sig columns: BH FDR < 0.05 within expanded sets and
within anchors separately; sets with negative Prop._h2 never outlined.

Rows are traits in the order and IMD / control blocks of Fig. 2.2D; columns
are anchors, regulators, targets x Rest, 8 h, 48 h. Grey squares: empty set
(light) or no estimate (with x: run failed, not run, or annotation dropped
by the 1,000-SNP filter). Reads heritability_long.csv from
01_collect_results.py.
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from hcommon import (EMPTY, FDR, FITTED, H_ANNOT_LABEL, H_ANNOTS, add_args,
                     flagged_name, load_diseases, load_long, resolve)
from common import (COND_SHORT, CONDS, PRINT_W, PT, axes_in, color_ticklabels,
                    disease_order, panel_label, print_style, save, shade_groups)

CMAP = plt.get_cmap("RdBu_r")
GREY = "#dcdbd6"
CELL_W, ROW_H, GROUP_GAP = 0.2, 0.165, 0.1   # inches
P_CAP = 10                                    # -log10 p at which dot size saturates
S_MIN, S_MAX = 4.0, 100.0                     # marker area (pt^2)
LIMITS = {"log2_enrichment": 3.0, "tau_star": 0.5}


def column_x():
    """x position (in cell units) of each (annot, cond) column, with gaps."""
    gap = GROUP_GAP / CELL_W
    return {(a, c): i * (len(CONDS) + gap) + j
            for i, a in enumerate(H_ANNOTS) for j, c in enumerate(CONDS)}


def dot_size(p):
    v = np.clip(-np.log10(p), 0, P_CAP) / P_CAP
    return S_MIN + (S_MAX - S_MIN) * v


def trait_order(long, diseases):
    anchors = (long[long.annot == "H4_genes"].groupby("study")["n_genes"].max())
    return disease_order(diseases, anchors)


def with_ld_flag(diseases, long):
    d = diseases.copy()
    d["ld_mismatch"] = long.groupby("study")["ld_mismatch"].first().reindex(d.index)
    return d


def dot_heatmap(ax, long, order, value, p_col, sig_col, norm, diseases, labels=True):
    xs = column_x()
    n = len(order)
    ypos = {s: n - 1 - i for i, s in enumerate(order)}
    shade_groups(ax, order[::-1], diseases, label=False)
    for _, r in long.iterrows():
        x, y = xs[(r.annot, r.condition)], ypos[r.study]
        if r.status != FITTED:
            ax.scatter(x, y, s=48, marker="s", color=GREY, lw=0, zorder=2)
            if r.status != EMPTY:
                ax.scatter(x, y, s=14, marker="x", color="#6f6e69", lw=0.6, zorder=3)
            continue
        v = r[value]
        if value == "log2_enrichment" and not np.isfinite(v):
            v = -np.inf   # enrichment <= 0: plotted at the depleted end
        colour = CMAP(norm(np.clip(v, norm.vmin, norm.vmax)))
        sig = bool(r[sig_col])
        ax.scatter(x, y, s=dot_size(r[p_col]), color=colour, zorder=3,
                   edgecolors="black" if sig else "none", linewidths=0.7 if sig else 0)
    x_max = max(xs.values())
    ax.set_xlim(-0.6, x_max + 0.6)
    ax.set_ylim(-0.5, n - 0.5)
    ax.set_yticks(range(n))
    if labels:
        ax.set_yticklabels([flagged_name(diseases, s) for s in order[::-1]])
        color_ticklabels(ax.get_yticklabels(), order[::-1], diseases)
    else:
        ax.tick_params(axis="y", labelleft=False)
    ax.tick_params(axis="y", length=0)
    ax.set_xticks([xs[k] for k in xs])
    ax.set_xticklabels([COND_SHORT[c] for _, c in xs], rotation=90)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=1.5)
    for s in ax.spines.values():
        s.set_visible(False)
    for t in ax.get_xticklabels():
        t.set_fontsize(PT["tick"] - 0.5)
    # Annotation group headers above the rotated condition labels.
    rise = 18 / 72   # inches above the axes top
    line_y = 1 + rise / (n * ROW_H)
    for a in H_ANNOTS:
        lo, hi = xs[(a, CONDS[0])], xs[(a, CONDS[-1])]
        ax.plot([lo - 0.4, hi + 0.4], [line_y, line_y], color="#52514e", lw=0.6,
                transform=ax.get_xaxis_transform(), clip_on=False)
        ax.text((lo + hi) / 2, line_y + 0.02 / (n * ROW_H), H_ANNOT_LABEL[a],
                transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                fontsize=PT["tick"])


def size_handles():
    out = []
    for p in [0.05, 1e-3, 1e-6, 10.0 ** -P_CAP]:
        lab = "0.05" if p == 0.05 else f"10$^{{{int(np.log10(p))}}}$"
        if p == 10.0 ** -P_CAP:
            lab = "≤ " + lab
        out.append(Line2D([], [], ls="", marker="o", color="#6f6e69",
                          markersize=np.sqrt(dot_size(p)), label=lab))
    return out


def status_handles():
    return [Line2D([], [], ls="", marker="o", markerfacecolor="white",
                   markeredgecolor="black", markeredgewidth=0.7, markersize=6,
                   label=f"FDR < {FDR}"),
            Patch(facecolor=GREY, edgecolor="none", label="Empty set"),
            Line2D([], [], ls="", marker="x", color="#6f6e69", markersize=4,
                   markeredgewidth=0.6, label="No estimate")]


def draw_figure(long, diseases):
    order = trait_order(long, diseases)
    n = len(order)
    xs = column_x()
    panel_w = (max(xs.values()) + 1.2) * CELL_W
    side = n * ROW_H
    names_w, gap = 1.62, 0.22
    W = PRINT_W
    ax_y = 0.92                              # matrix bottom (legends below)
    H = ax_y + side + 0.62
    fig = plt.figure(figsize=(W, H))
    xa = names_w
    xb = xa + panel_w + gap
    ax_a = axes_in(fig, xa, ax_y, panel_w, side)
    ax_b = axes_in(fig, xb, ax_y, panel_w, side)

    norm_e = Normalize(-LIMITS["log2_enrichment"], LIMITS["log2_enrichment"])
    norm_t = Normalize(-LIMITS["tau_star"], LIMITS["tau_star"])
    dot_heatmap(ax_a, long, order, "log2_enrichment", "enrichment_p", "enrichment_sig",
                norm_e, diseases, labels=True)
    dot_heatmap(ax_b, long, order, "tau_star", "tau_star_p", "tau_star_sig",
                norm_t, diseases, labels=False)
    panel_label(fig, 0.02, H - 0.04, "A", "Enrichment")
    panel_label(fig, xb - 0.2, H - 0.04, "B", "τ* (conditional effect)")

    # Colour bars under each panel, shared size / status legend below.
    cb_w = panel_w * 0.8
    for ax_x, norm, lab in [(xa, norm_e, "log$_2$ enrichment"), (xb, norm_t, "τ*")]:
        cax = axes_in(fig, ax_x + (panel_w - cb_w) / 2, ax_y - 0.2, cb_w, 0.07)
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cax,
                          orientation="horizontal", extend="both")
        cb.outline.set_linewidth(0.5)
        cb.ax.tick_params(labelsize=PT["tick"], length=2, width=0.5, pad=1)
        cb.set_label(lab, fontsize=PT["label"], labelpad=1)
    fig.legend(handles=size_handles(), title="p", loc="lower left", ncol=4,
               bbox_to_anchor=(xa / W, 0.02 / H), borderaxespad=0, handletextpad=0.2,
               columnspacing=0.9, title_fontsize=PT["tick"], alignment="left")
    fig.legend(handles=status_handles(), loc="lower left", ncol=3,
               bbox_to_anchor=(xb / W, 0.06 / H), borderaxespad=0, handletextpad=0.3,
               columnspacing=0.9, handlelength=1.0)
    return fig


def main():
    ap = argparse.ArgumentParser()
    add_args(ap)
    args = resolve(ap.parse_args())
    long = load_long(args.long)
    diseases = with_ld_flag(load_diseases(), long)
    with print_style():
        fig = draw_figure(long, diseases)
    save(fig, "fig_4_heritability", args.out_dir, exact=True)


if __name__ == "__main__":
    main()
