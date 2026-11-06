#!/usr/bin/env python
"""Jaccard overlap of gene sets between diseases (Fig. 2.1C, 2.6A, S2).

Outputs (in build/<VARIANT>/):
  * jaccard_anchors - coloc genes; one set per disease, condition-independent
    (drawn into Fig. 2.1C by 07_fig_2_1.py via draw_anchor_jaccard)
  * jaccard_{regulators|targets}_{cond|all}[_excess] - one condition
    (--condition) or all three side by side (default)
  * jaccard_{regulators|targets}_{cond}_pair - observed and excess side by
    side for one condition (--condition X --pair; Fig. 2.6A uses Stim8hr)

--excess plots observed - expected Jaccard, where the expectation for two
random sets of sizes a and b drawn from a universe of N genes is
    E = (ab/N) / (a + b - ab/N).
N is the regulator universe (knockdowns with >= 1 edge) for regulator sets
and the target universe (genes with >= 1 incoming edge) for target sets,
per condition (common.universe). Pairs involving an empty set have no
expectation and are drawn grey.

Ordering is average-linkage hierarchical clustering on 1 - observed Jaccard.
By default (--order shared) all conditions share one order, clustered on the
mean observed Jaccard across conditions, so cells line up between panels and
between the observed and excess versions; --order panel clusters each panel
separately and --order anchors reuses the anchor clustering.

Axis labels are disease display names coloured by IMD / non-IMD, with a
class strip on both axes. Observed and expected matrices are written to CSV.
"""

import argparse
import itertools
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

from matplotlib.patches import Rectangle

from common import (ANNOT_LABEL, COND_LABEL, CONDS, GROUP_COLORS,
                    GROUP_LABEL, GROUP_SHORT, GROUPS, PRINT_W, PT, add_variant_args,
                    axes_in, color_ticklabels, load_diseases, load_edges,
                    load_gene_sets, name, panel_label, print_style,
                    resolve_variant, save, universe)

CMAP = "Blues"
CMAP_EXCESS = "RdBu_r"
CMAP_EXCESS_SEQ = "Reds"   # print figure: excess from 0 to max
MASKED = "#f1f1ef"         # masked diagonal in print heatmaps
OUTLINE = "#0b0b0b"
UNDEFINED = "#d6d5d0"
LABEL_PT = 8
CELL_IN = 0.30          # heatmap cell size in inches
LABEL_IN = 2.0          # room for the longest disease name at LABEL_PT


def jaccard_matrix(sets):
    n = len(sets)
    m = np.eye(n)
    for i, j in itertools.combinations(range(n), 2):
        union = len(sets[i] | sets[j])
        m[i, j] = m[j, i] = len(sets[i] & sets[j]) / union if union else 0.0
    return m


def expected_jaccard(sizes, n_universe):
    a = np.asarray(sizes, float)[:, None]
    b = np.asarray(sizes, float)[None, :]
    inter = a * b / n_universe
    with np.errstate(invalid="ignore", divide="ignore"):
        e = inter / (a + b - inter)
    e[(a == 0) | (b == 0)] = np.nan
    np.fill_diagonal(e, 1.0)
    return e


def cluster_order(m):
    d = 1 - m
    np.fill_diagonal(d, 0)
    return leaves_list(linkage(squareform(d, checks=False), "average",
                               optimal_ordering=True))


def offdiag(m):
    return m[~np.eye(len(m), dtype=bool)]


def draw(ax, m, studies, diseases, title, vmin, vmax, cmap,
         show_ylabels=True, annotate=False):
    n = len(studies)
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad(UNDEFINED)
    shown = m.copy()
    np.fill_diagonal(shown, np.nan if vmin < 0 else 1.0)
    im = ax.imshow(np.ma.masked_invalid(shown), cmap=cm, vmin=vmin, vmax=vmax,
                   interpolation="nearest")
    if vmin < 0:  # diagonal is meaningless for excess; draw it neutral
        ax.imshow(np.where(np.eye(n, dtype=bool), 1.0, np.nan),
                  cmap=ListedColormap(["#52514e"]), vmin=0, vmax=1,
                  interpolation="nearest")
    names = [name(diseases, s) for s in studies]
    ax.set_xticks(range(n))
    ax.set_xticklabels(names, rotation=90, fontsize=LABEL_PT)
    color_ticklabels(ax.get_xticklabels(), studies, diseases)
    ax.set_yticks(range(n))
    if show_ylabels:
        ax.set_yticklabels(names, fontsize=LABEL_PT)
        color_ticklabels(ax.get_yticklabels(), studies, diseases)
    else:
        ax.set_yticklabels([])
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)

    # IMD / non-IMD strip along the top and left edges.
    codes = np.array([[GROUPS.index(diseases.loc[s, "group"]) for s in studies]])
    strip = ListedColormap([GROUP_COLORS[g] for g in GROUPS])
    ax.imshow(codes, cmap=strip, vmin=0, vmax=len(GROUPS) - 1, aspect="auto",
              extent=(-0.5, n - 0.5, -0.9, -0.6))
    ax.imshow(codes.T, cmap=strip, vmin=0, vmax=len(GROUPS) - 1, aspect="auto",
              extent=(-0.9, -0.6, n - 0.5, -0.5))
    ax.set_xlim(-0.95, n - 0.5)
    ax.set_ylim(n - 0.5, -0.95)

    if annotate:
        span = max(abs(vmin), vmax)
        for i, j in itertools.product(range(n), range(n)):
            v = m[i, j]
            if i != j and np.isfinite(v) and abs(v) >= 0.05:
                ax.text(j, i, f"{v:.2f}".replace("0.", "."), ha="center", va="center",
                        fontsize=6.5, color="white" if abs(v) > 0.6 * span else "#0b0b0b")
    ax.set_title(title, loc="left", fontsize=11)
    return im


def group_handles(undefined=False, short=False):
    labels = GROUP_SHORT if short else GROUP_LABEL
    handles = [Patch(color=GROUP_COLORS[g], label=labels[g]) for g in GROUPS]
    if undefined:
        handles.append(Patch(color=UNDEFINED, label="Undefined (empty set)"))
    return handles


def legend(fig, axes, undefined=False):
    """Group legend just above the heatmaps (not the figure top)."""
    handles = group_handles(undefined)
    top = max(ax.get_position().y1 for ax in axes)
    fig.legend(handles=handles, loc="lower center", ncol=len(handles),
               bbox_to_anchor=(0.5, top + 0.06))


def write_matrix(m, studies, diseases, path):
    lab = [name(diseases, s) for s in studies]
    pd.DataFrame(m, index=lab, columns=lab).to_csv(path, float_format="%.5g")


def anchor_jaccard(sets, studies):
    """(matrix, clustering order) for the condition-independent anchor sets."""
    m = jaccard_matrix([sets[(s, "Rest", "H4_genes")] for s in studies])
    return m, cluster_order(m)


def draw_anchor_jaccard(ax, sets, diseases, title=None):
    """Draw the anchor Jaccard heatmap onto ``ax``; returns the image."""
    studies = [s for s in diseases.index if (s, "Rest", "H4_genes") in sets]
    m, o = anchor_jaccard(sets, studies)
    return draw(ax, m[np.ix_(o, o)], [studies[i] for i in o], diseases,
                ANNOT_LABEL["H4_genes"] if title is None else title, 0, max(offdiag(m).max(), 0.05),
                CMAP, annotate=True)


def same_gwas_pairs(studies, diseases):
    """Index pairs (i, j), i != j, of traits from the same GWAS (shared samples)."""
    grp = [diseases.loc[s, "sample_group"] for s in studies]
    return [(i, j) for i, j in itertools.permutations(range(len(studies)), 2)
            if isinstance(grp[i], str) and grp[i] == grp[j]]


def draw_print(ax, m, studies, diseases, vmin, vmax, cmap, ylabels=True,
               outline=None, title=None):
    """Print-size heatmap: rows 'Name  k', columns numbered k, diagonal masked.

    ``outline`` is a list of (row, col) index pairs to box. Values outside
    [vmin, vmax] are clipped to the colour range.
    """
    n = len(studies)
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad(MASKED)
    shown = m.astype(float).copy()
    np.fill_diagonal(shown, np.nan)
    im = ax.imshow(np.ma.masked_invalid(shown), cmap=cm, vmin=vmin, vmax=vmax,
                   interpolation="nearest", aspect="auto")
    ax.set_xticks(range(n))
    # Numbers at the tick size; rotated so two digits fit a narrow cell.
    ax.set_xticklabels([str(k + 1) for k in range(n)], rotation=90)
    ax.set_yticks(range(n))
    if ylabels:
        ax.set_yticklabels([f"{name(diseases, s)}  {k + 1}" for k, s in enumerate(studies)])
        color_ticklabels(ax.get_yticklabels(), studies, diseases)
    else:
        ax.set_yticklabels([])
    color_ticklabels(ax.get_xticklabels(), studies, diseases)
    ax.tick_params(length=0, pad=1.5)
    for sp in ax.spines.values():
        sp.set_visible(False)
    codes = np.array([[GROUPS.index(diseases.loc[s, "group"]) for s in studies]])
    strip = ListedColormap([GROUP_COLORS[g] for g in GROUPS])
    ax.imshow(codes, cmap=strip, vmin=0, vmax=len(GROUPS) - 1, aspect="auto",
              extent=(-0.5, n - 0.5, -0.85, -0.6))
    ax.imshow(codes.T, cmap=strip, vmin=0, vmax=len(GROUPS) - 1, aspect="auto",
              extent=(-0.85, -0.6, n - 0.5, -0.5))
    ax.set_xlim(-0.9, n - 0.5)
    ax.set_ylim(n - 0.5, -0.9)
    for i, j in outline or []:
        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec=OUTLINE,
                               lw=0.8, zorder=5))
    if title:
        ax.set_title(title, loc="left", pad=3)
    return im


def draw_anchor_jaccard_print(ax, cax, sets, diseases):
    """Anchor Jaccard heatmap for Fig. 2.1C, with a colour bar in ``cax``."""
    studies = [s for s in diseases.index if (s, "Rest", "H4_genes") in sets]
    m, o = anchor_jaccard(sets, studies)
    im = draw_print(ax, m[np.ix_(o, o)], [studies[i] for i in o], diseases,
                    0, max(offdiag(m).max(), 0.05), CMAP)
    cb = ax.figure.colorbar(im, cax=cax)
    cb.set_label("Jaccard index")
    cb.outline.set_linewidth(0.5)
    return im


def pair_figure(obs, exc, studies, diseases, annot, cond):
    """Fig. 2.6A at print size: observed | observed - expected, shared order."""
    W, H = PRINT_W, 3.6
    n = len(studies)
    fig = plt.figure(figsize=(W, H))
    side, x0, gap, y0 = 2.14, 1.78, 0.2, 0.8
    ax1 = axes_in(fig, x0, y0, side, side)
    ax2 = axes_in(fig, x0 + side + gap, y0, side, side)
    pairs = same_gwas_pairs(studies, diseases)
    obs_max = offdiag(obs).max()
    exc_max = np.nanmax(offdiag(exc))
    im1 = draw_print(ax1, obs, studies, diseases, 0, obs_max, CMAP, outline=pairs,
                     title="Observed Jaccard")
    im2 = draw_print(ax2, exc, studies, diseases, 0, exc_max, CMAP_EXCESS_SEQ,
                     ylabels=False, outline=pairs, title="Observed − expected Jaccard")
    for im, x, label in [(im1, x0, "Jaccard index"),
                         (im2, x0 + side + gap, "Observed − expected")]:
        cax = axes_in(fig, x, 0.34, side, 0.07)
        extend = "min" if im is im2 and np.nanmin(offdiag(exc)) < 0 else "neither"
        cb = fig.colorbar(im, cax=cax, orientation="horizontal", extend=extend)
        cb.set_label(label, labelpad=1)
        cb.outline.set_linewidth(0.5)
    handles = group_handles() + [
        Patch(facecolor="white", edgecolor=OUTLINE, lw=0.8,
              label="Same GWAS (shared samples)")]
    fig.legend(handles=handles, loc="upper left", ncol=3, bbox_to_anchor=(0.0, 1.0),
               handlelength=1.2, columnspacing=1.2, borderaxespad=0.2)
    panel_label(fig, 0.02, H - 0.3, "A", f"{ANNOT_LABEL[annot]}, {COND_LABEL[cond]}")
    return fig


def heatmap_size(n_panels, n):
    side = CELL_IN * n
    return (LABEL_IN + n_panels * (side + 0.25) + 1.2, side + LABEL_IN + 0.8)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--order", choices=["shared", "panel", "anchors"], default="shared",
                    help="'shared': one clustering per annotation from the mean "
                         "observed Jaccard across conditions; 'panel': cluster "
                         "each heatmap separately; 'anchors': use the anchor order.")
    ap.add_argument("--condition", choices=CONDS, default=None,
                    help="Render one condition only (main figure: Stim8hr). "
                         "Default: all three side by side.")
    ap.add_argument("--excess", action="store_true",
                    help="Plot observed - expected Jaccard (diverging colour scale).")
    ap.add_argument("--pair", action="store_true",
                    help="With --condition: observed and excess side by side.")
    ap.add_argument("--annotate", action="store_true",
                    help="Print values in regulator/target cells too.")
    add_variant_args(ap)
    args = resolve_variant(ap.parse_args())
    if args.pair and not args.condition:
        ap.error("--pair needs --condition")
    out_dir = Path(args.out_dir)

    diseases = load_diseases()
    sets = load_gene_sets(args.gene_sets_dir)
    edges = load_edges(min_abs_logfc=args.min_abs_logfc)
    studies = [s for s in diseases.index if (s, "Rest", "H4_genes") in sets]
    n = len(studies)

    # Anchors: one set per disease (identical across conditions), observed only.
    m, anchor_order = anchor_jaccard(sets, studies)
    o = anchor_order
    write_matrix(m[np.ix_(o, o)], [studies[i] for i in o], diseases,
                 out_dir / "jaccard_anchors_observed.csv")
    fig, ax = plt.subplots(figsize=heatmap_size(1, n))
    im = draw_anchor_jaccard(ax, sets, diseases)
    fig.colorbar(im, ax=ax, shrink=0.6, label="Jaccard index")
    legend(fig, [ax])
    save(fig, "jaccard_anchors", out_dir)

    conds = [args.condition] if args.condition else CONDS
    for annot in ["regulators_H4", "targets_H4"]:
        short = annot.split("_")[0]
        obs = {c: jaccard_matrix([sets[(s, c, annot)] for s in studies]) for c in CONDS}
        exp = {}
        for c in CONDS:
            regs, tgts = universe(edges, c)
            n_univ = len(regs) if annot == "regulators_H4" else len(tgts)
            exp[c] = expected_jaccard([len(sets[(s, c, annot)]) for s in studies], n_univ)
        shared = (anchor_order if args.order == "anchors"
                  else cluster_order(np.mean(list(obs.values()), axis=0)))
        per_panel = args.order == "panel"

        for c in CONDS:
            for kind, mm in [("observed", obs[c]), ("expected", exp[c])]:
                write_matrix(mm[np.ix_(shared, shared)], [studies[i] for i in shared],
                             diseases, out_dir / f"jaccard_{short}_{c}_{kind}.csv")

        excess = {c: obs[c] - exp[c] for c in CONDS}
        span = max(np.nanmax(np.abs(offdiag(excess[c]))) for c in conds)
        obs_max = max(offdiag(obs[c]).max() for c in conds)
        observed_style = (0, obs_max, CMAP, "Jaccard index")
        excess_style = (-span, span, CMAP_EXCESS, "Observed − expected Jaccard")

        if args.pair:
            c = args.condition
            o = shared
            with print_style():
                fig = pair_figure(obs[c][np.ix_(o, o)], excess[c][np.ix_(o, o)],
                                  [studies[i] for i in o], diseases, annot, c)
            suffix = f"_{c}_pair" + ("" if args.order == "shared" else f"_{args.order}_order")
            save(fig, f"jaccard_{short}{suffix}", out_dir, exact=True)
            continue
        else:
            style = excess_style if args.excess else observed_style
            src = excess if args.excess else obs
            panels = [(COND_LABEL[c], src[c], style) for c in conds]

        w, h = heatmap_size(len(panels), n)
        if per_panel:
            w += LABEL_IN * (len(panels) - 1)
        fig, axes = plt.subplots(1, len(panels), figsize=(w, h), squeeze=False,
                                 gridspec_kw={"wspace": 0.75 if per_panel else 0.08})
        axes = axes[0]
        for k, (ax, (sub, mm, (vmin, vmax, cmap, cbar))) in enumerate(zip(axes, panels)):
            c = args.condition or conds[k]
            o = cluster_order(obs[c]) if per_panel else shared
            title = (f"{ANNOT_LABEL[annot]} – {COND_LABEL[c]}: {sub}" if args.pair
                     else f"{ANNOT_LABEL[annot]} – {sub}")
            im = draw(ax, mm[np.ix_(o, o)], [studies[i] for i in o], diseases,
                      title, vmin, vmax, cmap,
                      show_ylabels=(k == 0 or per_panel), annotate=args.annotate)
            if args.pair:
                fig.colorbar(im, ax=ax, shrink=0.55, label=cbar, pad=0.02)
        if not args.pair:
            fig.colorbar(im, ax=list(axes), shrink=0.6, label=panels[0][2][3], pad=0.02)
        undefined = any(np.isnan(p[1]).any() for p in panels if p[2][0] < 0)
        legend(fig, axes, undefined=undefined)
        if args.pair:
            suffix = f"_{args.condition}_pair"
        else:
            suffix = f"_{args.condition or 'all'}" + ("_excess" if args.excess else "")
        suffix += "" if args.order == "shared" else f"_{args.order}_order"
        save(fig, f"jaccard_{short}{suffix}", out_dir)


if __name__ == "__main__":
    main()
