#!/usr/bin/env python
"""fig_h1_overlap - does excluding hub regulators change cross-disease overlap?

Regulator sets only (targets: fig_s_h1_targets_summary; they change only
through anchor loss).

  A. Regulator-set Jaccard, Stim 8 h (condition of Fig. 2.6A), one heatmap
     per variant in the deg0 trait order; top observed, bottom observed -
     expected (size-only expectation of 04_jaccard_heatmaps.py: random sets of
     the same sizes from the variant's hub-filtered universe). Shared scale.
  B. Median (IQR) size-only excess Jaccard vs threshold by pair class, per
     condition. Related-GWAS pairs (ASTCO/ASTAO, RA EAS/EUR, IBD-CD, IBD-UC)
     are left out of the summaries and drawn as a faint dashed line.
  E. Same with the degree-aware expectation of 01b_degree_null.py (sets
     redrawn with matched out-degree profiles, 1,000 draws); filled =
     empirical one-sided p, BH FDR < 0.05 within each variant over the 9
     class x condition tests.
  C. Spearman correlation of each variant's pairwise size-only excess vector
     with deg0's (pairs defined in both), per condition.
  D. Median regulator-set size as % of the variant's regulator universe.

Writes data/h1b_excess_summary.csv (regulators and targets), h1c_spearman.csv,
h1d_set_size.csv (panel A: jaccard_pairs.csv; E: h1e_degree_null_summary.csv).
"""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import hxcommon  # noqa: F401  (sets sys.path, CH2_VARIANT)
from common import (COND_COLORS, CONDS, PRINT_W, PT, axes_in, panel_label,
                    print_style, save)
from hxcommon import (DOUT, MAIN, MAIN_COND, PAIR_CLASSES, PCT, VKEYS,
                      load_diseases, read)
from hxplot import XPOS, cond_handles, heatmap_block, threshold_axis

ROWS = [("regulators_H4", False, "Regulators"),
        ("targets_H4", False, "Targets, all sets"),
        ("targets_H4", True, "Targets, anchors unchanged")]
ROW_YLABEL = {"Regulators": "Regulators", "Targets, all sets": "Targets\n(all sets)",
              "Targets, anchors unchanged": "Targets\n(anchors\nunchanged)"}
CDODGE = {c: (i - 1) * 0.12 for i, c in enumerate(CONDS)}


def excess_summary(pairs):
    rows = []
    for annot, unchanged, label in ROWS:
        p = pairs[pairs.annot == annot]
        if unchanged:
            p = p[p.both_unchanged_10]
        for (v, c, k), g in p.groupby(["variant", "condition", "pair_class"]):
            x = g.excess.dropna()
            rows.append({"set_type": label, "variant": v, "pct_excluded": PCT[v],
                         "condition": c, "pair_class": k, "n_pairs": len(x),
                         "median": x.median(), "q25": x.quantile(0.25),
                         "q75": x.quantile(0.75)})
    return pd.DataFrame(rows)


def spearman_table(pairs):
    rows = []
    for annot, unchanged, label in ROWS:
        p = pairs[pairs.annot == annot]
        if unchanged:
            p = p[p.both_unchanged_10]
        for c in CONDS:
            key = ["study_a", "study_b"]
            base = p[(p.variant == "deg0") & (p.condition == c)].set_index(key).excess
            for v in VKEYS:
                x = p[(p.variant == v) & (p.condition == c)].set_index(key).excess
                ok = base.notna() & x.reindex(base.index).notna()
                rho = spearmanr(base[ok], x.reindex(base.index)[ok]).statistic
                rows.append({"set_type": label, "condition": c, "variant": v,
                             "pct_excluded": PCT[v], "n_pairs": int(ok.sum()),
                             "spearman_rho": rho})
    return pd.DataFrame(rows)


def size_table(sizes):
    r = sizes[sizes.annot == "regulators_H4"]
    return (r.groupby(["variant", "pct_excluded", "condition"])
            .agg(median_pct_of_universe=("pct_of_universe", "median"),
                 median_n_genes=("n_genes", "median"), universe=("universe", "first"),
                 n_traits=("study", "size"))
            .reset_index().sort_values(["condition", "pct_excluded"]))


def line(ax, x, y, c, lo=None, hi=None, filled=None, **kw):
    """Medians across thresholds for one condition; open markers where
    ``filled`` is False (significance), otherwise filled."""
    xs = np.array([XPOS[v] for v in x]) + CDODGE[c]
    ax.plot(xs, y, color=COND_COLORS[c], lw=0.8, **kw)
    f = np.ones(len(xs), bool) if filled is None else np.asarray(filled, bool)
    ax.scatter(xs[f], np.asarray(y)[f], s=5, color=COND_COLORS[c], zorder=3,
               alpha=kw.get("alpha", 1))
    ax.scatter(xs[~f], np.asarray(y)[~f], s=5, facecolor="white",
               edgecolor=COND_COLORS[c], lw=0.6, zorder=3)
    if lo is not None:
        ax.vlines(xs, lo, hi, color=COND_COLORS[c], lw=0.5)


def class_panels(fig, x0, top, rows, summ_rel, bw=0.98, bh=0.72, sharey=None,
                 xlabel=True):
    """Rows of pair-class panels (columns = PAIR_CLASSES); y shared by all.
    ``rows``: list of (summary table, set_type label, y label). Points are
    medians with IQR bars; filled when the table has ``sig`` True. Related
    pairs: faint dashed line in the IMD-IMD panel. Returns the bottom (in)."""
    first = sharey
    for i, (summ, label, ylab) in enumerate(rows):
        for j, k in enumerate(PAIR_CLASSES):
            ax = axes_in(fig, x0 + j * (bw + 0.1), top - (i + 1) * bh - i * 0.14, bw, bh,
                         sharey=first)
            first = first or ax
            s = summ[(summ.set_type == label) & (summ.pair_class == k)]
            for c in CONDS:
                g = s[s.condition == c].sort_values("pct_excluded")
                line(ax, g.variant, g["median"], c, g.q25, g.q75,
                     filled=g["sig"] if "sig" in g else None)
            if k == "IMD–IMD":
                r = summ_rel[i]
                for c in CONDS:
                    g = r[r.condition == c].sort_values("pct_excluded")
                    line(ax, g.variant, g["median"], c, ls="--", alpha=0.3)
            ax.axhline(0, color="#9a9994", lw=0.5, zorder=0)
            threshold_axis(ax, label=xlabel and (i == len(rows) - 1 and j == 1))
            if j > 0:
                ax.tick_params(axis="y", labelleft=False)
            if i < len(rows) - 1:
                ax.set_xticklabels([])
            if i == 0:
                ax.set_title(k, fontsize=PT["tick"], pad=3)
            if j == 0:
                ax.set_ylabel(ylab, fontsize=PT["tick"])
    return top - len(rows) * (bh + 0.14), first


def degree_summary():
    """fig_h1 E table in the column names of excess_summary."""
    d = read("h1e_degree_null_summary.csv").rename(
        columns={"median_excess_degree": "median"})
    d["set_type"] = "Regulators"
    return d


def main():
    diseases = load_diseases()
    pairs, sizes = read("jaccard_pairs.csv"), read("set_sizes.csv")
    summ = excess_summary(pairs)
    summ.to_csv(DOUT / "h1b_excess_summary.csv", index=False)
    rho = spearman_table(pairs)
    rho.to_csv(DOUT / "h1c_spearman.csv", index=False)
    size = size_table(sizes)
    size.to_csv(DOUT / "h1d_set_size.csv", index=False)
    deg = degree_summary()

    with print_style():
        W, H = PRINT_W, 7.15
        fig = plt_figure(W, H)
        panel_label(fig, 0.02, H - 0.05, "A", "Regulator-set overlap, Stim 8 h")
        side, gap = 1.30, 0.07
        bottom = heatmap_block(fig, 0.42, H - 0.45 - side, side, gap, pairs, sizes,
                               MAIN_COND, "regulators_H4", diseases)

        # B (size-only expectation) and E (degree-aware expectation), regulators
        top = bottom - 0.25
        x0 = 0.7
        panel_label(fig, 0.02, top, "B", "Excess Jaccard by pair class, regulators")
        reg = summ[summ.set_type == "Regulators"]
        rel = summ[summ.pair_class == "related"]
        ytop = top - 0.62
        mid, share = class_panels(fig, x0, ytop,
                                  [(reg, "Regulators", "Size-only\nexpectation")],
                                  [rel[rel.set_type == "Regulators"]], xlabel=False)
        panel_label(fig, 0.02, mid - 0.08, "E", "Degree-aware expectation "
                    "(filled: empirical FDR < 0.05)")
        class_panels(fig, x0, mid - 0.45,
                     [(deg[deg.pair_class != "related"], "Regulators",
                       "Degree-aware\nexpectation")],
                     [deg[deg.pair_class == "related"]], sharey=share)
        fig.legend(handles=cond_handles() + [plt_line()], loc="upper left", ncol=4,
                   bbox_to_anchor=(x0 / W, (top - 0.17) / H), fontsize=PT["tick"])

        # C: Spearman with deg0; D: regulator-set size
        cx, cw, ch = 4.65, 1.45, 0.78
        panel_label(fig, cx - 0.4, top, "C", "Structure retained")
        ax = axes_in(fig, cx, ytop - ch, cw, ch)
        for c in CONDS:
            g = rho[(rho.set_type == "Regulators") & (rho.condition == c)]
            line(ax, g.variant, g.spearman_rho, c)
        ax.set_ylabel("Spearman ρ of\nexcess vs 0 %", fontsize=PT["tick"])
        threshold_axis(ax)
        dy = ytop - ch - 0.8 - ch
        panel_label(fig, cx - 0.4, dy + ch + 0.3, "D", "Regulator-set size")
        ax = axes_in(fig, cx, dy, cw, ch)
        for c in CONDS:
            g = size[size.condition == c].sort_values("pct_excluded")
            line(ax, g.variant, g.median_pct_of_universe, c)
        ax.set_ylabel("Median size\n(% of regulators)", fontsize=PT["tick"])
        ax.set_ylim(bottom=0)
        threshold_axis(ax)
        save(fig, "fig_h1_overlap", MAIN, exact=True)


def plt_figure(w, h):
    import matplotlib.pyplot as plt
    return plt.figure(figsize=(w, h))


def plt_line():
    from matplotlib.lines import Line2D
    return Line2D([], [], color="#6f6e69", ls="--", alpha=0.4, label="Related GWAS")


if __name__ == "__main__":
    main()
