#!/usr/bin/env python
"""Supplementary scatter plots, one point per fitted gene set.

  fig_s_enrichment_vs_coverage: log2 enrichment against % reference SNPs in
      the annotation (log axis), colour = annotation, shape = IMD / control.
      Sets with enrichment <= 0 have no log2 value and are left out (listed
      in the caption).
  fig_s_enrichment_vs_tau: signed enrichment z against tau* z. Dashed lines
      at the FDR 5% thresholds of the expanded-set family (regulators +
      targets): the largest p with BH q < 0.05, or the rank-1 BH hurdle
      0.05 / m when no set passes. Enrichment is two-sided (lines at +/- z);
      tau* is one-sided (H1: tau* > 0), so its line is at +z only. Dotted:
      nominal p = 0.05. Anchors have their own BH family (not drawn).
"""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.stats import norm, t as tdist

from hcommon import (FDR, FITTED, H_ANNOT_LABEL, H_ANNOTS, add_args, load_long,
                     resolve)
from common import (ANNOT_COLORS, GROUP_SHORT, GROUPS, PT, axes_in, plain_log_axis,
                    print_style, save)

MARKER = {"IMD": "o", "non-IMD": "^"}
W, H = 4.45, 3.0


def legend_handles():
    hs = [Line2D([], [], ls="", marker="o", markersize=4, color=ANNOT_COLORS[a],
                 label=H_ANNOT_LABEL[a]) for a in H_ANNOTS]
    hs += [Line2D([], [], ls="", marker=MARKER[g], markersize=4, color="#52514e",
                  markerfacecolor="none", mew=0.7, label=GROUP_SHORT[g]) for g in GROUPS]
    return hs


def points(ax, df, x, y):
    for a in H_ANNOTS:
        for g in GROUPS:
            d = df[(df.annot == a) & (df.group == g)]
            ax.scatter(d[x], d[y], s=11, marker=MARKER[g], color=ANNOT_COLORS[a],
                       alpha=0.85, lw=0, zorder=3)


def coverage_figure(fit):
    df = fit[fit.enrichment > 0].assign(pct=lambda d: 100 * d.prop_snps)
    fig = plt.figure(figsize=(W, H))
    ax = axes_in(fig, 0.5, 0.42, 2.45, 2.45)
    ax.axhline(0, color="#8a8985", lw=0.5, ls="--", zorder=1)
    points(ax, df, "pct", "log2_enrichment")
    ax.set_xscale("log")
    ax.set_xlim(df.pct.min() / 1.6, df.pct.max() * 1.6)
    plain_log_axis(ax, "x")
    ax.set_xlabel("Reference SNPs in annotation (%)")
    ax.set_ylabel("log$_2$ enrichment")
    fig.legend(handles=legend_handles(), loc="upper left", bbox_to_anchor=(3.05 / W, 2.87 / H),
               borderaxespad=0, handletextpad=0.3)
    return fig, fit[fit.enrichment <= 0]


def fdr_p(fam, p_col, q_col):
    """p threshold for BH FDR 5% within one family."""
    sig = fam.loc[fam[q_col] < FDR, p_col]
    return sig.max() if len(sig) else FDR / fam[p_col].notna().sum()


def tau_figure(fit):
    fig = plt.figure(figsize=(W, H))
    ax = axes_in(fig, 0.5, 0.42, 2.45, 2.45)
    # Enrichment z comes from ldsc's two-sided t test (200 jackknife blocks);
    # tau* z is normal and tested one-sided.
    exp = fit[fit.fdr_family == "expanded"]
    ze = tdist(200).isf(fdr_p(exp, "enrichment_p", "enrichment_fdr") / 2)
    zt = norm.isf(fdr_p(exp, "tau_star_p", "tau_star_fdr"))
    for z, ls in [(1.96, ":"), (ze, "--")]:
        for s in (-1, 1):
            ax.axvline(s * z, color="#8a8985", lw=0.5, ls=ls, zorder=1)
    for z, ls in [(norm.isf(0.05), ":"), (zt, "--")]:
        ax.axhline(z, color="#8a8985", lw=0.5, ls=ls, zorder=1)
    points(ax, fit, "enrichment_z", "tau_star_z")
    ax.set_xlabel("Enrichment z")
    ax.set_ylabel("τ* z")
    hs = legend_handles() + [
        Line2D([], [], color="#8a8985", lw=0.6, ls="--", label=f"FDR {FDR} (expanded)"),
        Line2D([], [], color="#8a8985", lw=0.6, ls=":", label="p = 0.05")]
    fig.legend(handles=hs, loc="upper left", bbox_to_anchor=(3.05 / W, 2.87 / H),
               borderaxespad=0, handletextpad=0.3)
    return fig, ze, zt


def main():
    ap = argparse.ArgumentParser()
    add_args(ap)
    args = resolve(ap.parse_args())
    fit = load_long(args.long).query("status == @FITTED")
    with print_style():
        fig, dropped = coverage_figure(fit)
        save(fig, "fig_s_enrichment_vs_coverage", args.out_dir, exact=True)
        fig, ze, zt = tau_figure(fit)
        save(fig, "fig_s_enrichment_vs_tau", args.out_dir, exact=True)
    print("Not on coverage plot (enrichment <= 0):",
          dropped[["disease", "condition", "annot", "enrichment"]].to_dict("records"))
    print(f"FDR thresholds: |enrichment z| >= {ze:.2f}, |tau* z| >= {zt:.2f}")


if __name__ == "__main__":
    main()
