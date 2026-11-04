#!/usr/bin/env python
"""Supplementary figure: effect of the >= 75-DEG knockdown filter (no new analysis).

  A  out-degree CCDF per condition: knockdowns kept by deg75 (solid) vs the
     knockdowns only nodeg adds (dashed); edges as used for matching
     (|z| >= 2, |logFC| >= 0.5, no self-edges)
  B  real regulator / target set size, deg75 (x) vs nodeg (y)
  C  pseudo-anchor recovery of the real set (median, 95% range), deg75 vs nodeg
B and C: one point per trait x condition; colour = condition, filled = IMD,
open = non-IMD control. Edges are rebuilt with overlap_nodeg.knockdown_edges;
B and C read overlap_summary.csv.

Writes fig_s_deg_filter.{pdf,png} to data/permutations/overlap_nodeg/ and prints
the median out-degree of deg75 vs added knockdowns per condition.

  uv run python figures/fig_s_deg_filter.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parents[2] / "vis"), str(HERE.parents[1] / "scripts")]
from common import (COND_COLORS, COND_LABEL, CONDS, GROUP_SHORT, PRINT_W, PT,  # noqa: E402
                    axes_in, panel_label, plain_log_axis, print_style, save)
from overlap_nodeg import DATA, MIN_DEGS, OUT, condition_edges, knockdown_edges  # noqa: E402

SETS = ["regulators", "targets"]
SET_LABEL = {"regulators": "Regulators", "targets": "Targets"}
GROUP_FILLED = {"IMD": True, "non-IMD": False}  # filled vs open marker


def out_degrees() -> pl.DataFrame:
    """Out-degree per knockdown x condition, tagged kept (deg75) or added (nodeg only)."""
    meta = pl.read_csv(DATA / "DE_stats.suppl_table.csv", infer_schema_length=None).select(
        pl.col("target_contrast").alias("cis_ensg"), "culture_condition",
        "n_total_de_genes", "n_cells_target", "ontarget_significant")
    sig = pl.read_csv(DATA / "DE_results_sig.csv")
    deg = {v: pl.concat([condition_edges(knockdown_edges(sig, meta, k), c, 0.5)
                         .group_by("cis_ensg").len("out_degree").with_columns(pl.lit(c).alias("condition"))
                         for c in CONDS]) for v, k in MIN_DEGS.items()}
    added = deg["nodeg"].join(deg["deg75"], on=["cis_ensg", "condition"], how="anti")
    return pl.concat([deg["deg75"].with_columns(pl.lit("deg75").alias("kind")),
                      added.with_columns(pl.lit("added").alias("kind"))])


def ccdf(ax, deg: pl.DataFrame) -> None:
    for c in CONDS:
        for kind, ls in [("deg75", "-"), ("added", "--")]:
            v = np.sort(deg.filter(pl.col("condition") == c, pl.col("kind") == kind)["out_degree"].to_numpy())
            ax.step(v, 1 - np.arange(v.size) / v.size, where="post", color=COND_COLORS[c], ls=ls)
    ax.set(xscale="log", yscale="log", xlabel="Out-degree (genes changed)",
           ylabel="Fraction of knockdowns ≥ x")
    plain_log_axis(ax, "x")
    plain_log_axis(ax, "y")


def thin_decades(ax, axis) -> None:
    """Label every other decade when four or more decades share a small axis."""
    a = ax.xaxis if axis == "x" else ax.yaxis
    lo, hi = ax.get_xlim() if axis == "x" else ax.get_ylim()
    if np.log10(hi / lo) < 4:
        return
    fmt = a.get_major_formatter()
    a.set_major_formatter(lambda v, pos: fmt(v, pos) if round(np.log10(v)) % 2 == 0 else "")


def scatter(ax, w: pl.DataFrame, x: str, y: str, xerr=None, yerr=None) -> None:
    """deg75 (x) vs nodeg (y), one point per trait x condition, diagonal."""
    for (c, g), sub in w.group_by("condition", "group"):
        col = COND_COLORS[c]
        if xerr:
            ax.errorbar(sub[x], sub[y], fmt="none", ecolor=col, alpha=0.35, lw=0.5, zorder=1,
                        xerr=[sub[x] - sub[xerr[0]], sub[xerr[1]] - sub[x]],
                        yerr=[sub[y] - sub[yerr[0]], sub[yerr[1]] - sub[y]])
        ax.scatter(sub[x], sub[y], s=12, lw=0.7, edgecolor=col, zorder=3,
                   facecolor=col if GROUP_FILLED[g] else "white")


def diagonal(ax, lo, hi) -> None:
    ax.plot([lo, hi], [lo, hi], color="#b5b4ae", lw=0.6, ls=":", zorder=0)
    ax.set(xlim=(lo, hi), ylim=(lo, hi), aspect="equal")


def wide(summary: pl.DataFrame) -> pl.DataFrame:
    keys = ["study", "group", "condition", "set"]
    cols = ["real_size", "recovery_median", "recovery_q025", "recovery_q975"]
    d, n = (summary.filter(pl.col("variant") == v).select(keys + cols) for v in ("deg75", "nodeg"))
    return d.join(n, on=keys, suffix="_nodeg")


def main() -> None:
    deg = out_degrees()
    w = wide(pl.read_csv(OUT / "overlap_summary.csv"))

    with print_style():
        fig = plt.figure(figsize=(PRINT_W, 4.75))
        # A: CCDF
        a = axes_in(fig, 0.62, 2.95, 1.75, 1.45)
        ccdf(a, deg)
        panel_label(fig, 0.05, 4.68, "A", "Out-degree of knockdowns")
        # B: real set sizes (log), zero-size sets cannot be shown
        panel_label(fig, 2.75, 4.68, "B", "Real set size")
        for i, s in enumerate(SETS):
            ax = axes_in(fig, 3.2 + i * 1.65, 2.95, 1.2, 1.2)
            sub = w.filter(pl.col("set") == s, pl.col("real_size") > 0, pl.col("real_size_nodeg") > 0)
            scatter(ax, sub, "real_size", "real_size_nodeg")
            ax.set(xscale="log", yscale="log")
            diagonal(ax, 1, 10 ** np.ceil(np.log10(sub["real_size_nodeg"].max())))
            for axis in "xy":
                plain_log_axis(ax, axis)
                thin_decades(ax, axis)
            ax.set_title(SET_LABEL[s], fontsize=PT["title"])
            ax.set_xlabel("≥ 75 DEGs")
            ax.set_ylabel("No DEG filter" if i == 0 else "")
        # C: recovery
        panel_label(fig, 0.05, 2.45, "C", "Pseudo-anchor recovery of the real set")
        for i, s in enumerate(SETS):
            ax = axes_in(fig, 0.62 + i * 2.05, 0.4, 1.65, 1.65)
            sub = w.filter(pl.col("set") == s).drop_nulls("recovery_median")
            scatter(ax, sub, "recovery_median", "recovery_median_nodeg",
                    ("recovery_q025", "recovery_q975"), ("recovery_q025_nodeg", "recovery_q975_nodeg"))
            diagonal(ax, 0, 1)
            ax.set_title(SET_LABEL[s], fontsize=PT["title"])
            ax.set_xlabel("≥ 75 DEGs")
            ax.set_ylabel("No DEG filter" if i == 0 else "")
        # Shared legend: condition colour, IMD / control fill.
        fig.legend(handles=[Line2D([], [], color=COND_COLORS[c], lw=1.2, label=COND_LABEL[c]) for c in CONDS]
                   + [Line2D([], [], color="#52514e", ls="-", label="A: ≥ 75 DEGs"),
                      Line2D([], [], color="#52514e", ls="--", label="A: added without filter")]
                   + [Line2D([], [], ls="", marker="o", ms=4, mec="#52514e", mew=0.7,
                             mfc="#52514e" if GROUP_FILLED[g] else "white", label=GROUP_SHORT[g])
                      for g in GROUP_FILLED],
                   loc="lower left", bbox_to_anchor=(4.65 / PRINT_W, 0.5 / 4.75))
        save(fig, "fig_s_deg_filter", OUT, exact=True)

    print("Median out-degree per knockdown (|logFC| >= 0.5 edges)")
    with pl.Config(tbl_rows=-1):
        print(deg.group_by("condition", "kind").agg(pl.len().alias("knockdowns"),
                                                     pl.col("out_degree").median().alias("median_out_degree"))
              .pivot(on="kind", index="condition", values=["knockdowns", "median_out_degree"])
              .sort(pl.col("condition").replace_strict({c: i for i, c in enumerate(CONDS)})))
    zero = w.filter((pl.col("real_size") == 0) | (pl.col("real_size_nodeg") == 0))
    print(f"Panel B omits {zero.height} trait x condition x set rows with an empty real set")


if __name__ == "__main__":
    main()
