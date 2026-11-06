#!/usr/bin/env python
"""Supplementary figure: effect of the >= 75-DEG knockdown filter (no new analysis).

  A  knockdowns with >= 1 edge and edges per condition, stacked: kept by deg75
     (solid) and added without the filter (hatched); edges as used for matching
     (|z| >= 2, |logFC| >= 0.5, no self-edges)
  B  real regulator / target set size, deg75 (x) vs nodeg (y); target points
     more than 2x off the diagonal are labelled, and the anchors whose added
     knockdowns cause them are printed
  C  pseudo-anchor recovery of the real set (median, 95% range), deg75 vs nodeg
B and C: one point per trait x condition; colour = condition, filled = IMD,
open = non-IMD control. Edges are rebuilt with overlap_nodeg.knockdown_edges;
B and C read overlap_summary.csv.

Compares the two designs, so it does not follow CH2_VARIANT. Inputs are not
under results/ (left in place): data/DE_stats.suppl_table.csv,
data/DE_results_sig.csv, data/coloc_dat.tsv and
data/permutations/overlap_nodeg/overlap_summary.csv (from
permutation_analysis/scripts/overlap_nodeg.py). The per-design permutation
results are compared in table_s_deg_filter.py.

Writes fig_s_deg_filter.{pdf,png} to figures/comparison/ and prints the panel A
numbers, the median out-degree of deg75 vs added knockdowns, and the labelled
panel B points with their causal anchors.

  uv run python permutation/fig_s_deg_filter.py      (from vis/)
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parents[1]), str(HERE.parents[2] / "permutation_analysis" / "scripts")]
from common import (COND_COLORS, COND_LABEL, COND_SHORT, CONDS, FIGURES, GROUP_SHORT,  # noqa: E402
                    PRINT_W, PT, _plain_number, axes_in, panel_label, plain_log_axis,
                    print_style, save)
from overlap_nodeg import (DATA, MIN_DEGS, OUT, REPO, condition_edges,  # noqa: E402
                           knockdown_edges, load_coloc_genes)

SETS = ["regulators", "targets"]
SET_LABEL = {"regulators": "Regulators", "targets": "Targets"}
GROUP_FILLED = {"IMD": True, "non-IMD": False}  # filled vs open marker
KINDS = {"deg75": "≥ 75 DEGs", "added": "Added without filter"}
COMPARISON_DIR = FIGURES / "comparison"
OFF_DIAGONAL = 2.0  # label target points with nodeg / deg75 size beyond this factor


def load_edges():
    """nodeg edges per condition (deg75 edges are the subset whose knockdown is
    kept by deg75), each tagged with its knockdown's kind; and Ensembl -> symbol."""
    meta = pl.read_csv(DATA / "DE_stats.suppl_table.csv", infer_schema_length=None).select(
        pl.col("target_contrast").alias("cis_ensg"), "culture_condition",
        "n_total_de_genes", "n_cells_target", "ontarget_significant")
    sig = pl.read_csv(DATA / "DE_results_sig.csv")
    raw = {v: knockdown_edges(sig, meta, k) for v, k in MIN_DEGS.items()}
    edges = []
    for c in CONDS:
        kept = condition_edges(raw["deg75"], c, 0.5)["cis_ensg"].unique()
        edges.append(condition_edges(raw["nodeg"], c, 0.5).with_columns(
            pl.lit(c).alias("condition"),
            pl.when(pl.col("cis_ensg").is_in(kept.implode())).then(pl.lit("deg75"))
              .otherwise(pl.lit("added")).alias("kind")))
    symbols = dict(sig.select("cis_ensg", "cis_gene").unique("cis_ensg").iter_rows())
    return pl.concat(edges), symbols


def bars(axes, counts: pl.DataFrame) -> None:
    """Stacked bars per condition: kept (solid) under added (hatched)."""
    for ax, col, title in zip(axes, ["knockdowns", "edges"], ["Knockdowns", "Edges"]):
        for i, c in enumerate(CONDS):
            row = dict(counts.filter(pl.col("condition") == c).select("kind", col).iter_rows())
            ax.bar(i, row["deg75"], width=0.7, color=COND_COLORS[c], lw=0)
            ax.bar(i, row["added"], bottom=row["deg75"], width=0.7, facecolor="white",
                   edgecolor=COND_COLORS[c], hatch="//////", lw=0.6)
        ax.set_xticks(range(len(CONDS)), [COND_SHORT[c] for c in CONDS])
        ax.yaxis.set_major_formatter(_plain_number)
        ax.set_title(title, fontsize=PT["title"])


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


def off_diagonal(w: pl.DataFrame, edges: pl.DataFrame, symbols: dict) -> pl.DataFrame:
    """Target points > OFF_DIAGONAL x off the diagonal, with the anchors whose
    added knockdowns contribute the targets that only nodeg has."""
    traits = pl.read_csv(REPO / "vis/diseases.csv").select("study", "abbrev")
    flagged = (w.filter(pl.col("set") == "targets", pl.col("real_size") > 0)
                .with_columns((pl.col("real_size_nodeg") / pl.col("real_size")).alias("ratio"))
                .filter((pl.col("ratio") > OFF_DIAGONAL) | (pl.col("ratio") < 1 / OFF_DIAGONAL))
                .join(traits, on="study"))
    rows = []
    for r in flagged.iter_rows(named=True):
        anchors = load_coloc_genes(str(DATA / "coloc_dat.tsv"), "CD4", 0.8, r["study"])
        e = edges.filter(pl.col("condition") == r["condition"], pl.col("cis_ensg").is_in(anchors))
        new = set(e["trans_ensg"]) - set(e.filter(pl.col("kind") == "deg75")["trans_ensg"])
        for g, sub in e.filter(pl.col("kind") == "added").group_by("cis_ensg"):
            rows.append(r | {"anchor": g[0], "symbol": symbols.get(g[0]), "anchor_out_degree": sub.height,
                             "new_targets_from_anchor": len(new & set(sub["trans_ensg"])),
                             "new_targets_total": len(new)})
    return pl.DataFrame(rows)


def main() -> None:
    edges, symbols = load_edges()
    w = wide(pl.read_csv(OUT / "overlap_summary.csv"))
    counts = (edges.group_by("condition", "kind")
                   .agg(pl.col("cis_ensg").n_unique().alias("knockdowns"), pl.len().alias("edges")))
    deg = edges.group_by("condition", "kind", "cis_ensg").len("out_degree")
    flags = off_diagonal(w, edges, symbols)

    with print_style():
        fig = plt.figure(figsize=(PRINT_W, 4.75))
        # A: knockdowns and edges, kept vs added
        bars([axes_in(fig, 0.45, 2.95, 0.75, 1.2), axes_in(fig, 1.85, 2.95, 0.75, 1.2)], counts)
        panel_label(fig, 0.05, 4.68, "A", "Knockdowns and edges")
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
            if s == "targets" and flags.height:
                for r in flags.unique(["study", "condition"]).iter_rows(named=True):
                    ax.annotate(f"{r['abbrev']} {COND_SHORT[r['condition']]}",
                                (r["real_size"], r["real_size_nodeg"]), xytext=(0, 7),
                                textcoords="offset points", ha="center", fontsize=PT["tick"])
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
                   + [Patch(facecolor="#52514e", lw=0, label=f"A: {KINDS['deg75']}"),
                      Patch(facecolor="white", edgecolor="#52514e", hatch="//////", lw=0.6,
                            label=f"A: {KINDS['added'].lower()}")]
                   + [Line2D([], [], ls="", marker="o", ms=4, mec="#52514e", mew=0.7,
                             mfc="#52514e" if GROUP_FILLED[g] else "white", label=GROUP_SHORT[g])
                      for g in GROUP_FILLED],
                   loc="lower left", bbox_to_anchor=(4.65 / PRINT_W, 0.5 / 4.75))
        save(fig, "fig_s_deg_filter", COMPARISON_DIR, exact=True)

    order = pl.col("condition").replace_strict({c: i for i, c in enumerate(CONDS)})
    with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=200):
        print("Panel A: knockdowns with >= 1 edge and edges (|logFC| >= 0.5)")
        print(counts.pivot(on="kind", index="condition", values=["knockdowns", "edges"]).sort(order))
        print("Median out-degree per knockdown")
        print(deg.group_by("condition", "kind").agg(pl.len().alias("knockdowns"),
                                                     pl.col("out_degree").median().alias("median_out_degree"))
              .pivot(on="kind", index="condition", values=["knockdowns", "median_out_degree"])
              .sort(order))
        print(f"Panel B: target points > {OFF_DIAGONAL:g}x off the diagonal, by anchor with an added knockdown")
        print(flags.select("abbrev", "condition", "real_size", "real_size_nodeg", "anchor", "symbol",
                           "anchor_out_degree", "new_targets_from_anchor", "new_targets_total")
              if flags.height else "none")
    zero = w.filter((pl.col("real_size") == 0) | (pl.col("real_size_nodeg") == 0))
    print(f"Panel B omits {zero.height} trait x condition x set rows with an empty real set")


if __name__ == "__main__":
    main()
