#!/usr/bin/env python
"""Plot the S-LDSC permutation null distributions with the real value marked.

Reads permutation_ldsc_perperm.csv (from collect_permutation_ldsc.py) and, for
each condition, draws a grid of histograms: one row per annotation
(anchors = H4_genes / regulators_H4 / targets_H4) and one column per statistic
(enrichment, coefficient_z by default). Each panel shows the matched-anchor
null distribution as a histogram with the real (observed) value as a red
vertical line, annotated with the one-sided empirical p and n.

Empirical p = (#null >= real + 1) / (n_null + 1), matching
collect_permutation_ldsc.py (real should be larger than a matched-random set).

Reads every ``PERM_DIR/<study>_<cond>/permutation_ldsc_perperm.csv`` (PERM_DIR
follows CH2_VARIANT, see common.py) and writes
``FIG_DIR/permutation/<study>/permutation_null_<cond>.png``.

Usage (from Chapter_2/vis):
  uv run heritability/plot_permutation_null.py [--conditions Stim48hr]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless / cluster
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import polars as pl  # noqa: E402

from hcommon import FIG_DIR, PERM_DIR  # noqa: E402

# hcommon -> common sets the chapter rcParams on import; this figure keeps
# matplotlib's defaults, as it had in the pipeline.
matplotlib.rcdefaults()
matplotlib.use("Agg")

STAT_LABELS = {
    "enrichment": "Enrichment",
    "coefficient_z": "Coefficient z-score",
    "coefficient": "Coefficient (tau)",
    "enrichment_p": "Enrichment p",
    "prop_h2": "Prop. h2",
}
SET_LABELS = {"H4_genes": "anchors", "regulators_H4": "regulators",
              "targets_H4": "targets"}


def emp_p(real, null):
    null = np.asarray([x for x in null if x is not None and not np.isnan(x)])
    if real is None or np.isnan(real) or null.size == 0:
        return float("nan"), int(null.size)
    return float((np.sum(null >= real) + 1) / (null.size + 1)), int(null.size)


def panel(ax, real, null, stat, bins):
    null = np.asarray([x for x in null if x is not None and not np.isnan(x)])
    if null.size == 0:
        ax.text(0.5, 0.5, "no null data", ha="center", va="center",
                transform=ax.transAxes, color="grey")
        ax.set_yticks([])
        return
    ax.hist(null, bins=bins, color="#4C72B0", alpha=0.75,
            edgecolor="white", linewidth=0.3)
    p, n = emp_p(real, null)
    if real is not None and not np.isnan(real):
        ax.axvline(real, color="#C44E52", lw=2.2,
                   label=f"real = {real:.2f}")
        ax.legend(fontsize=8, loc="upper right", frameon=False)
        ptxt = "p < 0.001" if p < 1e-3 else f"p = {p:.3f}"
        ax.set_title(f"{ptxt}   (n={n})", fontsize=9)
    else:
        ax.set_title(f"no real run (n={n})", fontsize=9, color="grey")
    ax.set_ylabel("permutations", fontsize=8)


def plot_condition(df_cond, cond, annots, stats, bins, study, out_path):
    real_df = df_cond.filter(pl.col("perm") == "real")
    null_df = df_cond.filter(pl.col("perm") != "real")

    nrows, ncols = len(annots), len(stats)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 2.8 * nrows),
                             squeeze=False)
    for i, annot in enumerate(annots):
        for j, stat in enumerate(stats):
            ax = axes[i][j]
            null = (null_df.filter(pl.col("annot") == annot)[stat].to_list()
                    if stat in null_df.columns else [])
            r = real_df.filter(pl.col("annot") == annot)
            real = r[stat].item(0) if (r.height and stat in r.columns) else None
            panel(ax, real, null, stat, bins)
            if i == 0:
                ax.annotate(STAT_LABELS.get(stat, stat),
                            xy=(0.5, 1.18), xycoords="axes fraction",
                            ha="center", fontsize=12, fontweight="bold")
            if j == 0:
                ax.annotate(SET_LABELS.get(annot, annot), xy=(-0.28, 0.5), xycoords="axes fraction",
                            ha="center", va="center", rotation=90,
                            fontsize=11, fontweight="bold")
            if i == nrows - 1:
                ax.set_xlabel(STAT_LABELS.get(stat, stat), fontsize=9)

    title = f"Matched-anchor permutation null vs real — {cond}"
    if study:
        title += f"\n{study}"
    fig.suptitle(title, fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0.02, 0, 1, 0.94))
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--conditions", nargs="+", default=None,
                   help="Conditions to plot (default: all present).")
    p.add_argument("--annots", nargs="+",
                   default=["H4_genes", "regulators_H4", "targets_H4"])
    p.add_argument("--stats", nargs="+", default=["enrichment", "coefficient_z"])
    p.add_argument("--bins", type=int, default=30)
    p.add_argument("--out-dir", default=None,
                   help="Where to write figures (default: FIG_DIR/permutation).")
    args = p.parse_args()

    perperms = sorted(PERM_DIR.glob("*/permutation_ldsc_perperm.csv"))
    if not perperms:
        raise SystemExit(f"No permutation_ldsc_perperm.csv under {PERM_DIR} "
                         f"(run nf_ldsc_permutation for this design first)")
    out_root = Path(args.out_dir) if args.out_dir else FIG_DIR / "permutation"

    for perperm in perperms:
        df = pl.read_csv(perperm)
        for study in df["study"].unique(maintain_order=True).to_list():
            df_study = df.filter(pl.col("study") == study)
            conditions = args.conditions or df_study["condition"].unique().to_list()
            out_dir = out_root / study
            out_dir.mkdir(parents=True, exist_ok=True)
            for cond in conditions:
                df_cond = df_study.filter(pl.col("condition") == cond)
                if df_cond.height == 0:
                    continue
                plot_condition(df_cond, cond, args.annots, args.stats,
                               args.bins, study, out_dir / f"permutation_null_{cond}.png")


if __name__ == "__main__":
    main()
