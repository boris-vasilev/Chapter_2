#!/usr/bin/env python
"""Plot the drug-target enrichment permutation null (Fisher ORs) with the real OR.

Analogue of perturb_ldsc_permutation/scripts/plot_permutation_null.py, but for the
matched-anchor permutation in approved_drug_enrichment_permutation.py. Instead of
the S-LDSC per-perm CSV, it reads the two artefacts that script emits:

  * permutation_nulls.parquet  (requires --save-null-draws): one null OR per
    (disease, gene_set, condition, draw).
  * permutation_results.csv: the real_or per (disease, gene_set, condition).

For each disease it draws a grid of histograms: one row per gene set
(Coloc / Regulators / Targets, whichever are present) and one column per
condition (Rest / Stim8hr / Stim48hr). Each panel shows the null OR distribution
with the real OR as a red vertical line, annotated with the one-sided empirical p
and n (after dropping NaN null draws, where Fisher's OR is undefined).

Empirical p = (#null >= real + 1) / (n_null + 1), matching the permutation
script (real should sit in the upper tail of the matched-random null).

Usage:
  python perturb_ldsc_pipeline/vis/plot_drug_enrichment_null.py \
      --perm-dir plots/approved_drug_enrichment/permutation \
      --results permutation_results0.5.csv \
      --diseases ASTCO --conditions Stim48hr
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless / cluster
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import polars as pl  # noqa: E402

GENE_SET_ORDER = ["Coloc", "Regulators", "Targets"]
COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]


def valid_null(null) -> np.ndarray:
    """Drop only NaN (undefined Fisher OR, e.g. empty pseudo-set). Crucially
    KEEP +inf: it is a genuine OR (set ⊆ drug targets) and counts as >= real
    in the one-sided empirical p, exactly as in the permutation script."""
    arr = np.asarray(null, dtype=float)
    return arr[~np.isnan(arr)]


def emp_p(real, null: np.ndarray):
    """One-sided empirical p over the NaN-dropped null (inf retained)."""
    if real is None or np.isnan(real) or null.size == 0:
        return float("nan"), int(null.size)
    return float((np.sum(null >= real) + 1) / (null.size + 1)), int(null.size)


def panel(ax, real, null: np.ndarray, bins: int):
    # `null` is NaN-dropped but may contain +inf draws.
    if null.size == 0:
        ax.text(0.5, 0.5, "no null data", ha="center", va="center",
                transform=ax.transAxes, color="grey")
        ax.set_yticks([])
        return
    p, n = emp_p(real, null)
    finite = null[np.isfinite(null)]
    n_inf = int(np.isinf(null).sum())
    if finite.size:
        ax.hist(finite, bins=bins, color="#4C72B0", alpha=0.75,
                edgecolor="white", linewidth=0.3)
    if real is not None and np.isfinite(real):
        ax.axvline(real, color="#C44E52", lw=2.2, label=f"real OR = {real:.2f}")
    # Make the inf draws visible: a hatched bar pinned to the right edge.
    if n_inf:
        x0, x1 = ax.get_xlim()
        w = (x1 - x0) * 0.03
        ax.bar(x1, n_inf, width=w, align="edge", color="#C44E52", alpha=0.35,
               hatch="//", edgecolor="#C44E52", label=f"OR=∞ ×{n_inf}")
        ax.set_xlim(x0, x1 + w * 1.5)
    if (real is not None and np.isfinite(real)) or n_inf:
        ax.legend(fontsize=8, loc="upper right", frameon=False)
    if real is not None and np.isfinite(real):
        ptxt = "p < 0.001" if p < 1e-3 else f"p = {p:.3f}"
        ax.set_title(f"{ptxt}   (n={n})", fontsize=9)
    else:
        ax.set_title(f"no/undefined real OR (n={n})", fontsize=9, color="grey")
    ax.set_ylabel("permutations", fontsize=8)


def null_draws(nulls: pl.DataFrame, disease, gene_set, cond) -> np.ndarray:
    sub = nulls.filter(
        (pl.col("disease") == disease) & (pl.col("gene_set") == gene_set)
    )
    # Coloc draws are condition-agnostic (stored as condition == "all").
    if sub.filter(pl.col("condition") == cond).height:
        sub = sub.filter(pl.col("condition") == cond)
    elif sub.filter(pl.col("condition") == "all").height:
        sub = sub.filter(pl.col("condition") == "all")
    else:
        return np.array([])
    return valid_null(sub.get_column("null_or").to_list())


def real_or(results: pl.DataFrame, disease, gene_set, cond):
    r = results.filter(
        (pl.col("disease") == disease)
        & (pl.col("gene_set") == gene_set)
        & (pl.col("condition") == cond)
    )
    if r.height == 0:
        r = results.filter(
            (pl.col("disease") == disease) & (pl.col("gene_set") == gene_set)
        )
    return r.get_column("real_or").item(0) if r.height else None


def plot_disease(nulls, results, disease, display, gene_sets, conds, bins, out_path):
    nrows, ncols = len(gene_sets), len(conds)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 2.8 * nrows),
                             squeeze=False)
    for i, gs in enumerate(gene_sets):
        for j, cond in enumerate(conds):
            ax = axes[i][j]
            null = null_draws(nulls, disease, gs, cond)
            real = real_or(results, disease, gs, cond)
            panel(ax, real, null, bins)
            if i == 0:
                ax.annotate(cond, xy=(0.5, 1.18), xycoords="axes fraction",
                            ha="center", fontsize=12, fontweight="bold")
            if j == 0:
                ax.annotate(gs, xy=(-0.28, 0.5), xycoords="axes fraction",
                            ha="center", va="center", rotation=90,
                            fontsize=11, fontweight="bold")
            if i == nrows - 1:
                ax.set_xlabel("Odds ratio", fontsize=9)

    title = "Matched-anchor permutation null vs real — drug-target enrichment"
    title += f"\n{display or disease}"
    fig.suptitle(title, fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0.02, 0, 1, 0.93))
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--perm-dir", required=True,
                   help="Dir holding permutation_nulls.parquet and the results CSV.")
    p.add_argument("--nulls", default="permutation_nulls.parquet",
                   help="Null-draws parquet name (relative to --perm-dir or a path).")
    p.add_argument("--results", default="permutation_results.csv",
                   help="Results CSV name (relative to --perm-dir or a path).")
    p.add_argument("--diseases", nargs="+", default=None,
                   help="Disease codes to plot (default: all present).")
    p.add_argument("--gene-sets", nargs="+", default=None,
                   help="Gene sets to show as rows (default: present, in "
                        "Coloc/Regulators/Targets order).")
    p.add_argument("--conditions", nargs="+", default=None,
                   help="Conditions as columns (default: present, Rest/Stim8hr/Stim48hr).")
    p.add_argument("--bins", type=int, default=40)
    p.add_argument("--out-dir", default=None,
                   help="Where to write figures (default: --perm-dir).")
    args = p.parse_args()

    perm_dir = Path(args.perm_dir)
    nulls_path = Path(args.nulls) if Path(args.nulls).is_absolute() else perm_dir / args.nulls
    results_path = Path(args.results) if Path(args.results).is_absolute() else perm_dir / args.results
    if not nulls_path.exists():
        raise SystemExit(f"Null draws not found: {nulls_path} "
                         f"(re-run the permutation with --save-null-draws)")
    if not results_path.exists():
        raise SystemExit(f"Results CSV not found: {results_path}")

    nulls = pl.read_parquet(nulls_path)
    results = pl.read_csv(results_path)

    diseases = args.diseases or results["disease"].unique().to_list()
    present_gs = set(nulls["gene_set"].unique().to_list())
    gene_sets = args.gene_sets or [g for g in GENE_SET_ORDER if g in present_gs]
    present_conds = set(
        c for c in nulls["condition"].unique().to_list() if c != "all"
    )
    conds = args.conditions or [c for c in COND_ORDER if c in present_conds]

    displays = dict(
        zip(results["disease"].to_list(), results["display"].to_list())
    ) if "display" in results.columns else {}

    out_dir = Path(args.out_dir) if args.out_dir else perm_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    for disease in diseases:
        plot_disease(nulls, results, disease, displays.get(disease),
                     gene_sets, conds, args.bins,
                     out_dir / f"drug_enrichment_null_{disease}.png")


if __name__ == "__main__":
    main()
