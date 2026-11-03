#!/usr/bin/env python
"""Generate summary plots from parsed S-LDSC L2G results.

Reads the CSV produced by parse_ldsc_results.py and writes two PDFs:

  l2g_zscore_heatmap.pdf      Trait-mean coefficient z-score heatmap.
                              Rows: traits sorted by mean |z|.
                              Cols: (condition × annot) ordered by
                              Rest → Stim8hr → Stim48hr.

  l2g_enrichment_forest.pdf   Per-annot forest plot of enrichment ±
                              95% CI for each trait, with conditions
                              offset vertically.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]
COND_COLORS = {"Rest": "#1f77b4", "Stim8hr": "#ff7f0e", "Stim48hr": "#2ca02c"}


def trait_pivot(df, value_col, label_col, annot_order):
    """Mean of value_col per (label, condition, annot), wide format."""
    df_agg = df.group_by([label_col, "condition", "annot"]).agg(
        pl.col(value_col).mean()
    )
    pivot = (
        df_agg.with_columns(
            (pl.col("condition") + "_" + pl.col("annot")).alias("col_id")
        )
        .pivot(values=value_col, index=label_col, on="col_id")
        .to_pandas()
        .set_index(label_col)
    )
    desired = [f"{c}_{a}" for c in COND_ORDER for a in annot_order]
    pivot = pivot.reindex(columns=[c for c in desired if c in pivot.columns])
    return pivot


def heatmap(df, out_path, label_col, annot_order, anchor_set, vlim=4):
    pivot = trait_pivot(df, "coefficient_z", label_col, annot_order)
    # Sort rows by overall signal strength.
    pivot = pivot.reindex(
        index=pivot.abs().mean(axis=1).sort_values(ascending=False).index
    )

    # Truncate over-long trait names so y-axis stays readable.
    pivot.index = [s if len(s) <= 60 else s[:57] + "…" for s in pivot.index]

    fig, ax = plt.subplots(
        figsize=(
            max(8, 0.7 * pivot.shape[1] + 4),
            max(6, 0.35 * pivot.shape[0] + 2),
        )
    )

    im = ax.imshow(pivot.values, cmap="RdBu_r", vmin=-vlim, vmax=vlim, aspect="auto")

    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right")
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels(pivot.index, fontsize=8)

    # Vertical separators between conditions.
    n_a = len(annot_order)
    for x in range(n_a, pivot.shape[1], n_a):
        ax.axvline(x - 0.5, color="black", linewidth=1.2)

    # Annotate each cell with its value.
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.values[i, j]
            if np.isnan(v):
                continue
            ax.text(
                j,
                i,
                f"{v:.1f}",
                ha="center",
                va="center",
                fontsize=6,
                color="white" if abs(v) > 2.5 else "black",
            )

    plt.colorbar(im, ax=ax, label="Coefficient z-score (trait mean)")
    ax.set_title(f"S-LDSC partitioned heritability - {anchor_set} annotations")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close(fig)
    print(f"Wrote {out_path}")


def forest(df, out_path, label_col, annot_order, anchor_set):
    df_agg = (
        df.group_by([label_col, "condition", "annot"]).agg(
            [pl.col("enrichment").mean(), pl.col("enrichment_se").mean()]
        )
    ).to_pandas()
    df_agg["ci_low"] = df_agg["enrichment"] - 1.96 * df_agg["enrichment_se"]
    df_agg["ci_high"] = df_agg["enrichment"] + 1.96 * df_agg["enrichment_se"]

    # Sort labels by mean enrichment across all annots / conditions.
    order = (
        df_agg.groupby(label_col)["enrichment"]
        .mean()
        .sort_values(ascending=False)
        .index.tolist()
    )
    short = {t: (t if len(t) <= 50 else t[:47] + "…") for t in order}

    fig, axes = plt.subplots(
        1,
        len(annot_order),
        figsize=(15, max(6, 0.35 * len(order) + 2)),
        sharey=True,
    )

    for ax, annot in zip(axes, annot_order):
        sub = df_agg[df_agg["annot"] == annot]
        for i_cond, cond in enumerate(COND_ORDER):
            s = sub[sub["condition"] == cond].set_index(label_col).reindex(order)
            offset = (i_cond - 1) * 0.25
            y = np.arange(len(order)) + offset
            ax.errorbar(
                s["enrichment"],
                y,
                xerr=[s["enrichment"] - s["ci_low"], s["ci_high"] - s["enrichment"]],
                fmt="o",
                capsize=2,
                markersize=3.5,
                color=COND_COLORS[cond],
                ecolor=COND_COLORS[cond],
                alpha=0.85,
                label=cond,
            )

        ax.axvline(1, color="grey", linewidth=0.7, linestyle="--")
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels([short[t] for t in order], fontsize=8)
        ax.set_title(annot)
        ax.set_xlabel("Enrichment (log scale)")
        ax.set_xscale("log")
        ax.grid(True, axis="x", alpha=0.3)

    axes[-1].legend(title="Condition", loc="lower right", fontsize=8)
    plt.suptitle(
        f"S-LDSC enrichment by trait — mean across studies ({anchor_set})",
        fontsize=12,
        y=1.0,
    )
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close(fig)
    print(f"Wrote {out_path}")


def pick_label_col(df, requested):
    if requested:
        return requested
    if (
        "traitFromSource" in df.columns
        and df["traitFromSource"].null_count() < df.height
    ):
        return "traitFromSource"
    if "trait_short" in df.columns:
        return "trait_short"
    return "study"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--summary", required=True, help="CSV from parse_ldsc_results.py"
    )
    parser.add_argument(
        "--out-dir", default=None, help="Output dir (defaults to alongside --summary)"
    )
    parser.add_argument(
        "--annots",
        nargs="+",
        default=None,
        help="Annotation names + display order. "
        "Default: alphabetical from the data.",
    )
    parser.add_argument(
        "--label-col",
        default=None,
        help="Column to use as row labels "
        "(default: traitFromSource → trait_short → study).",
    )
    parser.add_argument(
        "--prefix",
        default="ldsc",
        help="Output filename prefix (default 'ldsc'; e.g. 'l2g' or 'coloc').",
    )
    args = parser.parse_args()

    df = pl.read_csv(args.summary)

    label_col = pick_label_col(df, args.label_col)
    annot_order = args.annots or sorted(df["annot"].unique().to_list())

    print(
        f"Loaded {df.height} rows from "
        f"{df[label_col].n_unique()} {label_col}s across "
        f"{df['study'].n_unique()} studies. "
        f"Annotations: {annot_order}."
    )

    out_dir = Path(args.out_dir) if args.out_dir else Path(args.summary).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    heatmap(
        df,
        out_dir / f"{args.prefix}_zscore_heatmap.pdf",
        label_col=label_col,
        annot_order=annot_order,
        anchor_set=args.prefix,
    )
    forest(
        df,
        out_dir / f"{args.prefix}_enrichment_forest.pdf",
        label_col=label_col,
        annot_order=annot_order,
        anchor_set=args.prefix,
    )


if __name__ == "__main__":
    main()
