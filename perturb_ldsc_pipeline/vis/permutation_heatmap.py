#!/usr/bin/env python
"""Plot a permutation-overlay version of the drug-target enrichment heatmap.

Same colour scale as the Fisher heatmap (`approved_drug_enrichment_all.py`):
log2(real OR) on RdBu_r, capped at ±4. The only thing that changes is the
*significance overlay*:

  - asterisks: empirical p thresholds (* < 0.05, ** < 0.01, *** < 0.001)
  - gold border: empirical p < 0.05

Sat next to the Fisher heatmap, the two figures show the same numerical
result; the difference is which cells survive each kind of test.

Reads `permutation_results.csv` (produced by
`approved_drug_enrichment_permutation.py`).
"""

import argparse
from pathlib import Path

import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

IMD_DISEASES = [
    "AST",
    "ASTAO",
    "ASTCO",
    "ASPH",
    "ATD",
    "AITD",
    "CD",
    "IBD",
    "UC",
    "RA",
    "SLE",
    "T1D",
]
NON_IMD = ["IS", "MYOI", "OA", "T2D", "ADFH", "SCZ"]

SHORT_NAMES: dict[str, str] = {
    "AST": "Asthma",
    "ASTAO": "Asthma (AO)",
    "ASTCO": "Asthma (CO)",
    "ATD": "Atopic derm.",
    "AITD": "Autoimm. thyroid",
    "CD": "Crohn's",
    "IBD": "IBD",
    "UC": "Ulcer. colitis",
    "RA": "Rheum. arthritis",
    "SLE": "SLE",
    "T1D": "T1D",
    "IS": "Ischemic stroke",
    "MYOI": "Myocard. infarct.",
    "OA": "Osteoarthritis",
    "T2D": "T2D",
    "ADFH": "Alzheimer's",
    "SCZ": "Schizophrenia",
}


def plot_heatmap(df: pl.DataFrame, out_path: Path, lfc_threshold=None) -> None:
    condition_order = ["Rest", "Stim8hr", "Stim48hr"]
    gene_sets_present = set(df.get_column("gene_set").unique().to_list())
    # Coloc (anchor) is condition-independent -> one collapsed row on top, if
    # present. Regulators / Targets get one row per condition.
    row_keys: list[tuple[str, object]] = []
    row_labels: list[str] = []
    if "Coloc" in gene_sets_present:
        row_keys.append(("Coloc", None))
        row_labels.append("Coloc genes")
    for gs in ["Regulators", "Targets"]:
        for cond in condition_order:
            row_keys.append((gs, cond))
            row_labels.append(f"{gs} | {cond}")

    present = df.get_column("disease").unique().to_list()
    disease_order = [d for d in IMD_DISEASES + NON_IMD if d in present]

    emp_p = np.full((len(row_keys), len(disease_order)), np.nan)
    or_matrix = np.full((len(row_keys), len(disease_order)), np.nan)
    for i, (gs, cond) in enumerate(row_keys):
        for j, disease in enumerate(disease_order):
            flt = (pl.col("gene_set") == gs) & (pl.col("disease") == disease)
            if cond is not None:  # Coloc row collapses across conditions
                flt = flt & (pl.col("condition") == cond)
            sub = df.filter(flt)
            if sub.height > 0:
                v = sub.get_column("empirical_p").item(0)
                if v is not None and not np.isnan(v):
                    emp_p[i, j] = v
                o = sub.get_column("real_or").item(0)
                if o is not None and not np.isnan(o):
                    or_matrix[i, j] = o

    # Same colour scale as the Fisher heatmap: log2(OR), capped at ±4
    # (inf → 16 before log2). Lets the two figures sit side-by-side without
    # the viewer having to recalibrate.
    capped = np.where(np.isinf(or_matrix), 16.0, or_matrix)
    log2_or = np.log2(np.clip(capped, 0.25, 16))

    fig_w = max(10, 0.85 * len(disease_order) + 4)
    fig, ax = plt.subplots(figsize=(fig_w, 7))
    im = ax.imshow(log2_or, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")

    # Mark cells with no defined OR as no-data dots (matches Fisher heatmap).
    for i in range(len(row_keys)):
        for j in range(len(disease_order)):
            if np.isnan(or_matrix[i, j]):
                ax.text(j, i, "·", ha="center", va="center", color="grey", fontsize=18)

    # Add asterisks for emp_p thresholds (white text with black stroke,
    # mirroring the Fisher heatmap's convention).
    stroke = [
        path_effects.Stroke(linewidth=2.5, foreground="black"),
        path_effects.Normal(),
    ]
    for i in range(len(row_keys)):
        for j in range(len(disease_order)):
            v = emp_p[i, j]
            if np.isnan(v):
                continue
            stars = "***" if v < 0.001 else "**" if v < 0.01 else "*" if v < 0.1 else ""
            if stars:
                txt = ax.text(
                    j,
                    i,
                    stars,
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=20,
                    fontweight="bold",
                )
                txt.set_path_effects(stroke)

    # IMD / non-IMD divider.
    n_imd = sum(1 for d in disease_order if d in IMD_DISEASES)
    if 0 < n_imd < len(disease_order):
        ax.axvline(n_imd - 0.5, color="black", linewidth=2.0)
        ax.text(
            n_imd / 2 - 0.5,
            -0.85,
            "Immune-mediated diseases",
            ha="center",
            fontsize=18,
            fontweight="bold",
        )
        ax.text(
            (n_imd + len(disease_order)) / 2 - 0.5,
            -0.85,
            "Other diseases",
            ha="center",
            fontsize=18,
            fontweight="bold",
        )

    ax.set_xticks(range(len(disease_order)))
    ax.set_xticklabels(
        [SHORT_NAMES.get(d, d) for d in disease_order],
        rotation=45,
        ha="right",
        fontsize=18,
    )
    ax.set_yticks(range(len(row_keys)))
    ax.set_yticklabels(row_labels, fontsize=18)

    cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("log2(OR), capped at ±4", fontsize=16)
    cbar.ax.tick_params(labelsize=12)

    if lfc_threshold is None:
        thr_note = ""
    elif lfc_threshold > 0:
        thr_note = f"  |logFC| ≥ {lfc_threshold:g}"
    else:
        thr_note = "  (no edge filter)"
    ax.set_title(
        "Drug-target enrichment with permutation-based significance "
        f"(size + property-matched){thr_note}\n"
        "* emp_p<0.05, ** emp_p<0.01, *** emp_p<0.001; "
        "; · = no data",
        pad=35,
        fontsize=16,
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=400, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-csv", required=True)
    parser.add_argument("--output", required=True, help="PNG output path")
    parser.add_argument(
        "--min-abs-logfc",
        type=float,
        default=None,
        help="|logFC| edge threshold for the title. Default: read from the "
        "CSV's min_abs_logfc column if present (written by "
        "approved_drug_enrichment_permutation.py).",
    )
    args = parser.parse_args()

    df = pl.read_csv(args.results_csv)
    lfc = args.min_abs_logfc
    if lfc is None and "min_abs_logfc" in df.columns:
        vals = [
            v
            for v in df.get_column("min_abs_logfc").unique().to_list()
            if v is not None
        ]
        lfc = vals[0] if len(vals) == 1 else None
    plot_heatmap(df, Path(args.output), lfc_threshold=lfc)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
