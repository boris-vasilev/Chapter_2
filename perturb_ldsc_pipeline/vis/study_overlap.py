#!/usr/bin/env python
"""Pairwise overlap of pipeline gene sets across studies / traits.

For each (condition, annot) faceted panel, compute the pairwise Jaccard
index and overlap coefficient between every pair of studies and plot as
heatmaps.

Works with both gene-set modules. Filenames are expected as
``{study}_{condition}_{annot}.txt`` where ``annot`` is two
underscore-joined tokens (e.g. ``H4_genes``, ``regulators_l2g``).

Examples
--------
Coloc (18 studies):

    python vis/study_overlap.py \
        --gene-sets-dir data/gene_sets/coloc \
        --study-info    data/GWAS_study_info.csv \
        --annots H4_genes regulators_H4 targets_H4 \
        --out-prefix    plots/gene_set_analysis/coloc

L2G (collapse to ~17 traits via traitFromSource):

    python vis/study_overlap.py \
        --gene-sets-dir data/gene_sets/l2g \
        --study-info    data/opentargets/GWAS_study_info.csv \
        --annots l2g_genes regulators_l2g targets_l2g \
        --label-col     traitFromSource \
        --out-prefix    plots/gene_set_analysis/l2g
"""

import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

CONDITIONS = ["Rest", "Stim8hr", "Stim48hr"]


def parse_filename(fname):
    """{study}_{condition}_{annot1}_{annot2}.txt → (study, condition, annot)."""
    stem = fname.replace(".txt", "")
    parts = stem.split("_")
    if len(parts) < 4:
        return None
    annot = "_".join(parts[-2:])
    condition = parts[-3]
    study = "_".join(parts[:-3])
    return study, condition, annot


def load_gene_sets(gene_sets_dir, annots_keep):
    """Walk the directory and return a dict (study, condition, annot) → set[str]."""
    out = {}
    for fname in sorted(os.listdir(gene_sets_dir)):
        if not fname.endswith(".txt"):
            continue
        parsed = parse_filename(fname)
        if parsed is None:
            continue
        study, condition, annot = parsed
        if annot not in annots_keep:
            continue
        with open(Path(gene_sets_dir) / fname) as f:
            genes = {line.strip() for line in f if line.strip()}
        out[(study, condition, annot)] = genes
    return out


def label_studies(gene_sets, study_info_path, label_col):
    """Return a dict (label, condition, annot) → set[str], merging studies with
    the same label by taking the gene-set union.
    """
    if study_info_path is None or label_col == "study":
        return {k: v for k, v in gene_sets.items()}

    info = pl.read_csv(study_info_path)
    if label_col not in info.columns:
        raise SystemExit(
            f"{label_col} not in {study_info_path}. " f"Available: {info.columns}"
        )
    mapping = dict(zip(info["study"].to_list(), info[label_col].to_list()))

    out = {}
    for (study, cond, annot), genes in gene_sets.items():
        label = mapping.get(study, study)
        if label is None:
            label = study  # fall back if traitFromSource is null
        key = (label, cond, annot)
        out[key] = out.get(key, set()) | genes
    return out


def overlap_matrices(records, labels, conditions, annots):
    """Return {(cond, annot): (jaccard_df, overlap_df)}."""
    n = len(labels)
    pos = {s: i for i, s in enumerate(labels)}
    result = {}
    for cond in conditions:
        for annot in annots:
            jac = np.full((n, n), np.nan)
            ovl = np.full((n, n), np.nan)
            for i, si in enumerate(labels):
                a = records.get((si, cond, annot))
                if a is None:
                    continue
                for j, sj in enumerate(labels):
                    b = records.get((sj, cond, annot))
                    if b is None:
                        continue
                    inter = len(a & b)
                    union = len(a | b)
                    jac[i, j] = inter / union if union else np.nan
                    smaller = min(len(a), len(b))
                    ovl[i, j] = inter / smaller if smaller else np.nan
            result[(cond, annot)] = (
                pd.DataFrame(jac, index=labels, columns=labels),
                pd.DataFrame(ovl, index=labels, columns=labels),
            )
    return result


def truncate(s, n=40):
    return s if len(s) <= n else s[: n - 1] + "…"


def plot_heatmaps(
    matrices, conditions, annots, labels, value_idx, title, cbar_label, out_path
):
    nrows = len(conditions)
    ncols = len(annots)
    side = max(1.0, 0.45 * len(labels))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(ncols * (side + 1.5), nrows * (side + 1.5)),
        squeeze=False,
    )

    short = [truncate(l) for l in labels]

    for r, cond in enumerate(conditions):
        for c, annot in enumerate(annots):
            ax = axes[r, c]
            mat = matrices[(cond, annot)][value_idx].copy()
            # cluster using 1 - Jaccard as distance
            dist = 1 - mat.fillna(0).values
            np.fill_diagonal(dist, 0)
            dist = np.clip(dist, 0, None)  # guard against float noise
            condensed = squareform(dist)
            Z = linkage(condensed, method="average")
            order = leaves_list(Z)

            mat = mat.iloc[order, :].iloc[:, order]
            mat.index = [truncate(labels[i]) for i in order]
            mat.columns = [truncate(labels[i]) for i in order]

            mask = np.eye(len(order), dtype=bool)
            sns.heatmap(
                mat,
                ax=ax,
                mask=mask,
                vmin=0,
                vmax=1,
                cmap="YlOrRd",
                linewidths=0.4,
                linecolor="white",
                annot=(len(labels) <= 25),
                fmt=".2f",
                annot_kws={"size": 7},
                cbar_kws={"label": cbar_label, "shrink": 0.7},
            )

    fig.suptitle(title, fontsize=14, y=1.01)
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gene-sets-dir",
        required=True,
        help="e.g. data/gene_sets/coloc or data/gene_sets/l2g",
    )
    parser.add_argument(
        "--study-info",
        default=None,
        help="CSV with at least a 'study' column and the column "
        "named in --label-col. Required unless --label-col=study.",
    )
    parser.add_argument(
        "--annots",
        nargs="+",
        required=True,
        help="Annotation names (e.g. H4_genes regulators_H4 targets_H4)",
    )
    parser.add_argument("--conditions", nargs="+", default=CONDITIONS)
    parser.add_argument(
        "--label-col",
        default="study",
        help="Column in --study-info to use as row/col labels. "
        "Studies sharing the same label are merged "
        "(gene-set union). Default: study.",
    )
    parser.add_argument(
        "--out-prefix",
        required=True,
        help="Output prefix; writes {prefix}_jaccard.png and " "{prefix}_overlap.png.",
    )
    args = parser.parse_args()

    raw = load_gene_sets(args.gene_sets_dir, set(args.annots))
    if not raw:
        raise SystemExit(f"No matching gene-set files in {args.gene_sets_dir}")

    records = label_studies(raw, args.study_info, args.label_col)
    labels = sorted({k[0] for k in records.keys()})
    print(f"Loaded {len(raw)} gene sets → {len(labels)} unique '{args.label_col}'")

    matrices = overlap_matrices(records, labels, args.conditions, args.annots)

    out_prefix = Path(args.out_prefix)
    plot_heatmaps(
        matrices,
        args.conditions,
        args.annots,
        labels,
        value_idx=0,
        title="Pairwise Jaccard index between studies",
        cbar_label="Jaccard index",
        out_path=out_prefix.with_name(out_prefix.name + "_jaccard.png"),
    )
    plot_heatmaps(
        matrices,
        args.conditions,
        args.annots,
        labels,
        value_idx=1,
        title="Pairwise overlap coefficient between studies",
        cbar_label="Overlap coefficient",
        out_path=out_prefix.with_name(out_prefix.name + "_overlap.png"),
    )


if __name__ == "__main__":
    main()
