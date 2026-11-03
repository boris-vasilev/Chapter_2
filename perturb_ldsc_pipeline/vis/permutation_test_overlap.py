#!/usr/bin/env python
"""Permutation test: do IMD anchor-expanded gene sets overlap more than
size-matched random anchors from the Perturb-seq universe?

For each (condition, expansion type ∈ {regulators, targets}) we:

  1. Compute the mean pairwise Jaccard between observed expansions
     across IMDs.
  2. Generate N permutations where each disease's anchors are replaced
     with a same-size random draw from the Perturb-seq gene pool
     (union of cis-perturbed and trans-measured genes in
     DE_results_sig.csv). Re-expand through the actual GRN.
  3. Empirical p-value = fraction of permutations where the null mean
     Jaccard ≥ observed mean Jaccard.

Outputs three figures + a summary CSV:

  <prefix>_null_distribution.png   Histograms of null mean Jaccard with
                                   observed value marked.
  <prefix>_excess_heatmap.png      Per-pair observed - median(null)
                                   Jaccard. Where the signal lives.
  <prefix>_hub_diagnostic.png      Cis out-degree + trans in-degree
                                   distributions per condition.
  <prefix>_summary.csv             obs/null/p_value per (cond × annot).

Examples
--------
Coloc, H4-anchored:

    python vis/permutation_test_overlap.py \
        --gene-sets-dir data/gene_sets/coloc \
        --de-results    data/DE_results_sig.csv \
        --anchor-annot  H4_genes \
        --out-prefix    plots/gene_set_analysis/coloc_h4

L2G, trait-collapsed via OpenTargets traitFromSource:

    python vis/permutation_test_overlap.py \
        --gene-sets-dir data/gene_sets/l2g \
        --de-results    data/DE_results_sig.csv \
        --anchor-annot  l2g_genes \
        --study-info    data/opentargets/GWAS_study_info.csv \
        --label-col     traitFromSource \
        --out-prefix    plots/gene_set_analysis/l2g
"""

import argparse
import os
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
import seaborn as sns

from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

CONDITIONS = ["Rest", "Stim8hr", "Stim48hr"]


# ---------- I/O ----------


def parse_filename(fname):
    """`{study}_{condition}_{annot1}_{annot2}.txt` → (study, cond, annot)."""
    stem = fname.replace(".txt", "")
    parts = stem.split("_")
    if len(parts) < 4:
        return None
    return "_".join(parts[:-3]), parts[-3], "_".join(parts[-2:])


def load_anchor_sets(gene_sets_dir, anchor_annot, conditions, label_map=None):
    """Return {(label, condition) → set[gene_id]} unioned across same-label studies."""
    out = defaultdict(set)
    for fname in sorted(os.listdir(gene_sets_dir)):
        if not fname.endswith(".txt"):
            continue
        parsed = parse_filename(fname)
        if parsed is None:
            continue
        study, condition, annot = parsed
        if annot != anchor_annot or condition not in conditions:
            continue
        label = label_map.get(study, study) if label_map else study
        if label is None:
            label = study
        with open(Path(gene_sets_dir) / fname) as f:
            out[(label, condition)] |= {line.strip() for line in f if line.strip()}
    return dict(out)


def build_grn(de_csv):
    """Per-condition GRN: returns (trans_to_cis, cis_to_trans, pool).

    trans_to_cis[cond][trans] = set of cis genes whose perturbation has trans as a DEG
    cis_to_trans[cond][cis]   = set of trans DEGs when cis is perturbed
    pool                      = union of all cis and trans genes across all conditions
    """
    de = pl.read_csv(de_csv).filter(pl.col("cis_ensg") != pl.col("trans_ensg"))
    trans_to_cis = defaultdict(lambda: defaultdict(set))
    cis_to_trans = defaultdict(lambda: defaultdict(set))
    pool = set()
    for row in de.iter_rows(named=True):
        cond = row["culture_condition"]
        cis = row["cis_ensg"]
        trans = row["trans_ensg"]
        trans_to_cis[cond][trans].add(cis)
        cis_to_trans[cond][cis].add(trans)
        pool.add(cis)
        pool.add(trans)
    return (
        {c: dict(d) for c, d in trans_to_cis.items()},
        {c: dict(d) for c, d in cis_to_trans.items()},
        pool,
    )


# ---------- Expansion ----------


def regulators_of(anchor_set, trans_to_cis_cond):
    out = set()
    for g in anchor_set:
        out |= trans_to_cis_cond.get(g, set())
    return out


def targets_of(anchor_set, cis_to_trans_cond):
    out = set()
    for g in anchor_set:
        out |= cis_to_trans_cond.get(g, set())
    return out


# ---------- Test statistics ----------


def pairwise_jaccard_matrix(sets_dict, labels):
    n = len(labels)
    out = np.full((n, n), np.nan)
    for i, li in enumerate(labels):
        a = sets_dict.get(li, set())
        for j, lj in enumerate(labels):
            b = sets_dict.get(lj, set())
            union = a | b
            if union:
                out[i, j] = len(a & b) / len(union)
    return out


def mean_upper_triangle(mat):
    """Mean of upper triangle excluding diagonal, NaN-safe."""
    n = mat.shape[0]
    iu = np.triu_indices(n, k=1)
    vals = mat[iu]
    return np.nanmean(vals) if np.any(~np.isnan(vals)) else 0.0


def run_test(anchors_cond, expand_fn, network_cond, pool_arr, n_perm, rng):
    """For one (condition, expansion-type): observed + null distribution."""
    labels = sorted(anchors_cond.keys())
    real_expanded = {l: expand_fn(anchors_cond[l], network_cond) for l in labels}
    obs_mat = pairwise_jaccard_matrix(real_expanded, labels)
    obs_mean = mean_upper_triangle(obs_mat)

    null_means = np.empty(n_perm)
    null_mats = np.empty((n_perm, len(labels), len(labels)))
    for p in range(n_perm):
        perm_expanded = {}
        for l in labels:
            k = len(anchors_cond[l])
            if k == 0:
                perm_expanded[l] = set()
                continue
            sampled = set(rng.choice(pool_arr, size=k, replace=False))
            perm_expanded[l] = expand_fn(sampled, network_cond)
        m = pairwise_jaccard_matrix(perm_expanded, labels)
        null_mats[p] = m
        null_means[p] = mean_upper_triangle(m)

    null_median_mat = np.nanmedian(null_mats, axis=0)
    return labels, obs_mean, obs_mat, null_means, null_median_mat


# ---------- Plotting ----------


def plot_null_hist(ax, obs, null, p_val, title):
    ax.hist(null, bins=40, color="lightgray", edgecolor="white")
    ax.axvline(
        obs, color="red", linewidth=2, label=f"observed = {obs:.3f}\np = {p_val:.3g}"
    )
    ax.set_xlabel("Mean pairwise Jaccard")
    ax.set_ylabel("Permutations")
    ax.set_title(title, fontsize=10)
    ax.legend(loc="upper right", fontsize=8)


def plot_excess(ax, obs_mat, null_med_mat, labels, title, vlim=0.3):
    excess = obs_mat - null_med_mat
    short = [l if len(l) <= 30 else l[:29] + "…" for l in labels]
    df = pd.DataFrame(excess, index=short, columns=short)
    mask = np.eye(len(labels), dtype=bool)
    sns.heatmap(
        df,
        mask=mask,
        ax=ax,
        cmap="RdBu_r",
        center=0,
        vmin=-vlim,
        vmax=vlim,
        cbar_kws={"label": "J_obs - median J_null", "shrink": 0.6},
    )
    ax.set_title(title, fontsize=10, pad=4)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=6)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=6)
    ax.set_xlabel("")
    ax.set_ylabel("")


def plot_excess_clustered(ax, obs_mat, null_med_mat, labels, title, vlim=0.3):
    excess = obs_mat - null_med_mat

    # Use observed Jaccard for clustering
    sim = obs_mat.copy()

    # Replace NaNs with 0 similarity
    sim = np.nan_to_num(sim, nan=0.0)

    # Force diagonal = perfect self-similarity
    np.fill_diagonal(sim, 1.0)

    # Convert to distance
    dist = 1.0 - sim

    # Numerical safety
    dist = np.clip(dist, 0, 1)

    # Condensed distance vector
    condensed = squareform(dist, checks=False)

    # Cluster
    Z = linkage(condensed, method="average")
    order = leaves_list(Z)

    # Reorder everything
    excess = excess[np.ix_(order, order)]
    ordered_labels = [labels[i] for i in order]

    short = [l if len(l) <= 30 else l[:29] + "…" for l in ordered_labels]

    df = pd.DataFrame(excess, index=short, columns=short)

    mask = np.eye(len(labels), dtype=bool)

    sns.heatmap(
        df,
        mask=mask,
        ax=ax,
        cmap="RdBu_r",
        center=0,
        vmin=-vlim,
        vmax=vlim,
        cbar_kws={"label": "J_obs - median J_null", "shrink": 0.6},
    )

    ax.set_title(title, fontsize=10, pad=4)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=6)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=6)
    ax.set_xlabel("")
    ax.set_ylabel("")


def hub_diagnostic(cis_to_trans, trans_to_cis, conditions, out_path):
    fig, axes = plt.subplots(
        2, len(conditions), figsize=(5 * len(conditions), 7), squeeze=False
    )
    for c, cond in enumerate(conditions):
        out_deg = [len(v) for v in cis_to_trans.get(cond, {}).values()]
        in_deg = [len(v) for v in trans_to_cis.get(cond, {}).values()]

        axes[0, c].hist(out_deg, bins=50, color="steelblue", edgecolor="white")
        axes[0, c].set_xscale("symlog")
        axes[0, c].set_xlabel("Out-degree (# trans DEGs)")
        axes[0, c].set_ylabel("# cis perturbations")
        axes[0, c].set_title(
            f"{cond} — cis out-degree "
            f"(median={int(np.median(out_deg))}, "
            f"max={max(out_deg)})",
            fontsize=9,
        )

        axes[1, c].hist(in_deg, bins=50, color="darkorange", edgecolor="white")
        axes[1, c].set_xscale("symlog")
        axes[1, c].set_xlabel("In-degree (# perturbations affecting this gene)")
        axes[1, c].set_ylabel("# trans genes")
        axes[1, c].set_title(
            f"{cond} — trans in-degree "
            f"(median={int(np.median(in_deg))}, "
            f"max={max(in_deg)})",
            fontsize=9,
        )

    fig.suptitle("Perturb-seq GRN hub diagnostic", y=1.01, fontsize=12)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


# ---------- Main ----------


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-sets-dir", required=True)
    parser.add_argument(
        "--de-results", required=True, help="DE_results_sig.csv from the pipeline"
    )
    parser.add_argument(
        "--anchor-annot",
        required=True,
        help="Anchor annot name (e.g. H4_genes, l2g_genes)",
    )
    parser.add_argument("--study-info", default=None)
    parser.add_argument("--label-col", default="study")
    parser.add_argument("--conditions", nargs="+", default=CONDITIONS)
    parser.add_argument("--n-perm", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    # --- Build GRN ---
    print(f"Building GRN from {args.de_results}")
    trans_to_cis, cis_to_trans, pool = build_grn(args.de_results)
    pool_arr = np.array(sorted(pool))
    print(f"  pool size: {len(pool)} genes")
    for cond in args.conditions:
        n_edges = sum(len(v) for v in cis_to_trans.get(cond, {}).values())
        n_cis = len(cis_to_trans.get(cond, {}))
        n_trans = len(trans_to_cis.get(cond, {}))
        print(f"  {cond:>8s}: {n_edges:>7d} edges, {n_cis} cis, {n_trans} trans")

    # --- Anchor sets ---
    print(f"\nLoading anchor sets ({args.anchor_annot}) from {args.gene_sets_dir}")
    label_map = None
    if args.study_info and args.label_col != "study":
        info = pl.read_csv(args.study_info)
        if args.label_col in info.columns:
            label_map = dict(
                zip(info["study"].to_list(), info[args.label_col].to_list())
            )
        else:
            raise SystemExit(f"--label-col {args.label_col!r} not in study_info")

    raw = load_anchor_sets(
        args.gene_sets_dir, args.anchor_annot, args.conditions, label_map
    )

    # Filter anchors to pool so real and permuted have matched effective sizes.
    anchors = {}
    dropped_total = 0
    kept_total = 0
    for key, genes in raw.items():
        in_pool = genes & pool
        anchors[key] = in_pool
        dropped_total += len(genes) - len(in_pool)
        kept_total += len(in_pool)
    print(f"  Loaded {len(raw)} (label × condition) anchor sets")
    print(
        f"  In-pool / total anchor genes: {kept_total} / "
        f"{kept_total + dropped_total} "
        f"({100*kept_total/(kept_total+dropped_total+1e-9):.1f}%)"
    )

    # Per-label anchor sizes
    by_label = defaultdict(dict)
    for (l, c), g in anchors.items():
        by_label[l][c] = len(g)
    print("\nAnchor sizes after pool intersection:")
    print(f"  {'label':<35s} " + "  ".join(f"{c:>8s}" for c in args.conditions))
    for label in sorted(by_label):
        sizes = "  ".join(f"{by_label[label].get(c, 0):>8d}" for c in args.conditions)
        short = label if len(label) <= 33 else label[:32] + "…"
        print(f"  {short:<35s} {sizes}")

    # --- Run tests ---
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    EXPAND = [
        ("regulators", regulators_of, trans_to_cis),
        ("targets", targets_of, cis_to_trans),
    ]

    n_cond = len(args.conditions)
    fig_h, axes_h = plt.subplots(2, n_cond, figsize=(5 * n_cond, 8), squeeze=False)
    # Excess heatmap figure: panel size scales with n labels
    n_lab = len({l for (l, _) in anchors.keys()})
    panel = max(4, 0.35 * n_lab + 1.5)
    fig_e, axes_e = plt.subplots(
        2, n_cond, figsize=(panel * n_cond + 2, panel * 2 + 1), squeeze=False
    )

    rows = []
    for c, cond in enumerate(args.conditions):
        anchors_cond = {
            l: anchors[(l, cond)]
            for (l, c2) in anchors.keys()
            if c2 == cond and anchors[(l, c2)]
        }
        if not anchors_cond:
            print(f"\n[skip] {cond}: no non-empty anchor sets")
            continue

        for r, (name, expand_fn, network) in enumerate(EXPAND):
            print(f"\n=== {cond} | {name} ({len(anchors_cond)} labels) ===")
            labels, obs, obs_mat, null_means, null_med_mat = run_test(
                anchors_cond,
                expand_fn,
                network.get(cond, {}),
                pool_arr,
                args.n_perm,
                rng,
            )
            p = (np.sum(null_means >= obs) + 1) / (args.n_perm + 1)
            print(f"  observed mean J = {obs:.4f}")
            print(
                f"  null    mean J = {np.mean(null_means):.4f} "
                f"± {np.std(null_means):.4f}"
            )
            print(f"  empirical p    = {p:.4g}")
            rows.append(
                {
                    "condition": cond,
                    "expansion": name,
                    "obs_mean": obs,
                    "null_mean": float(np.mean(null_means)),
                    "null_std": float(np.std(null_means)),
                    "p_value": p,
                    "n_labels": len(labels),
                }
            )

            plot_null_hist(axes_h[r, c], obs, null_means, p, f"{cond} | {name}")
            # plot_excess(axes_e[r, c], obs_mat, null_med_mat, labels, f"{cond} | {name}")
            plot_excess_clustered(
                axes_e[r, c], obs_mat, null_med_mat, labels, f"{cond} | {name}"
            )

    fig_h.suptitle(
        f"Permutation null vs observed mean pairwise Jaccard "
        f"({args.n_perm} perms, anchor = {args.anchor_annot})",
        y=1.02,
        fontsize=12,
    )
    fig_h.tight_layout()
    null_path = out_prefix.with_name(out_prefix.name + "_null_distribution.png")
    fig_h.savefig(null_path, dpi=150, bbox_inches="tight")
    plt.close(fig_h)
    print(f"\nWrote {null_path}")

    fig_e.suptitle(
        f"Per-pair overlap excess (J_obs - median J_null), "
        f"anchor = {args.anchor_annot}",
        y=1.02,
        fontsize=12,
    )
    fig_e.tight_layout()
    excess_path = out_prefix.with_name(
        out_prefix.name + "_excess_heatmap_clustered.png"
    )
    fig_e.savefig(excess_path, dpi=150, bbox_inches="tight")
    plt.close(fig_e)
    print(f"Wrote {excess_path}")

    hub_path = out_prefix.with_name(out_prefix.name + "_hub_diagnostic.png")
    hub_diagnostic(cis_to_trans, trans_to_cis, args.conditions, hub_path)

    summary_path = out_prefix.with_name(out_prefix.name + "_summary.csv")
    pl.DataFrame(rows).write_csv(summary_path)
    print(f"Wrote {summary_path}")

    print("\n=== Summary ===")
    print(
        f"{'cond':<10s} {'expansion':<12s} {'n':>3s} "
        f"{'obs':>8s} {'null':>8s} {'p':>10s}"
    )
    for r in rows:
        print(
            f"{r['condition']:<10s} {r['expansion']:<12s} "
            f"{r['n_labels']:>3d} {r['obs_mean']:>8.4f} "
            f"{r['null_mean']:>8.4f} {r['p_value']:>10.4g}"
        )


if __name__ == "__main__":
    main()
