#!/usr/bin/env python
"""Disease x programme S-LDSC heritability-enrichment heatmap.

Reads the tidy CSV written by collect_program_results.py and draws, one
figure per condition, a heatmap of the cell-type-specific S-LDSC
coefficient z-score for every (programme, disease) pair:

    rows    = Perturb-seq programmes (labelled by manual_annotation)
    columns = diseases (GWAS studies), immune-mediated grouped left
    colour  = Coefficient_z-score (RdBu_r, centred at 0)
    * / ** / *** = nominal one-sided p (P(Z>z)) < 0.05 / 0.01 / 0.001
    gold border  = BH-FDR < 0.05 (the global FDR from the collector)
    ·            = (programme, disease) not tested

A positive z means the programme's genes contribute more per-SNP
heritability than expected given the baseline model — i.e. that perturb-seq
module is enriched for the disease's polygenic signal.
"""

import argparse
from collections import defaultdict
from pathlib import Path

import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from matplotlib.patches import Rectangle

COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]

# Same immune-mediated grouping as approved_drug_enrichment_all.py, keyed on
# the short trait code (study.split("_")[0]). Used only for column ordering
# and the divider; unknown traits fall to the right.
IMD_DISEASES = [
    "AST", "ASTAO", "ASTCO", "ASPH", "ATD", "AITD",
    "CD", "IBD", "UC", "RA", "SLE", "T1D",
]
SHORT_NAMES = {
    "AST": "Asthma", "ASTAO": "Asthma (AO)", "ASTCO": "Asthma (CO)",
    "ATD": "Atopic derm.", "AITD": "Autoimm. thyroid", "CD": "Crohn's",
    "IBD": "IBD", "UC": "Ulcer. colitis", "RA": "Rheum. arthritis",
    "SLE": "SLE", "T1D": "T1D", "IS": "Ischemic stroke",
    "MYOI": "Myocard. infarct.", "OA": "Osteoarthritis", "T2D": "T2D",
    "ADFH": "Alzheimer's", "SCZ": "Schizophrenia",
}


def order_diseases(studies):
    """IMD-first, then alphabetical; returns (ordered_studies, n_imd)."""
    def short(s):
        return s.split("_")[0]
    imd = [s for s in studies if short(s) in IMD_DISEASES]
    non = [s for s in studies if short(s) not in IMD_DISEASES]
    imd.sort(key=lambda s: (IMD_DISEASES.index(short(s)), s))
    non.sort()
    return imd + non, len(imd)


def disease_labels(studies):
    """Short, de-duplicated column labels (two RA studies -> RA, RA#2)."""
    counts = defaultdict(int)
    labels = {}
    for s in studies:
        code = s.split("_")[0]
        base = SHORT_NAMES.get(code, code)
        counts[base] += 1
        labels[s] = base if counts[base] == 1 else f"{base}#{counts[base]}"
    return labels


def program_labels(rows_for_cond):
    """Programme row labels: manual_annotation, deduped by programme number."""
    ann_counts = defaultdict(int)
    for r in rows_for_cond:
        ann_counts[r["manual_annotation"]] += 1
    labels = {}
    for r in rows_for_cond:
        cl = r["cluster"]
        ann = r["manual_annotation"] or f"cluster {cl}"
        ann = ann if len(ann) <= 30 else ann[:29] + "…"
        labels[cl] = (f"{ann} (P{r['program']})"
                      if ann_counts[r["manual_annotation"]] > 1 else ann)
    return labels


def plot_condition(df_cond, cond, out_path, vmax, fdr_thr, only_sig):
    studies = df_cond.get_column("study").unique().to_list()
    studies, n_imd = order_diseases(studies)
    dlabels = disease_labels(studies)
    study_idx = {s: j for j, s in enumerate(studies)}

    # Programme rows: one per cluster present in this condition.
    cl_meta = {}
    for r in df_cond.iter_rows(named=True):
        cl_meta.setdefault(r["cluster"], r)
    clusters = list(cl_meta)

    # Optionally drop programmes never significant (nominal) in any disease.
    if only_sig:
        sig_cl = {
            r["cluster"] for r in df_cond.iter_rows(named=True)
            if r["coef_p_onesided"] < 0.05
        }
        clusters = [c for c in clusters if c in sig_cl]
        if not clusters:
            print(f"  [{cond}] no programme nominally significant; skipping")
            return

    plabels_map = program_labels([cl_meta[c] for c in clusters])
    # Order programmes by peak z-score across diseases (strongest on top).
    peak = defaultdict(lambda: -np.inf)
    for r in df_cond.iter_rows(named=True):
        if r["cluster"] in plabels_map:
            peak[r["cluster"]] = max(peak[r["cluster"]], r["coefficient_z"])
    clusters.sort(key=lambda c: peak[c], reverse=True)
    cl_idx = {c: i for i, c in enumerate(clusters)}

    shape = (len(clusters), len(studies))
    z = np.full(shape, np.nan)
    pmat = np.full(shape, np.nan)
    fmat = np.full(shape, np.nan)
    for r in df_cond.iter_rows(named=True):
        if r["cluster"] not in cl_idx or r["study"] not in study_idx:
            continue
        i, j = cl_idx[r["cluster"]], study_idx[r["study"]]
        z[i, j] = r["coefficient_z"]
        pmat[i, j] = r["coef_p_onesided"]
        fmat[i, j] = r["coef_fdr_bh"]

    if vmax is None:
        finite = np.abs(z[np.isfinite(z)])
        vmax = float(np.percentile(finite, 99)) if finite.size else 4.0
        vmax = max(vmax, 1.0)

    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("white")

    fig_w = max(8, 0.55 * len(studies) + 4)
    fig_h = max(5, 0.42 * len(clusters) + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(z, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")

    for i in range(len(clusters)):
        for j in range(len(studies)):
            p = pmat[i, j]
            if np.isnan(p):
                ax.text(j, i, "·", ha="center", va="center",
                        color="grey", fontsize=11)
                continue
            stars = ("***" if p < 0.001 else "**" if p < 0.01
                     else "*" if p < 0.05 else "")
            if stars:
                t = ax.text(j, i, stars, ha="center", va="center",
                            color="white", fontsize=12, fontweight="bold")
                t.set_path_effects([path_effects.Stroke(linewidth=2,
                                    foreground="black"), path_effects.Normal()])
            if not np.isnan(fmat[i, j]) and fmat[i, j] < fdr_thr:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1.0, 1.0,
                             fill=False, edgecolor="#FFD700", linewidth=3.0))

    if 0 < n_imd < len(studies):
        ax.axvline(n_imd - 0.5, color="black", linewidth=2.0)
        ax.text((n_imd - 1) / 2, -0.9, "Immune-mediated", ha="center",
                fontsize=12, fontweight="bold")
        ax.text((n_imd + len(studies) - 1) / 2, -0.9, "Other", ha="center",
                fontsize=12, fontweight="bold")

    ax.set_xticks(range(len(studies)))
    ax.set_xticklabels([dlabels[s] for s in studies], rotation=45,
                       ha="right", fontsize=11)
    ax.set_yticks(range(len(clusters)))
    ax.set_yticklabels([plabels_map[c] for c in clusters], fontsize=11)
    ax.tick_params(axis="both", length=0)
    ax.set_xlim(-0.5, len(studies) - 0.5)
    ax.set_ylim(len(clusters) - 0.5, -0.5)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("S-LDSC coefficient z-score", fontsize=12)
    cbar.ax.tick_params(labelsize=10)

    ax.set_title(
        f"Programme heritability enrichment — {cond}\n"
        f"* p<0.05  ** p<0.01  *** p<0.001 (one-sided);  "
        f"gold border: BH-FDR < {fdr_thr};  · = not tested",
        fontsize=13, fontweight="bold", pad=28,
    )
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results", required=True,
                   help="program_enrichment.csv from collect_program_results.py")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--conditions", nargs="+", default=COND_ORDER)
    p.add_argument("--vmax", type=float, default=None,
                   help="Symmetric z colour cap (default = 99th pct of |z|).")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--only-significant", action="store_true",
                   help="Drop programmes never nominally significant (p<0.05).")
    p.add_argument("--format", default="png", choices=["png", "pdf"])
    args = p.parse_args()

    df = pl.read_csv(args.results)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for cond in args.conditions:
        df_cond = df.filter(pl.col("condition") == cond)
        if df_cond.height == 0:
            print(f"  [{cond}] no rows; skipping")
            continue
        out_path = out_dir / f"program_heritability_{cond}.{args.format}"
        plot_condition(df_cond, cond, out_path, args.vmax, args.fdr,
                       args.only_significant)


if __name__ == "__main__":
    main()
