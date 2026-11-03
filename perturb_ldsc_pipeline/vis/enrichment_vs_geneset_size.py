#!/usr/bin/env python
"""Scatter: log2(Enrichment) vs gene-set size, one panel per annotation.

Joins the long-format LDSC summary (parse_ldsc_results.py) with the raw
gene-set files (one gene ID per line) under data/gene_sets/<module>/.
Plots three side-by-side panels — Coloc genes / Regulators / Targets —
with every (disease × condition) as one point.

Use it to sanity-check whether the bright cells in plot_enrichment_panels
are an artefact of gene-set size: if all the FDR-significant points are
also the smallest sets, the enrichment is likely a small-set artefact.

  Y = log2(Enrichment), clipped at log2_vmin (default 0 = depletion floor)
  X = log10(n_genes)
  colour = condition (Rest / 8h / 48h)
  filled + black outline + label = BH-FDR < --fdr across all visible cells
"""

import argparse
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

CONDS = ["Rest", "Stim8hr", "Stim48hr"]
COND_COLORS = {"Rest": "#1f77b4", "Stim8hr": "#ff7f0e", "Stim48hr": "#2ca02c"}


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan)
    mask = ~np.isnan(p)
    n = mask.sum()
    if n == 0:
        return out
    idx = np.where(mask)[0]
    order = np.argsort(p[idx])
    ranked = p[idx][order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out[idx[order]] = adj
    return out


def spearman_rho_p(x, y):
    """Spearman ρ + normal-approx two-sided p. NaN-tolerant."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = ~(np.isnan(x) | np.isnan(y))
    x, y = x[m], y[m]
    n = len(x)
    if n < 3:
        return float("nan"), float("nan")
    rx = np.empty(n)
    ry = np.empty(n)
    rx[np.argsort(x)] = np.arange(n) + 1
    ry[np.argsort(y)] = np.arange(n) + 1
    if rx.std() == 0 or ry.std() == 0:
        return 0.0, 1.0
    rho = float(np.corrcoef(rx, ry)[0, 1])
    z = rho * math.sqrt(n - 1)
    p = math.erfc(abs(z) / math.sqrt(2))
    return rho, p


def count_genes(path):
    with open(path) as f:
        return sum(1 for line in f if line.strip())


def discover_sizes(gene_sets_dir, annots):
    """{(study, cond, annot): n_genes}"""
    annot_set = set(annots)
    annot_re = "|".join(re.escape(a) for a in sorted(annots, key=len, reverse=True))
    pat = re.compile(rf"^(?P<study>.+)_(?P<cond>[^_]+)_(?P<annot>{annot_re})\.txt$")
    sizes = {}
    for p in sorted(Path(gene_sets_dir).glob("*.txt")):
        m = pat.match(p.name)
        if not m or m.group("annot") not in annot_set:
            continue
        if m.group("cond") not in CONDS:
            continue
        sizes[(m.group("study"), m.group("cond"), m.group("annot"))] = count_genes(p)
    return sizes


def plot(rows, annots, display, out_path, figsize, fdr_thr,
         log2_vmin, log2_vmax, annotate_sig):
    log2_floor = 2 ** log2_vmin

    # BH-FDR jointly across all visible cells
    p = np.array([r["p"] for r in rows], dtype=float)
    fdr = bh_fdr(p)
    for r, q in zip(rows, fdr):
        r["fdr"] = float(q)

    fig, axes = plt.subplots(1, len(annots), figsize=figsize, sharey=True)
    if len(annots) == 1:
        axes = [axes]

    for ax, annot in zip(axes, annots):
        sub = [r for r in rows if r["annot"] == annot]
        if not sub:
            ax.set_title(f"{display[annot]} (no data)", fontsize=14)
            continue

        xs = np.array([math.log10(max(r["size"], 1)) for r in sub])
        ys = np.array([math.log2(max(r["enrichment"], log2_floor)) for r in sub])
        sig = np.array([r["fdr"] < fdr_thr for r in sub])
        colours = [COND_COLORS[r["condition"]] for r in sub]

        # Non-significant points: small, semi-transparent, no edge
        ns = ~sig
        if ns.any():
            ax.scatter(xs[ns], ys[ns],
                       c=[colours[i] for i in np.where(ns)[0]],
                       s=60, alpha=0.45, edgecolor="none", zorder=2)
        # Significant points: bigger, opaque, black outline
        if sig.any():
            ax.scatter(xs[sig], ys[sig],
                       c=[colours[i] for i in np.where(sig)[0]],
                       s=130, alpha=0.95, edgecolor="black",
                       linewidth=1.0, zorder=3)

        if annotate_sig:
            for r, x, y, is_sig in zip(sub, xs, ys, sig):
                if not is_sig:
                    continue
                ax.annotate(r["label"], (x, y),
                            xytext=(6, 0), textcoords="offset points",
                            fontsize=10, va="center")

        rho, p_corr = spearman_rho_p(xs, ys)
        ax.axhline(0, color="grey", linewidth=0.8, linestyle="--", alpha=0.7)
        ax.set_title(
            f"{display[annot]}\nSpearman ρ = {rho:.2f}   p = {p_corr:.2g}   (n = {len(sub)})",
            fontsize=14,
        )
        ax.set_xlabel("log10(gene-set size)", fontsize=15)
        ax.tick_params(labelsize=12)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[0].set_ylabel("log2(Enrichment)", fontsize=15)
    if log2_vmax is not None:
        for ax in axes:
            ax.set_ylim(top=log2_vmax + 0.3)

    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="",
                   markerfacecolor=COND_COLORS[c], markeredgecolor="none",
                   markersize=9, label=c)
        for c in CONDS
    ]
    handles.append(plt.Line2D([0], [0], marker="o", linestyle="",
                              markerfacecolor="lightgrey",
                              markeredgecolor="black", markeredgewidth=1.0,
                              markersize=11, label=f"FDR < {fdr_thr}"))
    fig.legend(handles=handles, loc="upper center", frameon=False,
               fontsize=12, ncol=4, bbox_to_anchor=(0.5, 1.04))

    fig.suptitle("LDSC enrichment vs gene-set size",
                 fontsize=17, fontweight="bold", y=1.10)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", required=True,
                        help="CSV from parse_ldsc_results.py")
    parser.add_argument("--gene-sets-dir", required=True,
                        help="e.g. data/gene_sets/coloc")
    parser.add_argument("--anchor-annot", default="H4_genes")
    parser.add_argument("--regulators-annot", default="regulators_H4")
    parser.add_argument("--targets-annot", default="targets_H4")
    parser.add_argument("--anchor-label", default="Coloc genes")
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--log2-vmin", type=float, default=0.0,
                        help="Floor for log2(Enrichment). "
                             "Default 0 collapses depletion (default 0).")
    parser.add_argument("--log2-vmax", type=float, default=None,
                        help="Optional cap on the y-axis.")
    parser.add_argument("--no-labels", action="store_true",
                        help="Don't label significant points.")
    parser.add_argument("--label-col", default="trait_short",
                        help="Column to use for point labels "
                             "(default trait_short).")
    parser.add_argument("--figsize", nargs=2, type=float, default=(16.0, 5.5),
                        metavar=("W", "H"))
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    annots = [args.anchor_annot, args.regulators_annot, args.targets_annot]
    display = {
        args.anchor_annot: args.anchor_label,
        args.regulators_annot: "Regulators",
        args.targets_annot: "Targets",
    }

    df = pl.read_csv(args.summary).filter(pl.col("annot").is_in(annots))
    if df.height == 0:
        raise SystemExit(f"No rows matched annots {annots}")

    sizes = discover_sizes(args.gene_sets_dir, annots)

    rows = []
    for r in df.iter_rows(named=True):
        size = sizes.get((r["study"], r["condition"], r["annot"]))
        if not size:
            continue
        label = r.get(args.label_col) or r["study"].split("_")[0]
        rows.append({
            "study": r["study"],
            "condition": r["condition"],
            "annot": r["annot"],
            "enrichment": float(r["enrichment"]),
            "p": float(r["enrichment_p"]),
            "size": int(size),
            "label": label,
        })

    if not rows:
        raise SystemExit("No matched (enrichment × size) rows found")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plot(
        rows, annots, display, out_path,
        figsize=tuple(args.figsize),
        fdr_thr=args.fdr,
        log2_vmin=args.log2_vmin,
        log2_vmax=args.log2_vmax,
        annotate_sig=not args.no_labels,
    )


if __name__ == "__main__":
    main()
