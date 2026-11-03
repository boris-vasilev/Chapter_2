#!/usr/bin/env python
"""Compact gene-count stacked barplot per (study × condition × annot).

Each study gets a *group* of three thin stacked bars (Rest | 8h | 48h,
left to right). Each bar stacks the three H4 annotations:

    H4_genes (blue, bottom)  targets_H4 (orange)  regulators_H4 (green, top)

Single axis, single legend — fits the top or bottom half of a 16:9 slide.
"""

import argparse
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

CONDS = ["Rest", "Stim8hr", "Stim48hr"]
COND_SHORT = {"Rest": "R", "Stim8hr": "8h", "Stim48hr": "48h"}

DEFAULT_ANNOTS = ["H4_genes", "targets_H4", "regulators_H4"]
DEFAULT_COLORS = {
    "H4_genes": "#1f77b4",
    "targets_H4": "#ff7f0e",
    "regulators_H4": "#2ca02c",
    "l2g_genes": "#1f77b4",
    "targets_l2g": "#ff7f0e",
    "regulators_l2g": "#2ca02c",
}
DEFAULT_DISPLAY = {
    "H4_genes": "Coloc genes",
    "targets_H4": "Targets",
    "regulators_H4": "Regulators",
    "l2g_genes": "L2G genes",
    "targets_l2g": "Targets",
    "regulators_l2g": "Regulators",
}


def count_genes(path):
    with open(path) as f:
        return sum(1 for line in f if line.strip())


def discover(gene_sets_dir, annots):
    """dict[(study, cond)][annot] = gene count."""
    annot_set = set(annots)
    annot_re = "|".join(re.escape(a) for a in sorted(annots, key=len, reverse=True))
    pat = re.compile(rf"^(?P<study>.+)_(?P<cond>[^_]+)_(?P<annot>{annot_re})\.txt$")
    counts = defaultdict(dict)
    for p in sorted(Path(gene_sets_dir).glob("*.txt")):
        m = pat.match(p.name)
        if not m or m.group("annot") not in annot_set:
            continue
        if m.group("cond") not in CONDS:
            continue
        counts[(m.group("study"), m.group("cond"))][m.group("annot")] = count_genes(p)
    return counts


def plot(counts, annots, colors, display_names, out_path, figsize, log_y):
    studies = sorted({s for (s, c) in counts.keys()})
    labels = [s.split("_")[0] for s in studies]
    n = len(studies)

    fig, ax = plt.subplots(figsize=figsize)

    bar_w = 0.28
    spacing = 1.4  # >1 = extra gap between disease groups
    offsets = {c: (i - 1) * (bar_w + 0.015) for i, c in enumerate(CONDS)}
    centers = [s_idx * spacing for s_idx in range(n)]

    for s_idx, study in enumerate(studies):
        cx = centers[s_idx]
        for cond in CONDS:
            x = cx + offsets[cond]
            bottom = 0
            for annot in annots:
                h = counts.get((study, cond), {}).get(annot, 0)
                if h <= 0:
                    continue
                ax.bar(x, h, width=bar_w, bottom=bottom,
                       color=colors[annot], edgecolor="none")
                bottom += h

    ax.set_xticks(centers)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=15)
    ax.tick_params(axis="x", length=0)
    ax.tick_params(axis="y", labelsize=13)
    ax.set_ylabel("Number of genes", fontsize=16)
    if log_y:
        ax.set_yscale("log")

    ax.set_xlim(-spacing * 0.55, centers[-1] + spacing * 0.55)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    handles = [mpatches.Patch(color=colors[a], label=display_names[a])
               for a in annots]
    ax.legend(handles=handles, loc="upper right", frameon=False,
              fontsize=13, ncol=3, handlelength=1.4, columnspacing=1.2,
              bbox_to_anchor=(1.0, 1.12))

    ax.text(0.0, 1.04,
            f"Bars per disease (L→R): {' | '.join(COND_SHORT[c] for c in CONDS)}",
            transform=ax.transAxes, ha="left", va="bottom",
            fontsize=12, style="italic", color="#555555")

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-sets-dir", required=True,
                        help="e.g. data/gene_sets/coloc")
    parser.add_argument("--annots", nargs=3, default=DEFAULT_ANNOTS,
                        metavar=("ANCHOR", "TARGETS", "REGULATORS"),
                        help="Three annot names matching the gene-set "
                             "filenames (stack order = bottom→top). "
                             "Default: H4_genes targets_H4 regulators_H4")
    parser.add_argument("--figsize", nargs=2, type=float, default=(13.0, 5.2),
                        metavar=("W", "H"),
                        help="Figure size in inches (default 13×5.2 for half "
                             "of a 16:9 slide).")
    parser.add_argument("--legend-labels", nargs=3, default=None,
                        metavar=("ANCHOR", "TARGETS", "REGULATORS"),
                        help="Override the legend labels (in --annots order). "
                             "Default: 'Coloc genes' 'Targets' 'Regulators' "
                             "for H4, 'L2G genes' 'Targets' 'Regulators' for "
                             "L2G.")
    parser.add_argument("--log-y", action="store_true",
                        help="Use log scale on y (compresses dynamic range so "
                             "small annotations stay visible).")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    fallback = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    colors = {a: DEFAULT_COLORS.get(a, fallback[i])
              for i, a in enumerate(args.annots)}

    if args.legend_labels:
        display = dict(zip(args.annots, args.legend_labels))
    else:
        display = {a: DEFAULT_DISPLAY.get(a, a) for a in args.annots}

    counts = discover(args.gene_sets_dir, args.annots)
    if not counts:
        raise SystemExit("No (study, condition) gene-set files found.")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plot(counts, args.annots, colors, display, out_path,
         tuple(args.figsize), args.log_y)


if __name__ == "__main__":
    main()
