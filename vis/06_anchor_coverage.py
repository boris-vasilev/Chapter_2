#!/usr/bin/env python
"""Figure 2.1B - which anchors can be expanded, per disease and condition.

No list of measured or perturbed genes is available locally (the DE h5ad
lives on the cluster), so anchors are classified from the gene-set edge
list itself (common.load_edges):

  * both              - >= 1 incoming edge and >= 1 outgoing edge
  * regulators only   - >= 1 incoming edge (some knockdown changes it),
                        no outgoing edge (not a strong perturbation)
  * targets only      - >= 1 outgoing edge, no incoming edge
  * neither           - no edge in that condition

Horizontal stacked bars, one facet per condition, same disease order and
IMD / control shading as the other figures. Diseases with no targets in
any condition are marked with an asterisk; the asterisk is explained in the
caption (NO_TARGETS_NOTE), not on the figure. Also writes anchor_coverage.csv.

draw_anchor_coverage() draws the facets onto given Axes, so 07_fig_2_1.py
can place them as Fig. 2.1B.
"""

import argparse
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import polars as pl

from common import (ANNOT_COLORS, COND_LABEL, CONDS, add_variant_args,
                    color_ticklabels, disease_order, load_diseases,
                    load_edges, load_gene_sets, name, resolve_variant, save,
                    shade_groups)

CATS = ["both", "regulators only", "targets only", "neither"]
CAT_LABEL = {
    "both": "Both (≥1 incoming and ≥1 outgoing edge)",
    "regulators only": "Regulators only (≥1 incoming edge)",
    "targets only": "Targets only (≥1 outgoing edge)",
    "neither": "Neither (no edge)",
}
CAT_COLORS = {
    # Yellow: not used for anchors / regulators / targets anywhere else;
    # validated against its stack neighbour (aqua), CVD dE 9.1.
    "both": "#eda100",
    "regulators only": ANNOT_COLORS["regulators_H4"],
    "targets only": ANNOT_COLORS["targets_H4"],
    "neither": "#d6d5d0",
}
NO_TARGETS_NOTE = "* no targets in any condition (no anchor has an outgoing edge)"


def classify(edges, sets, studies):
    rows = []
    for cond in CONDS:
        e = edges.filter(pl.col("culture_condition") == cond)
        has_in, has_out = set(e["trans_ensg"]), set(e["cis_ensg"])
        for study in studies:
            for a in sets[(study, cond, "H4_genes")]:
                i, o = a in has_in, a in has_out
                cat = ("both" if i and o else "regulators only" if i
                       else "targets only" if o else "neither")
                rows.append({"study": study, "condition": cond, "anchor": a,
                             "category": cat})
    return pd.DataFrame(rows).sort_values(["study", "condition", "anchor"],
                                          ignore_index=True)


def draw_anchor_coverage(axes, sets, edges, diseases, label_groups=True):
    """Draw one stacked-bar facet per condition onto ``axes`` (len 3).

    Returns (classification table, legend handles).
    """
    cls = classify(edges, sets, list(diseases.index))
    counts = (cls.groupby(["study", "condition", "category"]).size()
              .unstack(fill_value=0).reindex(columns=CATS, fill_value=0))
    anchors = {s: len(sets[(s, "Rest", "H4_genes")]) for s in diseases.index}
    order = disease_order(diseases, anchors)[::-1]
    no_targets = {s for s in diseases.index
                  if all(not sets[(s, c, "targets_H4")] for c in CONDS)}
    y = list(range(len(order)))

    for k, (ax, cond) in enumerate(zip(axes, CONDS)):
        shade_groups(ax, order, diseases, label=label_groups and cond == CONDS[-1])
        left = [0] * len(order)
        for cat in CATS:
            w = [counts.loc[(s, cond), cat] if (s, cond) in counts.index else 0
                 for s in order]
            ax.barh(y, w, left=left, height=0.7, color=CAT_COLORS[cat],
                    edgecolor="white", linewidth=0.8, zorder=2)
            left = [a + b for a, b in zip(left, w)]
        ax.set_title(COND_LABEL[cond], loc="left")
        ax.set_xlabel("Anchors (coloc genes)")
        ax.grid(axis="x", zorder=1)
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", length=0)
        ax.set_ylim(-0.5, len(order) - 0.5)
        ax.set_yticks(y)
        if k == 0:
            ax.set_yticklabels([name(diseases, s) + (" *" if s in no_targets else "")
                                for s in order])
            color_ticklabels(ax.get_yticklabels(), order, diseases)
        else:
            ax.tick_params(axis="y", labelleft=False)
    xmax = max(ax.get_xlim()[1] for ax in axes)
    for ax in axes:
        ax.set_xlim(0, xmax)
    handles = [mpatches.Patch(color=CAT_COLORS[c], label=CAT_LABEL[c]) for c in CATS]
    return cls, counts, order, no_targets, handles


def main():
    ap = argparse.ArgumentParser()
    add_variant_args(ap)
    args = resolve_variant(ap.parse_args())

    diseases = load_diseases()
    sets = load_gene_sets(args.gene_sets_dir)
    edges = load_edges(min_abs_logfc=args.min_abs_logfc)

    fig, axes = plt.subplots(1, len(CONDS), figsize=(13, 0.32 * len(diseases) + 1.8))
    cls, counts, order, no_targets, handles = draw_anchor_coverage(
        axes, sets, edges, diseases)
    cls.to_csv(Path(args.out_dir) / "anchor_coverage.csv", index=False)

    fig.legend(handles=handles, loc="upper center", ncol=2, fontsize=9,
               bbox_to_anchor=(0.5, 1.07))
    fig.tight_layout()
    save(fig, "anchor_coverage", args.out_dir)

    print(counts.groupby(level="condition").sum().loc[CONDS].to_string())
    print("No targets in any condition:",
          ", ".join(name(diseases, s) for s in order[::-1] if s in no_targets))


if __name__ == "__main__":
    main()
