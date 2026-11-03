#!/usr/bin/env python
"""Fig. S4 - gene-set sizes per disease: anchors, then one panel per condition.

Coloc genes (anchors) are one set per disease and do not depend on the
condition, and they are 1-2 orders of magnitude smaller than the expanded
sets. They get their own leading panel on its own linear axis, labelled
with counts; the condition panels stack targets + regulators on a shared
linear axis. Nothing is log-scaled, so bar lengths stay comparable within
each panel. IMDs and non-IMD controls are separate, shaded blocks; within
each block diseases are sorted by coloc-gene count.
"""

import argparse

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

from common import (ANNOT_COLORS, ANNOT_LABEL, ANNOTS, COND_LABEL, CONDS,
                    add_variant_args, color_ticklabels, disease_order, load_diseases,
                    load_gene_sets, name, plain_log_axis, resolve_variant, save,
                    shade_groups, size_table)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log-x", action="store_true",
                    help="Log x axis on the condition panels (stack segments "
                         "then no longer read as proportions; use for a "
                         "dynamic-range overview only).")
    add_variant_args(ap)
    args = resolve_variant(ap.parse_args())

    diseases = load_diseases()
    sizes = size_table(load_gene_sets(args.gene_sets_dir), diseases)
    wide = sizes.pivot_table(index=["study", "condition"], columns="annot",
                             values="n_genes", fill_value=0)
    anchors = wide.xs("Rest", level="condition")["H4_genes"]
    order = disease_order(diseases, anchors)[::-1]  # bottom -> top
    y = range(len(order))

    fig, axes = plt.subplots(1, len(CONDS) + 1, figsize=(15, 0.32 * len(order) + 1.6),
                             sharey=True,
                             gridspec_kw={"width_ratios": [0.45] + [1] * len(CONDS)})
    ax0, cond_axes = axes[0], axes[1:]
    for ax in cond_axes[1:]:
        ax.sharex(cond_axes[0])

    # Anchors: condition-independent, own axis, direct value labels.
    shade_groups(ax0, order, diseases, label=False)
    a = [int(anchors[s]) for s in order]
    ax0.barh(list(y), a, height=0.7, color=ANNOT_COLORS["H4_genes"],
             edgecolor="white", linewidth=0.8, zorder=2)
    for yi, v in zip(y, a):
        ax0.text(v + max(a) * 0.03, yi, str(v), va="center", fontsize=8,
                 color="#52514e")
    ax0.set_xlim(0, max(a) * 1.25)
    ax0.set_title("Coloc genes\n(all conditions)", fontsize=11, loc="left")
    ax0.set_xlabel("Number of genes")

    for ax, cond in zip(cond_axes, CONDS):
        shade_groups(ax, order, diseases, label=(cond == CONDS[-1]))
        left = [0] * len(order)
        for annot in ANNOTS[1:]:
            w = [wide.loc[(s, cond), annot] if (s, cond) in wide.index else 0
                 for s in order]
            ax.barh(list(y), w, left=left, height=0.7, color=ANNOT_COLORS[annot],
                    edgecolor="white", linewidth=0.8, zorder=2)
            left = [a + b for a, b in zip(left, w)]
        ax.set_title(f"Expanded sets\n{COND_LABEL[cond]}", fontsize=11, loc="left")
        ax.set_xlabel("Number of genes")
        if args.log_x:
            ax.set_xscale("log")
    if args.log_x:
        for ax in cond_axes:
            plain_log_axis(ax, "x")
    for ax in axes:
        ax.grid(axis="x", zorder=1)
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", length=0)

    axes[0].set_yticks(list(y))
    axes[0].set_yticklabels([name(diseases, s) for s in order])
    color_ticklabels(axes[0].get_yticklabels(), order, diseases)
    axes[0].set_ylim(-0.5, len(order) - 0.5)

    handles = [mpatches.Patch(color=ANNOT_COLORS[a], label=ANNOT_LABEL[a]) for a in ANNOTS]
    fig.legend(handles=handles, loc="upper center", ncol=3,
               bbox_to_anchor=(0.5, 1.03))
    fig.tight_layout()
    save(fig, "supp_gene_set_sizes_by_condition" + ("_logx" if args.log_x else ""),
         args.out_dir)


if __name__ == "__main__":
    main()
