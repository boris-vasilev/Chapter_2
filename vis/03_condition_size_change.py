#!/usr/bin/env python
"""Fig. S3 - how regulator / target set sizes change with stimulation.

Top row: set size (log scale) per disease across Rest -> 8 h -> 48 h, one
line per disease, coloured and dashed by IMD / non-IMD, with the
across-disease median overlaid. Bottom row: log2 fold change relative to
Rest, so diseases with very different anchor counts are comparable.

Per annotation, diseases with an empty set in any condition are excluded
(their fold change is undefined); the panel titles give n and the excluded
diseases. Friedman / Wilcoxon results are not drawn on the figure; they are
written to supp_set_size_by_condition_tests.csv with a note that the
diseases are not independent observations.
"""

import argparse
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import friedmanchisquare, wilcoxon

from common import (ANNOT_LABEL, COND_LABEL, CONDS, GROUP_COLORS,
                    GROUP_LABEL, GROUP_LINESTYLE, GROUPS, load_diseases,
                    load_gene_sets, name, plain_log_axis, save, size_table,
                    add_variant_args, resolve_variant)

ANNOTS = ["regulators_H4", "targets_H4"]
TEST_NOTE = "diseases not independent (shared anchors, duplicated traits)"


def condition_tests(wide):
    """wide: rows = studies, columns = CONDS, no empty sets."""
    out = {"n_diseases": len(wide),
           "friedman_p": friedmanchisquare(*[wide[c] for c in CONDS]).pvalue}
    for c in CONDS[1:]:
        out[f"wilcoxon_{c}_vs_Rest_p"] = wilcoxon(wide[c], wide["Rest"]).pvalue
        out[f"median_log2fc_{c}"] = np.median(np.log2(wide[c] / wide["Rest"]))
    return out


def main():
    ap = argparse.ArgumentParser()
    add_variant_args(ap)
    args = resolve_variant(ap.parse_args())

    diseases = load_diseases()
    sizes = size_table(load_gene_sets(args.gene_sets_dir), diseases)
    x = np.arange(len(CONDS))

    fig, axes = plt.subplots(2, 2, figsize=(9, 7.5), sharex=True)
    stats, notes = [], []
    for j, annot in enumerate(ANNOTS):
        wide = sizes[sizes.annot == annot].pivot(index="study", columns="condition",
                                                 values="n_genes")[CONDS]
        empty = wide.index[(wide == 0).any(axis=1)]
        wide = wide.drop(empty)
        excluded = sorted(name(diseases, s) for s in empty)
        lfc = np.log2(wide.div(wide["Rest"], axis=0))
        for data, ax in [(wide, axes[0, j]), (lfc, axes[1, j])]:
            for study, row in data.iterrows():
                g = diseases.loc[study, "group"]
                ax.plot(x, row.values, color=GROUP_COLORS[g], ls=GROUP_LINESTYLE[g],
                        lw=1.2, alpha=0.55, marker="o", ms=3.5, zorder=2)
            ax.plot(x, data.median().values, color="#0b0b0b", lw=2.5, marker="o",
                    ms=7, zorder=3)
            ax.grid(axis="y")
            ax.set_axisbelow(True)
        axes[0, j].set_yscale("log")
        plain_log_axis(axes[0, j], "y")
        axes[1, j].axhline(0, color="#52514e", lw=0.8, zorder=1)

        axes[0, j].set_title(f"{ANNOT_LABEL[annot]} (n = {len(wide)} diseases)",
                             loc="left", fontsize=11)
        if excluded:
            notes.append(f"{ANNOT_LABEL[annot]}: excluded (empty set in ≥1 "
                         f"condition): " + ", ".join(excluded) + ".")

        stats.append({"annot": annot, **condition_tests(wide),
                      "excluded": "; ".join(excluded), "note": TEST_NOTE})

    axes[0, 0].set_ylabel("Genes in set (log scale)")
    axes[1, 0].set_ylabel("log$_2$ fold change vs Rest")
    for ax in axes[1]:
        ax.set_xticks(x)
        ax.set_xticklabels([COND_LABEL[c] for c in CONDS])
        ax.set_xlim(-0.25, len(CONDS) - 0.75)

    handles = [Line2D([], [], color=GROUP_COLORS[g], ls=GROUP_LINESTYLE[g],
                      marker="o", ms=3.5, label=GROUP_LABEL[g]) for g in GROUPS]
    handles.append(Line2D([], [], color="#0b0b0b", lw=2.5, marker="o", ms=7,
                          label="Median across diseases"))
    fig.legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.03))
    fig.tight_layout()
    if notes:
        fig.text(0.01, -0.01, "\n".join(textwrap.fill(t, 110) for t in notes),
                 va="top", fontsize=8.5, color="#52514e")
    save(fig, "supp_set_size_by_condition", args.out_dir)

    stats = pd.DataFrame(stats)
    print(stats.to_string(index=False))
    stats.to_csv(Path(args.out_dir) / "supp_set_size_by_condition_tests.csv", index=False)


if __name__ == "__main__":
    main()
