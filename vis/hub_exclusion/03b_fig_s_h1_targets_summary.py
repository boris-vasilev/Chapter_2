#!/usr/bin/env python
"""fig_s_h1_targets_summary - target-set overlap vs threshold (supplement to
fig_h1). Target sets change only through anchor loss: a target set is the
out-edges of the trait's anchor perturbations, so a hub that is not one of
its anchors leaves it unchanged.

  A. Median (IQR) size-only excess Jaccard by pair class, per condition:
     all target sets (top) and target sets whose seeding anchors are unchanged
     at every threshold (bottom; identical sets across variants, so only the
     smaller target universe moves them; n pairs in h1b_excess_summary.csv).
     Related-GWAS pairs: faint dashed line, clipped where off scale.
  B. Spearman correlation of pairwise excess vs 0 %; solid all sets, dotted
     anchors unchanged.
Data: data/h1b_excess_summary.csv, h1c_spearman.csv.
"""

import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt

import hxcommon  # noqa: F401  (sets sys.path, CH2_VARIANT)
from common import CONDS, PRINT_W, PT, axes_in, panel_label, print_style, save
from hxcommon import SUPP, read
from hxplot import cond_handles, threshold_axis

_s = importlib.util.spec_from_file_location(
    "h1", Path(__file__).with_name("02_fig_h1_overlap.py"))
h1 = importlib.util.module_from_spec(_s)
_s.loader.exec_module(h1)

ROWS = [("Targets, all sets", "All target sets"),
        ("Targets, anchors unchanged", "Anchors\nunchanged")]


def main():
    summ, rho = read("h1b_excess_summary.csv"), read("h1c_spearman.csv")
    main_ = summ[summ.pair_class != "related"]
    rel = summ[summ.pair_class == "related"]
    with print_style():
        W, H = PRINT_W, 2.75
        fig = plt.figure(figsize=(W, H))
        panel_label(fig, 0.02, H - 0.05, "A", "Target sets: excess Jaccard (size-only)")
        rows = [(main_, st, lab) for st, lab in ROWS]
        _, ax0 = h1.class_panels(fig, 0.7, H - 0.62, rows,
                                 [rel[rel.set_type == st] for st, _ in ROWS], bh=0.7)
        q = main_[main_.set_type.isin([s for s, _ in ROWS])]
        lo, hi = q.q25.min(), q.q75.max()
        ax0.set_ylim(lo - 0.08 * (hi - lo), hi + 0.08 * (hi - lo))
        ax0.text(0.03, 0.97, "related pairs clipped", transform=ax0.transAxes, va="top",
                 fontsize=PT["tick"] - 1, color="#6f6e69")
        fig.legend(handles=cond_handles() + [h1.plt_line()], loc="upper left", ncol=4,
                   bbox_to_anchor=(0.7 / W, (H - 0.2) / H), fontsize=PT["tick"])
        panel_label(fig, 4.25, H - 0.05, "B", "Structure retained")
        ax = axes_in(fig, 4.65, H - 0.62 - 1.55, 1.45, 1.55)
        for c in CONDS:
            for st, ls in [("Targets, all sets", "-"), ("Targets, anchors unchanged", ":")]:
                g = rho[(rho.set_type == st) & (rho.condition == c)]
                h1.line(ax, g.variant, g.spearman_rho, c, ls=ls)
        ax.set_ylabel("Spearman ρ of\nexcess vs 0 %", fontsize=PT["tick"])
        ax.text(0.03, 0.04, "dotted: anchors unchanged", transform=ax.transAxes,
                fontsize=PT["tick"] - 1, color="#52514e")
        threshold_axis(ax)
        save(fig, "fig_s_h1_targets_summary", SUPP, exact=True)


if __name__ == "__main__":
    main()
