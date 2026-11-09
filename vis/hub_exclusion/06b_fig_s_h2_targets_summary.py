#!/usr/bin/env python
"""fig_s_h2_targets_summary - target-set heritability vs threshold
(supplement to fig_h2). Target sets change only through anchor loss: a
target set is the out-edges of the trait's anchor perturbations.

  A. All target sets, B. target sets whose seeding anchors are unchanged at
     every threshold (7 IMD / 11 control trait x condition sets; few fitted):
     median log2 enrichment, median tau*, FDR-significant counts, IMD vs
     control, pooled over conditions.
  C. Precision of target sets: SE of enrichment and 95 % CI width of log2
     enrichment vs % of SNPs in the annotation, by variant.
Data: data/h2b_summary.csv, h2d_precision.csv.
"""

import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt

import hxcommon  # noqa: F401  (sets sys.path, CH2_VARIANT)
from common import (GROUP_COLORS, GROUP_LINESTYLE, GROUP_SHORT, PRINT_W, panel_label,
                    print_style, save)
from hxcommon import SUPP, read
from hxplot import variant_handles

_s = importlib.util.spec_from_file_location(
    "h2", Path(__file__).with_name("04_fig_h2_heritability.py"))
h2 = importlib.util.module_from_spec(_s)
_s.loader.exec_module(h2)


def main():
    summ, prec = read("h2b_summary.csv"), read("h2d_precision.csv")
    with print_style():
        W, H = PRINT_W, 5.95
        fig = plt.figure(figsize=(W, H))
        top = H - 0.75
        panel_label(fig, 0.02, H - 0.05, "A", "Targets, all sets")
        h2.b_panels(fig, 0.45, top, summ, "Targets, all sets", bw=1.0, bh=0.75)
        panel_label(fig, 3.25, H - 0.05, "B", "Targets, anchors unchanged")
        h2.b_panels(fig, 3.6, top, summ, "Targets, anchors unchanged", bw=1.0, bh=0.75)
        fig.legend(handles=[plt.Line2D([], [], color=GROUP_COLORS[g], ls=GROUP_LINESTYLE[g],
                                       label=GROUP_SHORT[g]) for g in ["IMD", "non-IMD"]],
                   loc="upper left", ncol=2, bbox_to_anchor=(0.45 / W, (H - 0.22) / H))
        ctop = top - 2 * 0.75 - 0.38 - 0.7
        panel_label(fig, 0.02, ctop + 0.25, "C", "Precision, target sets")
        h2.d_panels(fig, 0.75, ctop, prec, "targets_H4", w=1.5, hgt=0.9)
        fig.legend(handles=variant_handles(), loc="upper left", ncol=1,
                   bbox_to_anchor=(2.6 / W, ctop / H), title="% of regulators excluded")
        save(fig, "fig_s_h2_targets_summary", SUPP, exact=True)


if __name__ == "__main__":
    main()
