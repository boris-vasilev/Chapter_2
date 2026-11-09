#!/usr/bin/env python
"""fig_s_h1_overlap_<set>_<cond> - supplement to fig_h1 A for every other
set type x condition: observed (top) and observed - expected Jaccard
(bottom) per variant, deg0 trait order, shared scale per row. Grey = empty
set; open circle on the diagonal = target set that lost >= 1 seeding anchor
as a hub. Data: data/jaccard_pairs.csv.
"""

import matplotlib.pyplot as plt

import hxcommon  # noqa: F401  (sets sys.path, CH2_VARIANT)
from common import ANNOT_LABEL, COND_LABEL, CONDS, PRINT_W, panel_label, print_style, save
from hxcommon import MAIN_COND, SET_ANNOTS, SUPP, load_diseases, read
from hxplot import heatmap_block


def main():
    diseases = load_diseases()
    pairs, sizes = read("jaccard_pairs.csv"), read("set_sizes.csv")
    for annot in SET_ANNOTS:
        for cond in CONDS:
            if (annot, cond) == ("regulators_H4", MAIN_COND):
                continue      # main figure
            with print_style():
                W, H = PRINT_W, 3.85
                fig = plt.figure(figsize=(W, H))
                panel_label(fig, 0.02, H - 0.05, "",
                            f"{ANNOT_LABEL[annot]}, {COND_LABEL[cond]}")
                heatmap_block(fig, 0.42, H - 0.45 - 1.30, 1.30, 0.07, pairs, sizes,
                              cond, annot, diseases)
                short = annot.split("_")[0]
                save(fig, f"fig_s_h1_overlap_{short}_{cond}", SUPP, exact=True)


if __name__ == "__main__":
    main()
