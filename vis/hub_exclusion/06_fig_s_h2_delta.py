#!/usr/bin/env python
"""fig_s_h2_delta_<cond> - change vs 0 % per condition for regulator and
target sets: delta log2 enrichment and delta tau* with 95 % CI.
SE = sqrt(se_variant^2 + se_deg0^2), conservative: ignores positive
covariance between nested annotations. Filled = BH FDR < 0.05 of the
two-sided delta test within each variant over all regulator + target sets.
Diamond = target set that lost >= 1 seeding anchor. Regulator deltas also
reflect the target annotation changing within the same joint fit.
Data: data/h2c_delta.csv.
"""

import matplotlib.pyplot as plt

import hxcommon  # noqa: F401  (sets sys.path, CH2_VARIANT)
from common import COND_LABEL, CONDS, PRINT_W, PT, axes_in, panel_label, print_style, save
from hcommon import H_ANNOT_LABEL
from hxcommon import CONSERVATIVE, SET_ANNOTS, SUPP, VKEYS, load_diseases, read
from hxplot import forest, variant_handles

import importlib.util
from pathlib import Path
_s = importlib.util.spec_from_file_location(
    "h2", Path(__file__).with_name("04_fig_h2_heritability.py"))
h2 = importlib.util.module_from_spec(_s)
_s.loader.exec_module(h2)


def main():
    diseases = load_diseases()
    d = read("h2c_delta.csv")
    order = h2.trait_order(diseases)
    for cond in CONDS:
        with print_style():
            rowh, colw, gap, x0 = 0.165, 1.07, 0.06, 1.62
            side = len(order) * rowh
            W, H = PRINT_W, side + 1.45
            fig = plt.figure(figsize=(W, H))
            panel_label(fig, 0.02, H - 0.05, "", f"Change vs 0 % excluded, {COND_LABEL[cond]}")
            k = 0
            for annot in SET_ANNOTS:
                sub = d[(d.annot == annot) & (d.condition == cond)]
                for val, lab in [("delta_log2_enrichment", "Δ log$_2$ enrichment"),
                                 ("delta_tau_star", "Δτ*")]:
                    ax = axes_in(fig, x0 + k * (colw + gap), H - 0.55 - side, colw, side)
                    forest(ax, sub, order, val, val + "_se", val + "_sig", lab, diseases,
                           ylabels=(k == 0))
                    h2.xlim_from(ax, sub[val])
                    if val == "delta_log2_enrichment":
                        ax.set_title(H_ANNOT_LABEL[annot], loc="left", fontsize=PT["tick"])
                    k += 1
            fig.text(0.02, 0.02, f"95 % CI: {CONSERVATIVE}.", fontsize=PT["tick"] - 0.5,
                     color="#52514e")
            fig.legend(handles=variant_handles()[1:], ncol=3, loc="lower left",
                       bbox_to_anchor=(0.02, 0.2 / H), title="% of regulators excluded",
                       title_fontsize=PT["tick"], alignment="left")
            save(fig, f"fig_s_h2_delta_{cond}", SUPP, exact=True)


if __name__ == "__main__":
    main()
