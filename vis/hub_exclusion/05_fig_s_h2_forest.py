#!/usr/bin/env python
"""fig_s_h2_forest_<cond> - full S-LDSC forest per condition: anchors,
regulators and targets x (log2 enrichment, tau*), 95 % CI, one dodged point
per variant. Filled = FDR < 0.05 (BH within each variant: anchors and
regulators + targets as separate families, as Fig. 4; tau* one-sided).
Diamond = target set that lost >= 1 seeding anchor as a hub. Grey square =
no estimate. Anchor shifts across variants come only from the joint fit
(data/anchor_shifts.csv). Data: data/heritability_long.csv.
"""

import matplotlib.pyplot as plt

import hxcommon  # noqa: F401  (sets sys.path, CH2_VARIANT)
from common import COND_LABEL, CONDS, PRINT_W, PT, axes_in, panel_label, print_style, save
from hcommon import H_ANNOT_LABEL, H_ANNOTS
from hxcommon import SUPP, load_diseases, read
from hxplot import forest, variant_handles

import importlib.util
from pathlib import Path
_s = importlib.util.spec_from_file_location(
    "h2", Path(__file__).with_name("04_fig_h2_heritability.py"))
h2 = importlib.util.module_from_spec(_s)
_s.loader.exec_module(h2)


def main():
    diseases = load_diseases()
    h = read("heritability_long.csv")
    order = h2.trait_order(diseases)
    for cond in CONDS:
        with print_style():
            rowh, colw, gap, x0 = 0.165, 0.72, 0.05, 1.62
            side = len(order) * rowh
            W, H = PRINT_W, side + 1.3
            fig = plt.figure(figsize=(W, H))
            panel_label(fig, 0.02, H - 0.05, "", f"S-LDSC estimates, {COND_LABEL[cond]}")
            k = 0
            for annot in H_ANNOTS:
                sub = h[(h.annot == annot) & (h.condition == cond)]
                for val, se, sig, lab in [
                        ("log2_enrichment", "log2_enrichment_se", "enrichment_sig", "log$_2$ enr."),
                        ("tau_star", "tau_star_se", "tau_star_sig", "τ*")]:
                    ax = axes_in(fig, x0 + k * (colw + gap), H - 0.55 - side, colw, side)
                    forest(ax, sub, order, val, se, sig, lab, diseases, ylabels=(k == 0))
                    h2.xlim_from(ax, sub[val])
                    if val == "log2_enrichment":
                        ax.set_title(H_ANNOT_LABEL[annot], loc="left", fontsize=PT["tick"])
                    k += 1
            fig.legend(handles=variant_handles(), ncol=4, loc="lower left",
                       bbox_to_anchor=(0.02, 0.0), title="% of regulators excluded",
                       title_fontsize=PT["tick"], alignment="left")
            save(fig, f"fig_s_h2_forest_{cond}", SUPP, exact=True)


if __name__ == "__main__":
    main()
