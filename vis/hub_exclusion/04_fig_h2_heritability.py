#!/usr/bin/env python
"""fig_h2_heritability - do S-LDSC enrichment and tau* change with hub exclusion?

All estimates come from joint per-condition fits (anchors + regulators +
targets + baselineLD), so each tau* is conditional on the other two sets in
its run; a regulator's change also reflects its target set changing.

  A. Regulator sets, Stim 8 h: log2 enrichment and tau* with 95 % CI per
     trait, one dodged point per variant. Filled = FDR < 0.05 (BH within each
     variant over all regulator + target sets, as Fig. 4; tau* one-sided).
  C. Change vs 0 % for the same sets: delta log2 enrichment and delta tau*
     with 95 % CI. SE = sqrt(se_variant^2 + se_deg0^2) - conservative: ignores
     positive covariance between nested annotations (no LDSC delete values to
     pair the jackknife blocks). Filled = BH FDR < 0.05 of the two-sided delta
     test within each variant over all regulator + target sets.
  B. Regulator sets: medians and FDR-significant counts, IMD vs control, vs
     threshold, pooled over conditions.
  D. Precision of regulator sets: SE of enrichment and 95 % CI width of log2
     enrichment vs % of SNPs in the annotation, by variant; grey lines join the
     same set across variants. (SE(tau*) is not used: tau* is scaled by
     sd(annotation), so it does not compare power across annotation sizes.)

Diamonds = target sets that lost >= 1 seeding anchor as a hub (open/filled
still encodes significance). Grey squares = no estimate (set empty or below
the 1,000-SNP threshold). Other conditions and targets (which change only
through anchor loss): fig_s_h2_*.

Writes data/h2a_estimates.csv, h2b_summary.csv, h2c_delta.csv, h2d_precision.csv.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

import hxcommon  # noqa: F401  (sets sys.path, CH2_VARIANT)
from common import (plain_log_axis, GROUP_COLORS, GROUP_LINESTYLE, GROUP_SHORT, PRINT_W, PT,
                    axes_in, panel_label, print_style, save)
from hcommon import FDR
from hxcommon import (CONSERVATIVE, DOUT, MAIN, MAIN_COND, PCT, SET_ANNOTS, VCOLOR,
                      VKEYS, load_diseases, read)
from hxplot import XPOS, forest, threshold_axis, variant_handles

SET_TYPES = [("regulators_H4", False, "Regulators"),
             ("targets_H4", False, "Targets, all sets"),
             ("targets_H4", True, "Targets, anchors unchanged")]


def load():
    h = read("heritability_long.csv")
    return h[h.annot.isin(SET_ANNOTS)].copy()


def delta_table(h):
    key = ["study", "condition", "annot"]
    base = h[h.variant == "deg0"].set_index(key)
    rows = []
    for v in VKEYS[1:]:
        x = h[h.variant == v].set_index(key).reindex(base.index)
        d = pd.DataFrame(index=base.index)
        d["variant"], d["pct_excluded"] = v, PCT[v]
        for c in ["disease", "group", "lost_anchor", "n_lost_anchors", "unchanged_10"]:
            d[c] = x[c]
        d["status"], d["status_deg0"] = x["status"], base["status"]
        for val, se, name in [("log2_enrichment", "log2_enrichment_se", "log2_enrichment"),
                              ("tau_star", "tau_star_se", "tau_star")]:
            d[f"{name}_deg0"], d[name] = base[val], x[val]
            d[f"delta_{name}"] = x[val] - base[val]
            d[f"delta_{name}_se"] = np.sqrt(x[se] ** 2 + base[se] ** 2)
            d[f"delta_{name}_z"] = d[f"delta_{name}"] / d[f"delta_{name}_se"]
            d[f"delta_{name}_p"] = 2 * norm.sf(d[f"delta_{name}_z"].abs())
            ok = d[f"delta_{name}_p"].notna()
            d[f"delta_{name}_fdr"] = np.nan
            if ok.any():
                d.loc[ok, f"delta_{name}_fdr"] = multipletests(
                    d.loc[ok, f"delta_{name}_p"], method="fdr_bh")[1]
            d[f"delta_{name}_sig"] = d[f"delta_{name}_fdr"] < FDR
        d["se_method"] = CONSERVATIVE
        rows.append(d.reset_index())
    return pd.concat(rows, ignore_index=True)


def summary_table(h):
    rows = []
    for annot, unchanged, label in SET_TYPES:
        s = h[h.annot == annot]
        if unchanged:
            s = s[s.unchanged_10.astype(bool)]
        for cond in ["all"] + sorted(s.condition.unique()):
            sc = s if cond == "all" else s[s.condition == cond]
            for (v, g), x in sc.groupby(["variant", "group"]):
                f = x[x.status == "fitted"]
                rows.append({"set_type": label, "condition": cond, "variant": v,
                             "pct_excluded": PCT[v], "group": g, "n_sets": len(x),
                             "n_fitted": len(f),
                             "median_log2_enrichment": f.log2_enrichment.median(),
                             "median_tau_star": f.tau_star.median(),
                             "n_sig_enrichment": int(f.enrichment_sig.astype(bool).sum()),
                             "n_sig_tau_star": int(f.tau_star_sig.astype(bool).sum())})
    return pd.DataFrame(rows)


def trait_order(diseases):
    """Fig. 4 order: IMDs then controls, each by anchor count (descending)."""
    from common import disease_order
    h = read("heritability_long.csv")
    a = h[(h.variant == "deg0") & (h.annot == "H4_genes")].groupby("study").n_genes.max()
    return [s for s in disease_order(diseases, a) if s in a.index]


def xlim_from(ax, vals, pad=0.12):
    v = pd.Series(vals).replace([np.inf, -np.inf], np.nan).dropna()
    lo, hi = min(v.min(), 0), max(v.max(), 0)
    ax.set_xlim(lo - pad * (hi - lo), hi + pad * (hi - lo))


def forest_row(fig, x0, y0, colw, gap, rowh, h, d, order, diseases):
    """Four columns: log2 enrichment, tau*, delta log2 enrichment, delta tau*."""
    side = len(order) * rowh
    cols = [(h, "log2_enrichment", "log2_enrichment_se", "enrichment_sig", "log$_2$ enrichment"),
            (h, "tau_star", "tau_star_se", "tau_star_sig", "τ*"),
            (d, "delta_log2_enrichment", "delta_log2_enrichment_se",
             "delta_log2_enrichment_sig", "Δ log$_2$ enrichment"),
            (d, "delta_tau_star", "delta_tau_star_se", "delta_tau_star_sig", "Δτ*")]
    axes = []
    for k, (src, val, se, sig, lab) in enumerate(cols):
        ax = axes_in(fig, x0 + k * (colw + gap), y0, colw, side)
        forest(ax, src, order, val, se, sig, lab, diseases, ylabels=(k == 0))
        xlim_from(ax, src[val])
        axes.append(ax)
    return axes


B_METRICS = [("median_log2_enrichment", "Median log$_2$ enrichment"),
             ("median_tau_star", "Median τ*"),
             ("n_sig_enrichment", "Sets FDR < 0.05 (enrichment)"),
             ("n_sig_tau_star", "Sets FDR < 0.05 (τ*)")]


def b_panels(fig, x0, top, summ, set_type, bw=1.25, bh=0.8):
    """2 x 2 grid of IMD vs control summaries vs threshold for one set type,
    pooled over conditions."""
    s = summ[(summ.condition == "all") & (summ.set_type == set_type)]
    for i, (m, title) in enumerate(B_METRICS):
        r, c = divmod(i, 2)
        ax = axes_in(fig, x0 + c * (bw + 0.45), top - (r + 1) * bh - r * 0.38, bw, bh)
        for g in ["IMD", "non-IMD"]:
            x = s[s.group == g].sort_values("pct_excluded")
            ax.plot([XPOS[v] for v in x.variant], x[m], color=GROUP_COLORS[g],
                    ls=GROUP_LINESTYLE[g], marker="o", ms=2.2, lw=0.8)
        ax.axhline(0, color="#9a9994", lw=0.5, zorder=0)
        threshold_axis(ax, label=(r == 1))
        if r == 0:
            ax.set_xticklabels([])
        ax.set_title(title, fontsize=PT["tick"], pad=3)
    return top - 2 * bh - 0.38


D_METRICS = [("enrichment_se", "SE(enrichment)"),
             ("ci_width_log2_enrichment", "95 % CI width,\nlog$_2$ enrichment")]


def precision_table(h):
    p = h[["variant", "pct_excluded", "study", "disease", "group", "condition", "annot",
           "status", "n_genes", "prop_snps", "enrichment", "enrichment_se",
           "log2_enrichment", "log2_enrichment_se", "tau_star", "tau_star_se",
           "lost_anchor"]].copy()
    p["ci_width_enrichment"] = 2 * 1.96 * p.enrichment_se
    p["ci_width_log2_enrichment"] = 2 * 1.96 * p.log2_enrichment_se
    return p


def d_panels(fig, x0, top, p, annot, w=1.5, hgt=1.05):
    """Enrichment SE and log2-enrichment CI width vs % of SNPs in the
    annotation, by variant; grey lines join the same set across variants."""
    f = p[(p.annot == annot) & (p.status == "fitted")]
    for i, (m, ylab) in enumerate(D_METRICS):
        ax = axes_in(fig, x0, top - (i + 1) * hgt - i * 0.42, w, hgt)
        for _, g in f.groupby(["study", "condition"]):
            g = g.set_index("variant").reindex(VKEYS).dropna(subset=["prop_snps", m])
            ax.plot(g.prop_snps * 100, g[m], color="#c9c8c3", lw=0.4, zorder=1)
        for v in VKEYS:
            x = f[f.variant == v]
            ax.scatter(x.prop_snps * 100, x[m], s=5, color=VCOLOR[v], lw=0,
                       zorder=2 + VKEYS.index(v))
        ax.set_xscale("log")
        ax.set_yscale("log")
        plain_log_axis(ax, "x")
        plain_log_axis(ax, "y")
        ax.set_ylabel(ylab, fontsize=PT["tick"])
        ax.set_xlabel("% of SNPs in annotation" if i == len(D_METRICS) - 1 else "")


def main():
    diseases = load_diseases()
    h = load()
    d = delta_table(h)
    d.to_csv(DOUT / "h2c_delta.csv", index=False)
    summ = summary_table(h)
    summ.to_csv(DOUT / "h2b_summary.csv", index=False)
    h.to_csv(DOUT / "h2a_estimates.csv", index=False)
    prec = precision_table(h)
    prec.to_csv(DOUT / "h2d_precision.csv", index=False)

    order = trait_order(diseases)
    sel = lambda df: df[(df.annot == "regulators_H4") & (df.condition == MAIN_COND)]  # noqa: E731
    with print_style():
        W, H = PRINT_W, 7.7
        fig = plt.figure(figsize=(W, H))
        rowh, top = 0.165, H - 0.45
        side = len(order) * rowh
        panel_label(fig, 0.02, H - 0.05, "A", "Regulator sets, Stim 8 h")
        panel_label(fig, 1.62 + 2 * 1.13, H - 0.05, "C", "Change vs 0 %")
        forest_row(fig, 1.62, top - side, 1.07, 0.06, rowh, sel(h), sel(d), order, diseases)
        fig.text(0.02, (top - side - 0.42) / H,
                 f"Δ 95 % CI: {CONSERVATIVE}.", fontsize=PT["tick"] - 0.5, color="#52514e")
        sig = [plt.Line2D([], [], ls="", marker="o", mfc="#52514e", mec="#52514e", ms=3,
                          label="FDR < 0.05"),
               plt.Line2D([], [], ls="", marker="o", mfc="white", mec="#52514e", ms=3,
                          label="n.s."),
               plt.Line2D([], [], ls="", marker="s", color="#dcdbd6", ms=3, label="No estimate")]
        fig.legend(handles=variant_handles() + sig, ncol=7, loc="upper left",
                   bbox_to_anchor=(0.02, (top - side - 0.48) / H), fontsize=PT["tick"],
                   columnspacing=0.8, handletextpad=0.2, title="% of regulators excluded",
                   title_fontsize=PT["tick"], alignment="left")

        btop = top - side - 0.95
        panel_label(fig, 0.02, btop, "B", "Regulator sets, IMD vs control, pooled over conditions")
        b_panels(fig, 0.45, btop - 0.6, summ, "Regulators")
        fig.legend(handles=[plt.Line2D([], [], color=GROUP_COLORS[g], ls=GROUP_LINESTYLE[g],
                                       label=GROUP_SHORT[g]) for g in ["IMD", "non-IMD"]],
                   loc="upper left", ncol=2, bbox_to_anchor=(0.45 / W, (btop - 0.12) / H))
        panel_label(fig, 3.95, btop, "D", "Precision, regulator sets")
        d_panels(fig, 4.6, btop - 0.35, prec, "regulators_H4")
        save(fig, "fig_h2_heritability", MAIN, exact=True)


if __name__ == "__main__":
    main()
