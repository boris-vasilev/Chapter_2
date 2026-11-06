#!/usr/bin/env python
"""Numbers for the chapter text: writes the "Heritability" section of
numbers_for_text.md (00_build_all.py copies it to figures/data/).

All counts are over fitted sets. Significant = the *_sig columns of
heritability_long.csv: BH FDR < 0.05 within the expanded-set family
(regulators + targets) and within anchors separately, for enrichment and for
one-sided tau* (H1: tau* > 0); sets with negative Prop._h2 never count.
"""

import argparse

import numpy as np

from hcommon import (FDR, FITTED, H_ANNOT_LABEL, H_ANNOTS, add_args, load_long,
                     resolve)
from common import COND_SHORT, CONDS, GROUPS, GROUP_SHORT, VARIANT

GROUP_COL = {"IMD": "IMD", "non-IMD": "control"}


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return out


def section(long):
    fit = long[long.status == FITTED].copy()
    fit["sig_e"] = fit.enrichment_sig.astype(bool)
    fit["sig_t"] = fit.tau_star_sig.astype(bool)
    neg = fit[fit.neg_prop_h2.astype(bool)]
    L = ["## Heritability", "",
         f"S-LDSC, design {VARIANT} (baselineLD v2.2 + anchors, regulators and "
         f"targets fitted jointly per trait x condition). {len(fit)} fitted sets of "
         f"{len(long)} trait x condition x annotation cells. Significant = BH FDR < "
         f"{FDR}, corrected within expanded sets (regulators + targets, "
         f"{(fit.annot != 'H4_genes').sum()} tests) and within anchors "
         f"({(fit.annot == 'H4_genes').sum()} tests) separately, for enrichment "
         f"(two-sided) and for τ* (one-sided, H1: τ* > 0).", "",
         "Flagged, never counted as significant (negative Prop._h2): " +
         ("; ".join(f"{r.disease} {COND_SHORT[r.condition]} "
                    f"{H_ANNOT_LABEL[r.annot].lower()} (Prop._h2 = {r.prop_h2:.4f}, "
                    f"enrichment = {r.enrichment:.1f})" for r in neg.itertuples())
          or "none") + ".", ""]

    L += ["### Sets significant, by annotation and condition (significant / fitted)", ""]
    rows = []
    for a in H_ANNOTS:
        for c in CONDS:
            d = fit[(fit.annot == a) & (fit.condition == c)]
            cells = []
            for g in GROUPS:
                dg = d[d.group == g]
                cells += [f"{dg.sig_e.sum()}/{len(dg)}", f"{dg.sig_t.sum()}/{len(dg)}"]
            rows.append([H_ANNOT_LABEL[a], COND_SHORT[c], *cells])
    L += md_table(["Annotation", "Condition", "Enrichment, IMD", "τ*, IMD",
                   "Enrichment, control", "τ*, control"], rows) + [""]

    L += ["### Median log2 enrichment (all conditions pooled)", ""]
    rows = []
    for a in H_ANNOTS:
        cells = []
        for g in GROUPS:
            d = fit[(fit.annot == a) & (fit.group == g)]
            e = d.enrichment
            cells.append(f"{np.median(np.log2(e[e > 0])):.2f} (n = {(e > 0).sum()}; "
                         f"median enrichment {e.median():.2f})")
        rows.append([H_ANNOT_LABEL[a], *cells])
    L += md_table(["Annotation", "IMD", "Control"], rows)
    L += ["", "log2 medians exclude sets with enrichment ≤ 0 (" +
          (", ".join(f"{r.disease} {COND_SHORT[r.condition]} {H_ANNOT_LABEL[r.annot].lower()}"
                     f" = {r.enrichment:.1f}" for r in fit[fit.enrichment <= 0].itertuples())
           or "none") + "); the plain median includes them.", ""]

    L += ["### Marginal vs conditional", ""]
    rows = []
    for a in H_ANNOTS:
        d = fit[fit.annot == a]
        rows.append([H_ANNOT_LABEL[a], len(d), (d.sig_e & d.sig_t).sum(),
                     (d.sig_e & ~d.sig_t).sum(), (~d.sig_e & d.sig_t).sum(),
                     (~d.sig_e & ~d.sig_t).sum()])
    L += md_table(["Annotation", "Fitted", "Both", "Enrichment only", "τ* only",
                   "Neither"], rows) + [""]
    exp = fit[fit.annot != "H4_genes"]
    for g in GROUPS:
        d = exp[exp.group == g]
        L.append(f"- Expanded sets (regulators + targets), {GROUP_SHORT[g][0].lower() + GROUP_SHORT[g][1:]}: "
                 f"{d.sig_e.sum()}/{len(d)} enriched, {d.sig_t.sum()}/{len(d)} τ*-significant, "
                 f"{(d.sig_e & ~d.sig_t).sum()} enriched but not τ*.")
    L.append("")

    L += ["### Anchors per trait (reference point)", "",
          "Anchor estimates differ slightly between conditions only because the "
          "joint model's other annotations change; Rest shown, range over "
          "conditions in brackets.", ""]
    rows = []
    anc = fit[fit.annot == "H4_genes"]
    for s, d in anc.groupby("study", sort=False):
        r = d[d.condition == "Rest"].iloc[0] if (d.condition == "Rest").any() else d.iloc[0]
        rows.append([r.disease, GROUP_COL[r.group], r.n_genes,
                     f"{r.enrichment:.1f} ± {r.enrichment_se:.1f} "
                     f"[{d.enrichment.min():.1f}–{d.enrichment.max():.1f}]",
                     f"{r.enrichment_p:.1e}", "yes" if r.sig_e else "no",
                     f"{r.tau_star:.2f} ± {r.tau_star_se:.2f} "
                     f"[{d.tau_star.min():.2f}–{d.tau_star.max():.2f}]",
                     f"{r.tau_star_p:.1e}", "yes" if r.sig_t else "no"])
    L += md_table(["Trait", "Class", "Anchors", "Enrichment ± SE", "p", "FDR < 0.05",
                   "τ* ± SE", "p", "FDR < 0.05"], rows) + [""]

    L += ["### Missing runs and empty sets", ""]
    miss = long[long.status != FITTED]
    for status, d in miss.groupby("status"):
        items = [f"{r.disease} {COND_SHORT[r.condition]} {H_ANNOT_LABEL[r.annot].lower()}"
                 f" (n = {r.n_genes})" for r in d.itertuples()]
        L.append(f"- **{status}** ({len(d)}): " + "; ".join(items))
    L.append("")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser()
    add_args(ap)
    args = resolve(ap.parse_args())
    long = load_long(args.long)
    out = args.out_dir / "numbers_for_text.md"
    out.write_text("# Numbers for the text — Chapter 2\n\n" + section(long))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
