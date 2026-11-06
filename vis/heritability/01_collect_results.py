#!/usr/bin/env python
"""Collect S-LDSC results into one tidy table (heritability_long.csv).

One row per trait x condition x annotation (anchors, regulators, targets),
including cells with no estimate (status column says why). Reads
``LDSC_DIR/ldsc_results/<study>/<cond>_combined.{results,log}`` and the
gene sets; never re-runs S-LDSC.

tau* (standardised effect size, Gazal et al. 2017):

    tau* = tau * sd(a) * M / h2g,     sd(a) = sqrt(p (1 - p)),  p = Prop._SNPs

The M_5_50 files are not stored locally, so M / h2g is taken from the
.results file itself: ldsc computes h2g = sum_k tau_k M_k and writes
Prop._SNPs_k = M_k / M (ldscore/regressions.py, Hsq._tot and _overlap_output),
hence h2g / M = sum_k tau_k Prop._SNPs_k over all categories (baseline +
custom). The implied M = h2g(log) / sum_k tau_k p_k is written as a check.
SE(tau*) = SE(tau) * sd(a) * M / h2g (h2g and M treated as fixed); the
tau* p-value is one-sided (H1: tau* > 0), p = Phi(-z) with z = tau / SE(tau);
the two-sided p and its FDR are kept as *_two_sided columns.

FDR: Benjamini-Hochberg within two families, expanded sets (regulators +
targets) and anchors, separately for enrichment and tau*. Sets with negative
Prop._h2 are flagged (neg_prop_h2) and never called significant
(*_sig columns).

Also writes ldsc_runs.csv (one row per study x condition: LD reference,
h2g, implied M, status) and the Step 0 file inventory.
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, t as tdist
from statsmodels.stats.multitest import multipletests

from hcommon import (EMPTY, FAILED, FDR, FITTED, H_ANNOTS, LONG_CSV, NOT_RUN,
                     add_args, dropped_status, load_diseases, resolve)
from common import CONDS, load_gene_sets

# Runs where only one custom annotation survived combine_annots.py
# (--min-snps 1000). ldsc then names the LD-score column plain "L2", so the
# .results category is "L2_1" with no set name. Identified from the gene sets
# (see NOTES_heritability.md, Step 0): the other non-empty sets in these runs
# span < 1,000 1000G SNPs (gene body +/- 10 kb), and for osteoarthritis
# Prop._SNPs tracks the regulator count across conditions (108 / 111 / 49
# genes -> 0.22 / 0.21 / 0.09 %).
SINGLE_COLUMN = {
    ("ADFH_Jansen_30617256_1", "Stim48hr"): "targets_H4",
    ("OA_McDonald_36411363_3", "Rest"): "regulators_H4",
    ("OA_McDonald_36411363_3", "Stim8hr"): "regulators_H4",
    ("OA_McDonald_36411363_3", "Stim48hr"): "regulators_H4",
}

LD_REF = {  # --frqfile-chr prefix -> (panel label, ancestry)
    "GRCh37/1000G_Phase3_plink/1000G.EUR.QC.": ("1000G EUR (GRCh37)", "EUR"),
    "GRCh38/1000G_Phase3_plink/1000G.EUR.hg38.": ("1000G EUR (GRCh38)", "EUR"),
    "GRCh37_EAS/1000G_Phase3_plink/1000G.EAS.QC.": ("1000G EAS (GRCh37)", "EAS"),
}
N_JK_BLOCKS = 200  # ldsc default; Enrichment_p uses a t test with this df


def parse_log(path):
    txt = path.read_text()
    frq = re.search(r"--frqfile-chr \S*/data/(\S+)", txt)
    h2 = re.search(r"Total Observed scale h2: (\S+) \((\S+)\)", txt)
    cats = re.search(r"^Categories: (.*)$", txt, re.M)
    return {
        "frq": frq.group(1) if frq else None,
        "h2_obs": float(h2.group(1)) if h2 else np.nan,
        "h2_obs_se": float(h2.group(2)) if h2 else np.nan,
        "custom": [c for c in cats.group(1).split() if c.endswith("L2_1")] if cats else [],
        "error": txt.splitlines()[-4].strip() if "Traceback" in txt else "",
    }


def kept_annots(combine_txt):
    """Annotations combine_annots.py kept ("[KEEP] <annot>: ...") in this run."""
    return re.findall(r"\[KEEP\] (\S+):", combine_txt.read_text())


def category_map(study, cond, res, run_dir):
    """annot -> .results row for the custom annotations in this run.

    An unnamed "L2_1" column (single surviving annotation) is identified from
    the pipeline's <cond>_combined.combine.txt when present (newer runs), else
    from SINGLE_COLUMN (older runs, which did not write that file)."""
    custom = res[res.Category.str.endswith("L2_1")]
    out = {}
    for _, r in custom.iterrows():
        if r.Category == "L2_1":
            combine = Path(run_dir) / f"{cond}_combined.combine.txt"
            kept = kept_annots(combine) if combine.exists() else []
            if len(kept) == 1:
                out[kept[0]] = r
            elif (study, cond) in SINGLE_COLUMN:
                out[SINGLE_COLUMN[(study, cond)]] = r
            else:
                raise SystemExit(f"{study} {cond}: unnamed L2_1 column; not in "
                                 f"SINGLE_COLUMN and {combine.name} lists {kept}")
        else:
            out[r.Category[: -len("L2_1")]] = r
    return out


def collect(ldsc_dir, sets, diseases):
    rows, runs = [], []
    for study in diseases.index:
        for cond in CONDS:
            d = Path(ldsc_dir) / study
            res_f, log_f = d / f"{cond}_combined.results", d / f"{cond}_combined.log"
            log = parse_log(log_f) if log_f.exists() else None
            panel, ref_anc = LD_REF.get(log["frq"], (log["frq"], None)) if log else (None, None)
            run = {"study": study, "condition": cond, "ld_reference": panel,
                   "ld_reference_ancestry": ref_anc,
                   "h2_obs": log["h2_obs"] if log else np.nan,
                   "h2_obs_se": log["h2_obs_se"] if log else np.nan,
                   "files": ", ".join(sorted(p.suffix for p in d.glob(f"{cond}_combined.*"))),
                   "error": log["error"] if log else ""}
            if res_f.exists():
                res = pd.read_csv(res_f, sep="\t")
                h2_over_m = float((res["Coefficient"] * res["Prop._SNPs"]).sum())
                run.update(status="results", h2_over_M=h2_over_m,
                           M_implied=run["h2_obs"] / h2_over_m,
                           n_custom=int(res.Category.str.endswith("L2_1").sum()))
                fitted = category_map(study, cond, res, d)
            else:
                run.update(status="run failed" if log else "not run",
                           h2_over_M=np.nan, M_implied=np.nan, n_custom=0)
                fitted = {}
            runs.append(run)

            for annot in H_ANNOTS:
                n = len(sets.get((study, cond, annot), ()))
                row = {"study": study, "condition": cond, "annot": annot, "n_genes": n,
                       "ld_reference": panel}
                if annot in fitted:
                    r = fitted[annot]
                    p = r["Prop._SNPs"]
                    sd = np.sqrt(p * (1 - p))
                    scale = sd / run["h2_over_M"]
                    row.update(
                        status=FITTED, category=r.Category,
                        prop_snps=p, prop_h2=r["Prop._h2"],
                        prop_h2_se=r["Prop._h2_std_error"],
                        enrichment=r["Enrichment"], enrichment_se=r["Enrichment_std_error"],
                        enrichment_p=r["Enrichment_p"],
                        tau=r["Coefficient"], tau_se=r["Coefficient_std_error"],
                        tau_z=r["Coefficient_z-score"], annot_sd=sd,
                        tau_star=r["Coefficient"] * scale,
                        tau_star_se=r["Coefficient_std_error"] * scale,
                        tau_star_p=norm.sf(r["Coefficient_z-score"]),
                        tau_star_p_two_sided=2 * norm.sf(abs(r["Coefficient_z-score"])))
                elif n == 0:
                    row["status"] = EMPTY
                elif run["status"] == "results":
                    row["status"] = dropped_status(d / f"{cond}_combined.combine.txt", annot)
                else:
                    row["status"] = FAILED if run["status"] == "run failed" else NOT_RUN
                rows.append(row)
    return pd.DataFrame(rows), pd.DataFrame(runs)


def add_fdr(df):
    """BH within two families - expanded sets (regulators + targets) and
    anchors - separately for enrichment and tau*. Estimates with negative
    Prop._h2 stay in the correction but are never called significant."""
    df["fdr_family"] = np.where(df["annot"] == "H4_genes", "anchors", "expanded")
    df["neg_prop_h2"] = df["prop_h2"] < 0
    for p_col, q_col in [("enrichment_p", "enrichment_fdr"), ("tau_star_p", "tau_star_fdr"),
                         ("tau_star_p_two_sided", "tau_star_fdr_two_sided")]:
        for _, idx in df[df[p_col].notna()].groupby("fdr_family").groups.items():
            df.loc[idx, q_col] = multipletests(df.loc[idx, p_col], method="fdr_bh")[1]
    for q_col, sig_col in [("enrichment_fdr", "enrichment_sig"),
                           ("tau_star_fdr", "tau_star_sig"),
                           ("tau_star_fdr_two_sided", "tau_star_sig_two_sided")]:
        df[sig_col] = (df[q_col] < FDR) & ~df["neg_prop_h2"]
    # Signed z for enrichment, recovered from ldsc's two-sided t test
    # (df = jackknife blocks) so it matches Enrichment_p exactly.
    df["enrichment_z"] = (np.sign(df["enrichment"] - 1)
                          * tdist.isf(df["enrichment_p"] / 2, N_JK_BLOCKS))
    df["tau_star_z"] = df["tau_z"]
    return df


def main():
    ap = argparse.ArgumentParser()
    add_args(ap)
    args = resolve(ap.parse_args())

    diseases = load_diseases()
    sets = load_gene_sets(args.gene_sets_dir)
    long, runs = collect(args.ldsc_dir, sets, diseases)
    long = add_fdr(long)

    anc = diseases["ancestry"]
    runs["gwas_ancestry"] = runs["study"].map(anc)
    runs["ld_mismatch"] = (runs["ld_reference_ancestry"].notna()
                           & (runs["gwas_ancestry"] != runs["ld_reference_ancestry"]))
    mismatch = runs.groupby("study")["ld_mismatch"].any()
    long.insert(1, "disease", long["study"].map(diseases["display_name"]))
    long.insert(2, "group", long["study"].map(diseases["group"].astype(str)))
    long.insert(3, "gwas_ancestry", long["study"].map(anc))
    long.insert(4, "ld_mismatch", long["study"].map(mismatch))

    cols = ["study", "disease", "group", "gwas_ancestry", "ld_reference", "ld_mismatch",
            "condition", "annot", "status", "n_genes", "category", "prop_snps",
            "prop_h2", "prop_h2_se", "neg_prop_h2", "fdr_family", "enrichment", "enrichment_se", "enrichment_p",
            "enrichment_z", "enrichment_fdr", "enrichment_sig", "tau", "tau_se", "tau_z",
            "annot_sd", "tau_star", "tau_star_se", "tau_star_z", "tau_star_p",
            "tau_star_fdr", "tau_star_sig", "tau_star_p_two_sided",
            "tau_star_fdr_two_sided", "tau_star_sig_two_sided"]
    long[cols].to_csv(args.out_dir / LONG_CSV, index=False)
    runs.to_csv(args.out_dir / "ldsc_runs.csv", index=False)
    print(f"Wrote {args.out_dir / LONG_CSV} ({len(long)} rows, "
          f"{(long.status == FITTED).sum()} fitted) and ldsc_runs.csv")


if __name__ == "__main__":
    main()
