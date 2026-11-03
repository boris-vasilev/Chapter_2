#!/usr/bin/env python
"""Collect cell-type-specific S-LDSC programme results into one tidy CSV.

Each LDSC run for this module fits ``baseline + one programme``, so the
programme is the single non-baseline annotation in the ``.results`` file.
LDSC names that thin annotation ``ANNOT`` (suffixed ``L2_0``); we pull that
row's coefficient and enrichment.

The programme coefficient z-score is the cell-type-specific enrichment
statistic (Finucane 2018): a positive z means the programme's genes
contribute more per-SNP heritability than expected given the baseline. We
report the one-sided p = P(Z > z) and a BH-FDR across every (study,
condition, programme) test.

Output columns:
  study, trait_short, condition, cluster, program, manual_annotation, n_genes,
  prop_h2, prop_h2_se, enrichment, enrichment_se, enrichment_p,
  coefficient, coefficient_se, coefficient_z, coef_p_onesided, coef_fdr_bh,
  ancestry, build
"""

import argparse
import re
from pathlib import Path

import polars as pl
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

FNAME_RE = re.compile(r"^(?P<cond>.+)_program(?P<cluster>\d+)\.results$")


def parse_one(results_file):
    """Return the programme (non-baseline) annotation row, or None."""
    df = pl.read_csv(results_file, separator="\t")
    # The programme annotation is the thin one named 'ANNOT' (+ L2 suffix).
    matches = [c for c in df["Category"].to_list() if c.startswith("ANNOT")]
    if not matches:
        return None
    r = df.filter(pl.col("Category") == matches[0]).row(0, named=True)
    return {
        "prop_h2": r["Prop._h2"],
        "prop_h2_se": r["Prop._h2_std_error"],
        "enrichment": r["Enrichment"],
        "enrichment_se": r["Enrichment_std_error"],
        "enrichment_p": r["Enrichment_p"],
        "coefficient": r["Coefficient"],
        "coefficient_se": r["Coefficient_std_error"],
        "coefficient_z": r["Coefficient_z-score"],
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--results-dir", required=True)
    p.add_argument("--study-info", required=True)
    p.add_argument("--programs-index", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    rows = []
    for study_dir in sorted(Path(args.results_dir).iterdir()):
        if not study_dir.is_dir():
            continue
        for f in sorted(study_dir.glob("*_program*.results")):
            m = FNAME_RE.match(f.name)
            if not m:
                continue
            parsed = parse_one(f)
            if parsed is None:
                print(f"  WARNING: no programme row in {f}")
                continue
            rows.append({
                "study": study_dir.name,
                "condition": m.group("cond"),
                "cluster": int(m.group("cluster")),
                **parsed,
            })

    if not rows:
        raise SystemExit("No programme results parsed.")

    df = pl.DataFrame(rows)

    # One-sided p from the coefficient z-score (enrichment = z > 0).
    z = df["coefficient_z"].to_numpy()
    df = df.with_columns(
        pl.Series("coef_p_onesided", norm.sf(z))
    )
    fdr = multipletests(df["coef_p_onesided"].to_numpy(), method="fdr_bh")[1]
    df = df.with_columns(pl.Series("coef_fdr_bh", fdr))

    # Labels.
    index = pl.read_csv(args.programs_index)
    info_raw = pl.read_csv(args.study_info)
    info = info_raw.select(
        [c for c in ["study", "ancestry", "build", "trait"] if c in info_raw.columns]
    )

    out = (
        df.join(index, on=["condition", "cluster"], how="left")
        .join(info, on="study", how="left")
        .with_columns(pl.col("study").str.split("_").list.first().alias("trait_short"))
        .sort(["condition", "study", "coef_p_onesided"])
    )

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(args.out)
    n_sig = int((out["coef_fdr_bh"] < 0.05).sum())
    print(f"Wrote {args.out}: {out.height} rows, "
          f"{out['study'].n_unique()} studies, "
          f"{n_sig} programme tests at BH-FDR<0.05")


if __name__ == "__main__":
    main()
