#!/usr/bin/env python
"""Wrap LDSC's munge_sumstats.py with auto-detected column mappings.

Handles three sumstats formats seen in our OpenTargets-derived input:

  1. FinnGen R12  — `#chrom`, `pos`, `ref`, `alt`, `rsids`, `pval`,
                    `beta`, `sebeta`, ...
  2. GWAS Catalog "old harmonised"  — has `hm_rsid` and a full set of
                    `hm_*` columns (the harmonised values). Use those
                    and ignore the un-prefixed originals.
  3. GWAS Catalog "new harmonised" / unharmonised  — uses the bare
                    columns (`chromosome`, `base_pair_location`,
                    `effect_allele`, `other_allele`, `beta` or
                    `odds_ratio`, `standard_error`, `p_value`) plus a
                    marker column `hm_coordinate_conversion`/`hm_code`.

Within (3), the chosen `--snp` is `rsid` if present, else `variant_id`.
The chosen `--signed-sumstats` is `beta,0` if `beta` is present, else
`odds_ratio,1`. Other duplicates are added to `--ignore` so LDSC
doesn't pick them up.
"""

import argparse
import gzip
import subprocess

import polars as pl


def _header(file):
    with gzip.open(file, "rt") as f:
        return f.readline()


def _columns(file):
    return {c.strip() for c in _header(file).rstrip("\n").split("\t")}


def detect_format_args(file):
    """Return the format-specific ``munge_sumstats.py`` flags for ``file``."""
    cols = _columns(file)

    # ----- FinnGen -----------------------------------------------------
    if {"#chrom", "rsids", "sebeta"} <= cols:
        return [
            "--snp", "rsids",
            "--a1",  "alt",
            "--a2",  "ref",
            "--p",   "pval",
            "--signed-sumstats", "beta,0",
        ]

    # ----- GWAS Catalog "old harmonised" (hm_* family) ----------------
    if "hm_rsid" in cols:
        signed = "hm_beta,0" if "hm_beta" in cols else "hm_odds_ratio,1"
        kept = {
            "hm_rsid", "hm_effect_allele", "hm_other_allele",
            "p_value", signed.split(",")[0],
        }
        # LDSC auto-detects un-prefixed names; tell it explicitly which
        # ambiguous columns to skip.
        candidates = {
            "rsid", "variant_id",
            "effect_allele", "other_allele",
            "beta", "odds_ratio",
            "standard_error",
            "hm_beta", "hm_odds_ratio",
        }
        ignore = sorted((candidates & cols) - kept)
        args = [
            "--snp", "hm_rsid",
            "--a1",  "hm_effect_allele",
            "--a2",  "hm_other_allele",
            "--p",   "p_value",
            "--signed-sumstats", signed,
        ]
        if ignore:
            args += ["--ignore", ",".join(ignore)]
        return args

    # ----- GWAS Catalog "new harmonised" / unharmonised ---------------
    # SNP: prefer rsid; fall back to variant_id (often itself an rsID).
    snp = "rsid" if "rsid" in cols else "variant_id"
    # Signed effect: prefer beta over odds_ratio.
    signed = "beta,0" if "beta" in cols else "odds_ratio,1"

    ignore = []
    if snp == "rsid" and "variant_id" in cols:
        ignore.append("variant_id")
    elif snp == "variant_id" and "rsid" in cols:
        ignore.append("rsid")
    if signed.startswith("beta") and "odds_ratio" in cols:
        ignore.append("odds_ratio")
    elif signed.startswith("odds_ratio") and "beta" in cols:
        ignore.append("beta")

    args = [
        "--snp", snp,
        "--a1",  "effect_allele",
        "--a2",  "other_allele",
        "--p",   "p_value",
        "--signed-sumstats", signed,
    ]
    if ignore:
        args += ["--ignore", ",".join(ignore)]
    return args


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gwas",       required=True)
    parser.add_argument("--study",      required=True)
    parser.add_argument("--study-info", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--munge",      required=True,
                        help="Path to ldsc munge_sumstats.py")
    args = parser.parse_args()

    info = pl.read_csv(args.study_info).filter(pl.col("study") == args.study)
    if info.height == 0:
        raise SystemExit(f"No row in study info for {args.study}")

    n     = info["N"][0]
    n_cas = info["N_cases"][0]
    n_con = info["N_controls"][0]

    cmd = [
        "python", args.munge,
        "--sumstats", args.gwas,
        "--out",      args.out_prefix,
    ]
    cmd += detect_format_args(args.gwas)

    if n is not None:
        cmd += ["--N", str(n)]
    if n_cas is not None:
        cmd += ["--N-cas", str(n_cas)]
    if n_con is not None:
        cmd += ["--N-con", str(n_con)]

    print("Running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
