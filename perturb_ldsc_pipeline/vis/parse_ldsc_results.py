#!/usr/bin/env python
"""Collect every {study}/{condition}_combined.results into a tidy CSV.

Reads each LDSC ``.results`` file under --results-dir, extracts the
rows for the module's annotations (which end in ``L2_0``), joins on
GWAS_study_info.csv to bring in disease + trait labels, and writes a
long-format CSV ready for plotting.

Output columns:
  study, condition, annot,
  prop_snps, prop_h2, prop_h2_se,
  enrichment, enrichment_se, enrichment_p,
  coefficient, coefficient_se, coefficient_z,
  ancestry, build, diseaseIds, traitFromSource
"""

import argparse
import sys
from pathlib import Path

import polars as pl


def parse_one(results_file, annot_types):
    df = pl.read_csv(results_file, separator="\t")
    categories = df["Category"].to_list()
    rows = []
    for annot in annot_types:
        # LDSC appends L2_N suffixes; match by prefix to tolerate any N.
        matches = [c for c in categories if c.startswith(annot)]
        if not matches:
            continue
        # If multiple matches, take the first (LDSC typically only writes one).
        r = df.filter(pl.col("Category") == matches[0]).row(0, named=True)
        rows.append(
            {
                "annot": annot,
                "category": matches[0],
                "prop_snps": r["Prop._SNPs"],
                "prop_h2": r["Prop._h2"],
                "prop_h2_se": r["Prop._h2_std_error"],
                "enrichment": r["Enrichment"],
                "enrichment_se": r["Enrichment_std_error"],
                "enrichment_p": r["Enrichment_p"],
                "coefficient": r["Coefficient"],
                "coefficient_se": r["Coefficient_std_error"],
                "coefficient_z": r["Coefficient_z-score"],
            }
        )
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--results-dir", required=True, help="data/ldsc_results/<module>/"
    )
    parser.add_argument("--study-info", required=True)
    parser.add_argument(
        "--annots", nargs="+", default=["regulators_l2g", "targets_l2g", "l2g_genes"]
    )
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rows = []
    for study_dir in sorted(Path(args.results_dir).iterdir()):
        if not study_dir.is_dir():
            continue
        for f in sorted(study_dir.glob("*_combined.results")):
            cond = f.name.replace("_combined.results", "")
            for r in parse_one(f, args.annots):
                rows.append({"study": study_dir.name, "condition": cond, **r})

    if not rows:
        raise SystemExit(
            "\nNo annotation rows matched. Compare the prefixes above to "
            "what you passed via --annots."
        )

    long = pl.DataFrame(rows)
    print(
        f"Parsed {long.height} (study, condition, annot) rows from "
        f"{long['study'].n_unique()} studies."
    )

    info_raw = pl.read_csv(args.study_info)
    wanted = ["study", "ancestry", "build", "diseaseIds", "traitFromSource"]
    available = [c for c in wanted if c in info_raw.columns]
    info = info_raw.select(available)

    out = long.join(info, on="study", how="left").with_columns(
        # Short trait label derived from the study name
        # (e.g. "UC_DeLange_..." → "UC"). Useful when traitFromSource is absent.
        pl.col("study").str.split("_").list.first().alias("trait_short")
    )
    out.write_csv(args.out)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
