#!/usr/bin/env python
"""Compile a GWAS_study_info.csv from OpenTargets parquet data.

Filters applied (with row counts printed):
  1. studyType == "gwas"
  2. hasSumstats == True
  3. At least one diseaseId matches --imd-list
  4. Has at least one credible set in credible_set/*.parquet
  (Optional 5.) Primary ldPopulation is in --ancestry-filter

Ancestry is *not* filtered by default — instead three columns are
written so you can curate by hand:
  ancestry       — first (primary) ldPopulation
  n_populations  — how many populations the study spans
  populations    — comma-joined list of all populations

A summary of population distribution is printed at the end. If your
data uses unexpected values (e.g. "nfe" instead of "EUR"), they'll
show up there.

Output columns:
  study,N,N_cases,N_controls,
  ancestry,n_populations,populations,
  build,summarystatsLocation,traitFromSource,pubmedId,diseaseIds
"""

import argparse
from pathlib import Path

import polars as pl

# OpenTargets uses gnomAD population codes; map them onto the
# coarser ancestry buckets that the rest of the pipeline expects
# (EUR / EAS / AFR / SAS / AMR / Multi). Anything not listed here
# falls through unmapped — it'll be obvious in the output.
GNOMAD_TO_LDSC = {
    "nfe": "EUR",  # Non-Finnish European
    "fin": "EUR",  # Finnish European (uses standard EUR LD reference)
    "eas": "EAS",  # East Asian
    "afr": "AFR",
    "amr": "AMR",
    "sas": "SAS",
}


def n(lf):
    return lf.select(pl.len()).collect().item()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--opentargets-dir", required=True)
    parser.add_argument(
        "--imd-list",
        required=True,
        help="CSV with a 'diseaseId' column (output of step 1)",
    )
    parser.add_argument("--out", required=True)
    parser.add_argument(
        "--ancestry-filter",
        nargs="+",
        default=None,
        help="Optional: keep only studies whose first "
        "ldPopulation is in this list (e.g. EUR EAS). "
        "Default: keep all populations.",
    )
    args = parser.parse_args()

    ot_dir = Path(args.opentargets_dir)
    imd_set = set(pl.read_csv(args.imd_list)["diseaseId"].to_list())
    print(f"IMDs from {args.imd_list}: {len(imd_set)}")

    study_lazy = pl.scan_parquet(str(ot_dir / "study" / "*.parquet"))
    print(f"Total studies in OT: {n(study_lazy)}")

    # 1. GWAS
    f1 = study_lazy.filter(pl.col("studyType") == "gwas")
    print(f"  studyType == gwas:                 {n(f1)}")

    # 2. hasSumstats
    f2 = f1.filter(pl.col("hasSumstats"))
    print(f"  hasSumstats == True:               {n(f2)}")

    f2_df = f2.collect()

    # 3. Disease in IMD list
    matching_ids = (
        f2_df.select(["studyId", "diseaseIds"])
        .explode("diseaseIds")
        .filter(pl.col("diseaseIds").is_in(list(imd_set)))["studyId"]
        .unique()
        .to_list()
    )
    f3_df = f2_df.filter(pl.col("studyId").is_in(matching_ids))
    print(f"  at least one IMD diseaseId:        {f3_df.height}")

    # 4. Has credible set
    cs_studies = (
        pl.scan_parquet(str(ot_dir / "credible_set" / "*.parquet"))
        .select("studyId")
        .unique()
        .collect()["studyId"]
        .to_list()
    )
    f4_df = f3_df.filter(pl.col("studyId").is_in(cs_studies))
    print(f"  has >= 1 credible set:             {f4_df.height}")

    # Enrich with ancestry summaries (no hard filter).
    enriched = f4_df.with_columns(
        [
            pl.col("ldPopulationStructure").list.len().alias("n_populations"),
            pl.col("ldPopulationStructure")
            .list.first()
            .struct.field("ldPopulation")
            .alias("ancestry_raw"),
            pl.col("ldPopulationStructure")
            .list.eval(pl.element().struct.field("ldPopulation"))
            .list.join(",")
            .alias("populations"),
        ]
    ).with_columns(
        # Map gnomAD → LDSC ancestry bucket (the column the rest of the
        # pipeline reads). Unmapped codes are left as-is so they're visible.
        pl.col("ancestry_raw")
        .replace(GNOMAD_TO_LDSC)
        .alias("ancestry")
    )

    # Optional ancestry filter (5).
    if args.ancestry_filter:
        before = enriched.height
        enriched = enriched.filter(pl.col("ancestry").is_in(args.ancestry_filter))
        print(
            f"  ancestry in {args.ancestry_filter}: {enriched.height} (from {before})"
        )

    # Build output table.
    out = enriched.select(
        [
            pl.col("studyId").alias("study"),
            pl.col("nSamples").alias("N"),
            pl.col("nCases").alias("N_cases"),
            pl.col("nControls").alias("N_controls"),
            pl.col("ancestry"),  # mapped: EUR / EAS / AFR / ...
            pl.col("ancestry_raw"),  # original gnomAD code (nfe/fin/eas/...)
            pl.col("n_populations"),
            pl.col("populations"),
            pl.lit("GRCh38").alias("build"),
            pl.col("summarystatsLocation"),
            pl.col("traitFromSource"),
            pl.col("pubmedId"),
            pl.col("diseaseIds").list.join("|").alias("diseaseIds"),
        ]
    ).sort("study")

    out.write_csv(args.out)
    print(f"\nWrote {out.height} studies to {args.out}")

    # Diagnostic: unique ldPopulation values seen in the output.
    all_pops = (
        out.select(pl.col("populations").str.split(",").alias("p"))
        .explode("p")["p"]
        .unique()
        .drop_nulls()
        .to_list()
    )
    print(f"\nUnique ldPopulation values seen: {sorted(all_pops)}")

    # Diagnostic: top (ancestry, n_populations) combinations.
    pop_dist = (
        out.group_by(["ancestry", "n_populations"])
        .agg(pl.len().alias("n_studies"))
        .sort("n_studies", descending=True)
        .head(20)
    )
    print(f"\nTop (ancestry, n_populations) combinations:")
    print(pop_dist.to_pandas().to_string(index=False))


if __name__ == "__main__":
    main()
