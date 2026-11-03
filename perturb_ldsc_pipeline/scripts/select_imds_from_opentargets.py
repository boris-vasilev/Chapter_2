#!/usr/bin/env python
"""Select IMDs (immune-mediated diseases) from the OpenTargets `disease` table.

Two modes:

  --inspect
      Print every therapeutic area found in the disease table, along
      with how many diseases sit under it. Use this to pick the EFO
      ID you want to filter on.

  (default — no --inspect)
      Filter the disease table to entries whose `therapeuticAreas`
      array contains --therapeutic-area, and write a CSV of
      candidate IMDs to --out.

Output CSV columns:
  diseaseId         OT disease identifier (== id)
  name              human-readable disease label
  therapeuticAreas  pipe-separated list of therapeutic areas
"""

import argparse
from pathlib import Path

import polars as pl


DEFAULT_AREA = "EFO_0000540"   # immune system disease


def disease_lazy(opentargets_dir):
    return pl.scan_parquet(
        str(Path(opentargets_dir) / "disease" / "*.parquet")
    )


def inspect(opentargets_dir):
    df = (
        disease_lazy(opentargets_dir)
        .select(["id", "name", "therapeuticAreas"])
        .collect()
    )
    counts = (
        df.explode("therapeuticAreas")
        .group_by("therapeuticAreas")
        .agg(pl.len().alias("n_diseases"))
        .sort("n_diseases", descending=True)
    )
    # Look up the therapeutic area name in the same disease table
    counts = counts.join(
        df.select([pl.col("id"), pl.col("name").alias("therapeuticAreaName")]),
        left_on="therapeuticAreas",
        right_on="id",
        how="left",
    ).rename({"therapeuticAreas": "therapeuticAreaId"})

    print(counts.to_pandas().to_string(index=False))


def select(opentargets_dir, therapeutic_area, out):
    df = (
        disease_lazy(opentargets_dir)
        .filter(pl.col("therapeuticAreas").list.contains(therapeutic_area))
        .select([
            pl.col("id").alias("diseaseId"),
            pl.col("name"),
            pl.col("therapeuticAreas").list.join("|").alias("therapeuticAreas"),
        ])
        .sort("diseaseId")
        .collect()
    )
    df.write_csv(out)
    print(f"Selected {df.height} diseases under {therapeutic_area}")
    print(f"Wrote: {out}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--opentargets-dir", required=True)
    parser.add_argument("--inspect", action="store_true",
                        help="List therapeutic areas with disease counts")
    parser.add_argument("--therapeutic-area", default=DEFAULT_AREA,
                        help=f"EFO/MONDO ID to filter on (default {DEFAULT_AREA})")
    parser.add_argument("--out",
                        help="Output CSV (required unless --inspect)")
    args = parser.parse_args()

    if args.inspect:
        inspect(args.opentargets_dir)
    else:
        if not args.out:
            raise SystemExit("--out is required unless --inspect is used")
        select(args.opentargets_dir, args.therapeutic_area, args.out)


if __name__ == "__main__":
    main()
