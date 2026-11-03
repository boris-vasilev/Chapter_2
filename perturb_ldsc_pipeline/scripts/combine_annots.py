#!/usr/bin/env python
"""Combine per-annot binary annotation files for chr 1..22.

For each chromosome reads ``{in_dir}/{condition}_{name}.{chr}.annot.gz``
for every annotation name, renames the ``ANNOT`` column to that name,
horizontally concatenates the result, and writes
``{out_dir}/{condition}_combined.{chr}.annot``.

Annotations whose total SNP count across chr 1..22 is below
``--min-snps`` are dropped before writing the combined files — this
avoids degenerate (all-zero or near-empty) columns that produce NaN
LD scores and singular-matrix failures in h². Pass ``--min-snps 0``
to disable.

If every annotation falls below the threshold the script exits with
an error so the slurm job for that (study, condition) fails loudly
instead of producing a useless empty combined annot.
"""

import argparse
import sys
from pathlib import Path

import polars as pl


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--in-dir",    required=True)
    parser.add_argument("--out-dir",   required=True)
    parser.add_argument("--condition", required=True)
    parser.add_argument("--names",     nargs="+", required=True,
                        help="Annotation names matching {condition}_{name}.{chr}.annot.gz")
    parser.add_argument("--min-snps",  type=int, default=1000,
                        help="Drop annotations whose total SNP count summed "
                             "across chr 1..22 is below this (default 1000). "
                             "0 disables filtering.")
    args = parser.parse_args()

    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- First pass: total SNP count per annotation across all chroms ----
    totals = {name: 0 for name in args.names}
    for chrom in range(1, 23):
        for name in args.names:
            f = in_dir / f"{args.condition}_{name}.{chrom}.annot.gz"
            df = pl.read_csv(f, separator="\t")
            totals[name] += int(df["ANNOT"].sum())

    keep = [n for n in args.names if totals[n] >= args.min_snps]
    drop = [n for n in args.names if totals[n] <  args.min_snps]

    print("Annotation SNP totals across chr 1..22:", flush=True)
    for name in args.names:
        flag = "KEEP" if name in keep else "DROP"
        print(f"  [{flag}] {name}: {totals[name]} SNPs", flush=True)

    if not keep:
        raise SystemExit(
            f"All annotations have < {args.min_snps} SNPs total — refusing "
            f"to write an empty combined annot for {args.condition}."
        )
    if drop:
        print(
            f"\nDropping {len(drop)} sparse annotation(s): {drop}",
            file=sys.stderr, flush=True,
        )

    # ---- Second pass: write combined files using only kept columns ----
    for chrom in range(1, 23):
        dfs = []
        for name in keep:
            f = in_dir / f"{args.condition}_{name}.{chrom}.annot.gz"
            df = pl.read_csv(f, separator="\t").rename({"ANNOT": name})
            dfs.append(df)

        out_file = out_dir / f"{args.condition}_combined.{chrom}.annot"
        pl.concat(dfs, how="horizontal").write_csv(out_file, separator="\t")
        print(f"Wrote {out_file}", flush=True)


if __name__ == "__main__":
    main()
