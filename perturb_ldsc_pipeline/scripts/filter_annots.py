#!/usr/bin/env python
"""Drop low-SNP-count annotation columns from combined .annot files.

Reads a single chromosome's ``.l2.ldscore.gz`` (typically chr 1) to
identify annotation columns that are degenerate (``var == 0``) or
have a total tagged-SNP count below ``--min-snps``, then rewrites
``{in_prefix}.{chr}.annot`` → ``{out_prefix}.{chr}.annot`` for chr
1..22 keeping only the surviving columns.
"""

import argparse
from pathlib import Path

import polars as pl


def find_good_annotations(ldscore_file, min_snp_count):
    df = pl.read_csv(ldscore_file, separator="\t")
    good, bad = [], []
    for col in df.columns:
        if not col.endswith("L2"):
            continue
        base = col.replace("L2", "")
        var = df[col].var()
        total = df[col].sum()
        if var is None or var == 0 or total < min_snp_count:
            bad.append(base)
        else:
            good.append(base)
    print(f"Keeping ({len(good)}): {good}", flush=True)
    print(f"Dropping ({len(bad)}): {bad}", flush=True)
    return good


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ldscore", required=True,
                        help="Reference l2.ldscore.gz used to detect bad columns")
    parser.add_argument("--in-prefix",  required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--min-snps", type=int, default=1000)
    args = parser.parse_args()

    keep = find_good_annotations(args.ldscore, args.min_snps)
    if not keep:
        raise SystemExit("No annotations to keep after filtering")

    out_dir = Path(args.out_prefix).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    for chrom in range(1, 23):
        fin  = f"{args.in_prefix}.{chrom}.annot"
        fout = f"{args.out_prefix}.{chrom}.annot"
        df = pl.read_csv(fin, separator="\t")
        df.select(keep).write_csv(fout, separator="\t")
        print(f"Wrote {fout}", flush=True)


if __name__ == "__main__":
    main()
