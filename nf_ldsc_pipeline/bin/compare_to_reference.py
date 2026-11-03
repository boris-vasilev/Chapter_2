#!/usr/bin/env python
"""Compare a pipeline run with existing outputs; exits 1 on any difference.

  --de NEW REF          DE tables: same rows in any order, numbers within --tol
  --gene-sets NEW REF   directories: each file in NEW equals REF as a set
  --results NEW REF     directories: each *.results in NEW equals the same
                        relative path in REF, numbers within --tol

Gene-set files are compared as sets because make_gene_sets_coloc.py writes
them in polars .unique() order, which differs between runs.
"""

import argparse
from pathlib import Path

import numpy as np
import polars as pl


def frames_differ(a, b, tol):
    """Name of the first differing column, or None. Row order must match."""
    if a.columns != b.columns or a.height != b.height:
        return f"shape {a.shape} vs {b.shape}"
    for c in a.columns:
        if a[c].dtype.is_numeric():
            x, y = a[c].to_numpy().astype(float), b[c].to_numpy().astype(float)
            if not np.isclose(x, y, rtol=tol, atol=0, equal_nan=True).all():
                return c
        elif a[c].to_list() != b[c].to_list():
            return c
    return None


def compare_de(new, ref, tol):
    keys = ["culture_condition", "cis_ensg", "trans_ensg"]
    cols = pl.read_csv(ref, n_rows=0).columns
    a, b = (pl.read_csv(f, infer_schema_length=100000).select(cols).sort(keys) for f in (new, ref))
    bad = frames_differ(a, b, tol)
    return [f"DE table: {bad}"] if bad else []


def read_set(path):
    return {line.strip() for line in open(path) if line.strip()}


def compare_gene_sets(new, ref, tol):
    files = sorted(Path(new).glob("*.txt"))
    print(f"{len(files)} gene-set files")
    return [f.name for f in files if read_set(f) != read_set(Path(ref) / f.name)]


def compare_results(new, ref, tol):
    files = sorted(Path(new).rglob("*.results"))
    print(f"{len(files)} .results files")
    bad = []
    for f in files:
        rel = f.relative_to(new)
        a, b = (pl.read_csv(p, separator="\t", null_values=["NA"])
                for p in (f, Path(ref) / rel))
        diff = frames_differ(a, b, tol)
        if diff:
            bad.append(f"{rel}: {diff}")
    return bad


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--de", nargs=2, metavar=("NEW", "REF"))
    p.add_argument("--gene-sets", nargs=2, metavar=("NEW", "REF"))
    p.add_argument("--results", nargs=2, metavar=("NEW", "REF"))
    p.add_argument("--tol", type=float, default=1e-6, help="relative tolerance")
    args = p.parse_args()

    bad = []
    for pair, check in [(args.de, compare_de), (args.gene_sets, compare_gene_sets),
                        (args.results, compare_results)]:
        if pair:
            bad += check(*pair, args.tol)
    print("\n".join(bad) if bad else "OK: no differences")
    raise SystemExit(1 if bad else 0)


if __name__ == "__main__":
    main()
