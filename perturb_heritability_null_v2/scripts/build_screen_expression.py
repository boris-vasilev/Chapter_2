#!/usr/bin/env python
"""Per-condition baseline expression of every gene measured in the CD4 screen.

The matched null needs an expression value for every pool gene AND every real
anchor. DE_results_strong.csv only carries base_DESeq2 for genes that are
significantly DE somewhere, so it cannot supply one for anchors with no edges.
Here we read the DESeq2 baseMean layer of the full DE h5ad (the same quantity
DE_results_strong's base_DESeq2 column is copied from) and take, per gene and
culture condition, the median over all perturbation contrasts of that
condition. baseMean is ~constant within gene x condition (median relative SD
0.7% in DE_results_strong), so the median is a stable per-condition value.

Runs outside the LDSC container (needs h5py):
  uv run python scripts/build_screen_expression.py \
      --h5ad ../data/GWCD4i.DE_stats.h5ad \
      --out ../data/permutations/inputs/screen_expression.tsv

Output columns: gene (Ensembl, version stripped), condition, base_mean.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import h5py
import numpy as np
import polars as pl


def read_categorical(group):
    cats = [c.decode() if isinstance(c, bytes) else c for c in group["categories"][:]]
    return np.asarray(cats, dtype=object)[group["codes"][:]]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--block", type=int, default=2000,
                   help="Rows of the baseMean layer read per block.")
    args = p.parse_args()

    with h5py.File(args.h5ad, "r") as f:
        genes = [g.decode().split(".")[0] for g in f["var/gene_ids"][:]]
        cond = read_categorical(f["obs/culture_condition"])
        layer = f["layers/baseMean"]
        n_obs, n_var = layer.shape
        conditions = sorted(set(cond))
        # Collect the rows of each condition block-wise (the layer is a
        # contiguous dense float64 matrix, so whole-row blocks read fast).
        per_cond = {c: [] for c in conditions}
        for start in range(0, n_obs, args.block):
            stop = min(start + args.block, n_obs)
            block = layer[start:stop, :]
            for c in conditions:
                sel = cond[start:stop] == c
                if sel.any():
                    per_cond[c].append(block[sel])
        rows = []
        for c in conditions:
            mat = np.vstack(per_cond[c])
            med = np.nanmedian(mat, axis=0)
            rows.append(pl.DataFrame({
                "gene": genes,
                "condition": [c] * n_var,
                "base_mean": med,
                "n_contrasts": [mat.shape[0]] * n_var,
            }))
            print(f"{c}: {mat.shape[0]} contrasts, {np.isnan(med).sum()} genes NaN")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pl.concat(rows).write_csv(out, separator="\t")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
