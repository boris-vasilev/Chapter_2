#!/usr/bin/env python
"""Stage 1 re-check after the out-degree k-NN change: target-set sizes only.

For AITD and ASTCO at one condition, draw N pseudo-anchor sets with the
current matching.Matcher and compare the real target-set size with the null
median and 2.5-97.5% range. Also reports the k-NN relaxation log.

  uv run python scripts/target_size_check.py \
      --out ../data/permutations/matching_check/target_size_check_knn.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from matching import (  # noqa: E402
    build_matcher, convergent_targets, load_coloc_genes, perm_rng,
)

DATA = Path(__file__).resolve().parents[2] / "data"
STUDIES = {
    "AITD_Saevarsdottir_32581359_1": ("AITD", "gene_coord_hg38.txt",
                                      "gene_coord_protein_coding_hg38.txt"),
    "ASTCO_Ferreira_30929738_1": ("ASTCO", "gene_coord.txt",
                                  "gene_coord_protein_coding.txt"),
}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=True)
    p.add_argument("--condition", default="Rest")
    p.add_argument("--n-draws", type=int, default=500)
    p.add_argument("--min-abs-logfc", type=float, default=0.5)
    p.add_argument("--knn-k", type=int, default=20)
    args = p.parse_args()
    cond = args.condition

    de = pl.read_csv(DATA / "DE_results_strong.csv",
                     columns=["cis_ensg", "trans_ensg", "culture_condition", "logFC"])
    coloc = str(DATA / "coloc_dat.tsv")
    all_coloc = load_coloc_genes(coloc, "CD4", 0.8)
    rows = []
    for study, (short, coord, pc) in STUDIES.items():
        m, edges = build_matcher(de, cond, args.min_abs_logfc,
                                 str(DATA / "permutations/inputs/screen_expression.tsv"),
                                 all_coloc, str(DATA / coord), str(DATA / pc),
                                 knn_k=args.knn_k)
        real = load_coloc_genes(coloc, "CD4", 0.8, study)
        knocked = [g for g in real if m.binning.out_deg.get(g, 0) > 0]
        sizes, relaxed, max_d, n_fb = [], 0, 0, 0
        for perm in range(1, args.n_draws + 1):
            anchors, log = m.sample(real, perm_rng(42, perm, cond))
            sizes.append(len(convergent_targets(edges, anchors)))
            relaxed += log["n_knn_relaxed"]
            n_fb += log["n_fallback"]
            max_d = max([max_d] + log["knn_relax_dist"])
        sizes = np.array(sizes)
        real_size = len(convergent_targets(edges, real))
        lo, hi = np.quantile(sizes, [0.025, 0.975])
        rows.append({
            "study": short, "condition": cond, "knn_k": args.knn_k,
            "n_draws": args.n_draws,
            "knocked_down_anchors": ";".join(
                f"{g}:{m.binning.out_deg[g]}" for g in sorted(knocked)),
            "real_targets": real_size,
            "null_median": float(np.median(sizes)),
            "null_q025": float(lo), "null_q975": float(hi),
            "null_frac_ge_real": float(np.mean(sizes >= real_size)),
            "real_within_95": bool(lo <= real_size <= hi),
            "knn_draws_relaxed": relaxed,
            "knn_max_relax_dist": max_d,
            "n_fallback_total": n_fb,
        })
    df = pl.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(args.out)
    with pl.Config(tbl_cols=-1, tbl_width_chars=250):
        print(df)


if __name__ == "__main__":
    main()
