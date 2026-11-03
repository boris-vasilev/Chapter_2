#!/usr/bin/env python
"""Stage 1 matching check redone with the current (k-NN) matcher.

matching_check.py was run before the out-degree k-NN change, so its FIXED
draws are the quartile-bin matching. This re-draws, at one condition, the
pilot (OLD) sampler and the current matching.Matcher with the Stage 2 seeds
(perm_rng(42, perm, cond)), and records per draw the regulator / target set
sizes and their overlap with the real sets, plus per-anchor expression and
degrees. It also lists, for each knocked-down real anchor, the knocked-down
pool genes in its (in-degree bin, expression bin) cell and its k-NN
neighbours. Draws 1..n_check are checked against the Stage 2 H4 gene sets.

Outputs in --out-dir:
  matching_check_knn_draws.csv        one row per study x method x draw
  matching_check_knn_anchor_values.csv one row per study x method x draw x anchor
  knn_neighbours.csv                  knocked-down anchor cells and neighbours

  uv run python scripts/matching_check_knn.py \
      --out-dir ../data/permutations/figures/source_data
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from matching import (  # noqa: E402
    build_matcher, convergent_regulators, convergent_targets,
    load_coloc_genes, perm_rng,
)
from matching_check import DATA, SHORT, STUDIES, OldSampler, jaccard, load_pilot_module, read_set  # noqa: E402

STAGE2 = DATA / "permutations/heritability"


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--condition", default="Rest")
    p.add_argument("--n-draws", type=int, default=500)
    p.add_argument("--min-abs-logfc", type=float, default=0.5)
    p.add_argument("--knn-k", type=int, default=20)
    p.add_argument("--n-check", type=int, default=200,
                   help="Check k-NN draws 1..n against the Stage 2 H4 gene sets.")
    args = p.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cond = args.condition

    de_full = pl.read_csv(DATA / "DE_results_strong.csv",
                          columns=["cis_ensg", "trans_ensg", "culture_condition",
                                   "logFC", "base_DESeq2"])
    pilot = load_pilot_module()
    all_coloc = load_coloc_genes(str(DATA / "coloc_dat.tsv"), "CD4", 0.8)

    draws, values, neigh = [], [], []
    for study, (coord, pc_coord, _) in STUDIES.items():
        matcher, edges = build_matcher(
            de_full, cond, args.min_abs_logfc,
            str(DATA / "permutations/inputs/screen_expression.tsv"), all_coloc,
            str(DATA / coord), str(DATA / pc_coord), knn_k=args.knn_k)
        b = matcher.binning
        real = set(load_coloc_genes(str(DATA / "coloc_dat.tsv"), "CD4", 0.8, study))
        real_regs = set(convergent_regulators(edges, real))
        real_tgts = set(convergent_targets(edges, real))
        old = OldSampler(pilot, de_full, study)

        for g in sorted(x for x in real if matcher.is_measured(x) and b.out_deg.get(x, 0) > 0):
            cands, d = matcher._knn_candidates(g, set())
            _, ib, eb = b.key(g)
            for x in matcher.knock_cells.get((ib, eb), []):
                neigh.append({"study": SHORT[study], "anchor": g,
                              "anchor_out_degree": b.out_deg[g], "gene": x,
                              "out_degree": b.out_deg.get(x, 0), "in_cell": True,
                              "is_neighbour": x in cands, "relax_dist": d})
            for x in sorted(set(cands) - set(matcher.knock_cells.get((ib, eb), []))):
                neigh.append({"study": SHORT[study], "anchor": g,
                              "anchor_out_degree": b.out_deg[g], "gene": x,
                              "out_degree": b.out_deg.get(x, 0), "in_cell": False,
                              "is_neighbour": True, "relax_dist": d})

        def record(method, perm, anchors):
            regs = set(convergent_regulators(edges, anchors))
            tgts = set(convergent_targets(edges, anchors))
            draws.append({
                "study": SHORT[study], "method": method, "perm": perm,
                "n_regulators": len(regs), "n_targets": len(tgts),
                "reg_overlap_frac_of_real": len(regs & real_regs) / max(1, len(real_regs)),
                "tgt_overlap_frac_of_real": len(tgts & real_tgts) / max(1, len(real_tgts)),
                "reg_jaccard": jaccard(regs, real_regs),
                "tgt_jaccard": jaccard(tgts, real_tgts),
            })
            for g in anchors:
                values.append({"study": SHORT[study], "method": method, "perm": perm,
                               "gene": g, "out_degree": b.out_deg.get(g, 0),
                               "in_degree": b.in_deg.get(g, 0),
                               "expression": matcher.expr.get(g)})

        record("real", 0, real)
        n_match = 0
        for perm in range(1, args.n_draws + 1):
            record("old", perm, old.draw(perm))
            anchors, _ = matcher.sample(sorted(real), perm_rng(42, perm, cond))
            record("knn", perm, anchors)
            if perm <= args.n_check:
                saved = read_set(STAGE2 / study / f"perm{perm}/gene_sets"
                                 / f"{study}_{cond}_H4_genes.txt")
                n_match += anchors == saved
        print(f"[{SHORT[study]}] k-NN draws 1..{args.n_check} identical to Stage 2 "
              f"anchor sets: {n_match}/{args.n_check}", file=sys.stderr)

    pl.DataFrame(draws).write_csv(out / "matching_check_knn_draws.csv")
    pl.DataFrame(values).write_csv(out / "matching_check_knn_anchor_values.csv")
    pl.DataFrame(neigh).write_csv(out / "knn_neighbours.csv")
    ddf = pl.DataFrame(draws).filter(pl.col("method") != "real")
    with pl.Config(tbl_cols=-1, tbl_width_chars=200, float_precision=3):
        print(ddf.group_by("study", "method").agg(
            pl.col("reg_overlap_frac_of_real").median().alias("reg_ov_med"),
            pl.col("reg_overlap_frac_of_real").quantile(0.025).alias("reg_ov_lo"),
            pl.col("reg_overlap_frac_of_real").quantile(0.975).alias("reg_ov_hi"),
            pl.col("tgt_overlap_frac_of_real").median().alias("tgt_ov_med"),
            pl.col("tgt_overlap_frac_of_real").quantile(0.025).alias("tgt_ov_lo"),
            pl.col("tgt_overlap_frac_of_real").quantile(0.975).alias("tgt_ov_hi"),
            pl.col("n_targets").median().alias("tgt_med")).sort("study", "method"))


if __name__ == "__main__":
    main()
