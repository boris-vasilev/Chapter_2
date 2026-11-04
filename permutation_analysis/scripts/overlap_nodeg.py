#!/usr/bin/env python
"""Fig. 5C overlap with and without the >= 75-DEG knockdown filter (no LDSC).

Two edge tables are built from the same source (DE_results_sig.csv, padj < 0.05,
joined to DE_stats.suppl_table.csv on knockdown x condition):
  deg75  on-target significant, n_cells_target >= 50, n_total_de_genes >= 75
  nodeg  the same without the DEG filter
then |z| >= 2 (z = logFC / lfcSE), |logFC| >= 0.5, no self-edges. Before the
|logFC| filter, deg75 must reproduce DE_results_strong.csv exactly.

For each variant x trait x condition: real regulator / target sets from the
coloc anchors, and n_draws matched pseudo-anchor sets (matching.Matcher, k-NN
for knocked-down anchors, seeds perm_rng(42, draw, cond)) expanded with the same
edges. Per draw: fraction of the real set recovered, Jaccard, pseudo set size.

Anchors not measured in the screen have no edges, so their replacement cannot
change the overlap. They are kept in the real sets and left out of the draws,
because the gene-coordinate files their matching needs are not available locally.

Outputs in data/permutations/overlap_nodeg/: edge_counts.csv, overlap_draws.csv,
overlap_summary.csv.

  uv run python scripts/overlap_nodeg.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import polars as pl

# matching.py lives with the pipeline that uses it.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "nf_ldsc_permutation" / "bin"))
from matching import (  # noqa: E402
    Matcher, condition_edges, convergent_regulators, convergent_targets,
    degrees, load_coloc_genes, load_expression, perm_rng,
)

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "data"
OUT = DATA / "permutations/overlap_nodeg"
CONDS = ["Rest", "Stim8hr", "Stim48hr"]
KEY = ["culture_condition", "cis_ensg", "trans_ensg"]
MIN_DEGS = {"deg75": 75, "nodeg": 0}


def knockdown_edges(sig: pl.DataFrame, meta: pl.DataFrame, min_degs: int) -> pl.DataFrame:
    """Edges of knockdowns passing the extract_de_results.py filters, |z| >= 2,
    no self-edges (no |logFC| filter yet)."""
    kd = meta.filter(pl.col("ontarget_significant"), pl.col("n_cells_target") >= 50,
                     pl.col("n_total_de_genes") >= min_degs)
    return (sig.join(kd, on=["cis_ensg", "culture_condition"], how="inner")
               .with_columns((pl.col("logFC") / pl.col("lfcSE")).alias("zscore"))
               .filter(pl.col("zscore").abs() >= 2, pl.col("cis_ensg") != pl.col("trans_ensg")))


def same_rows(a: pl.DataFrame, b: pl.DataFrame) -> bool:
    cols = b.columns
    return a.select(cols).sort(KEY).equals(b.select(cols).sort(KEY))


def frac(pseudo: set, real: set) -> float:
    return len(pseudo & real) / len(real) if real else float("nan")


def jaccard(a: set, b: set) -> float:
    return len(a & b) / len(a | b) if (a or b) else float("nan")


def load_inputs():
    """Suppl metadata, the sig edge table and the raw (pre-|logFC|) edges per variant."""
    meta = pl.read_csv(DATA / "DE_stats.suppl_table.csv", infer_schema_length=None).select(
        pl.col("target_contrast").alias("cis_ensg"), "culture_condition",
        "n_total_de_genes", "n_cells_target", "ontarget_significant")
    sig = pl.read_csv(DATA / "DE_results_sig.csv")
    return meta, sig, {v: knockdown_edges(sig, meta, k) for v, k in MIN_DEGS.items()}


def load_traits():
    """The 18 traits, each trait's anchors, and all coloc genes (excluded from the pool)."""
    traits = pl.read_csv(REPO / "vis/diseases.csv").select("study", "trait", "group")
    coloc = str(DATA / "coloc_dat.tsv")
    real_anchors = {s: set(load_coloc_genes(coloc, "CD4", 0.8, s)) for s in traits["study"]}
    return traits, real_anchors, set(load_coloc_genes(coloc, "CD4", 0.8))


def neighbours(e: pl.DataFrame):
    """trans -> its regulators and cis -> its targets. With min_k = 1, regulator /
    target sets are unions of these."""
    regs_of, tgts_of = {}, {}
    for cis, trans in e.select("cis_ensg", "trans_ensg").iter_rows():
        regs_of.setdefault(trans, set()).add(cis)
        tgts_of.setdefault(cis, set()).add(trans)
    return regs_of, tgts_of


def expand(anchors, nb) -> set:
    return set().union(*(nb.get(g, set()) for g in anchors))


def matched_draws(e, cond, traits, real_anchors, all_coloc, n_draws, knn_k, seed, tag):
    """Real sets and n_draws matched pseudo-anchor draws for every trait on one
    edge set; degrees and matcher come from these edges. Returns (draw rows,
    {study: (real regulators, real targets)}); `tag` is added to every row."""
    out_deg, in_deg = degrees(e)
    expr = load_expression(str(DATA / "screen_expression.tsv"), cond)
    matcher = Matcher([g for g in expr if g not in all_coloc], {True: [], False: []},
                      set(), out_deg, in_deg, expr, knn_k=knn_k)
    regs_of, tgts_of = neighbours(e)
    rows, real_sets = [], {}
    for study, trait, group in traits.iter_rows():
        real = real_anchors[study]
        real_regs, real_tgts = expand(real, regs_of), expand(real, tgts_of)
        assert real_regs == set(convergent_regulators(e, real))
        assert real_tgts == set(convergent_targets(e, real))
        unmeasured_with_edges = [g for g in real if g not in expr and (g in regs_of or g in tgts_of)]
        assert not unmeasured_with_edges, f"unmeasured anchors with edges: {unmeasured_with_edges}"
        real_sets[study] = (real_regs, real_tgts)
        measured = sorted(g for g in real if matcher.is_measured(g))
        base = tag | {"study": study, "trait": trait, "group": group, "condition": cond,
                      "n_anchors": len(real), "n_anchors_drawn": len(measured),
                      "real_regulators": len(real_regs), "real_targets": len(real_tgts)}
        for d in range(1, n_draws + 1):
            pseudo, _ = matcher.sample(measured, perm_rng(seed, d, cond))
            regs, tgts = expand(pseudo, regs_of), expand(pseudo, tgts_of)
            rows.append(base | {
                "draw": d,
                "reg_recovery": frac(regs, real_regs), "tgt_recovery": frac(tgts, real_tgts),
                "reg_jaccard": jaccard(regs, real_regs), "tgt_jaccard": jaccard(tgts, real_tgts),
                "pseudo_regulators": len(regs), "pseudo_targets": len(tgts)})
    return rows, real_sets


def summarise(draws: pl.DataFrame, ids: list) -> pl.DataFrame:
    """One row per ids x set: real size, recovery median and 95% range, Jaccard, pseudo size."""
    return pl.concat([
        draws.group_by(ids, maintain_order=True).agg(
            pl.lit(s).alias("set"),
            pl.col(f"real_{s}").first().alias("real_size"),
            pl.col(f"{p}_recovery").median().alias("recovery_median"),
            pl.col(f"{p}_recovery").quantile(0.025).alias("recovery_q025"),
            pl.col(f"{p}_recovery").quantile(0.975).alias("recovery_q975"),
            pl.col(f"{p}_jaccard").median().alias("jaccard_median"),
            pl.col(f"pseudo_{s}").median().alias("pseudo_size_median"))
        for s, p in [("regulators", "reg"), ("targets", "tgt")]])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--n-draws", type=int, default=500)
    p.add_argument("--knn-k", type=int, default=20)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    # --- 1. edge tables -----------------------------------------------------
    meta, sig, raw = load_inputs()
    strong = pl.read_csv(DATA / "DE_results_strong.csv")
    if not same_rows(raw["deg75"], strong):
        raise SystemExit("deg75 does not reproduce DE_results_strong.csv: stopping")
    print(f"deg75 reproduces DE_results_strong.csv ({strong.height:,} rows)")
    del strong

    # Optional: the same filters applied to DE_results_sig2.csv.
    sig2 = knockdown_edges(pl.read_csv(DATA / "DE_results_sig2.csv"), meta, 0)
    print(f"nodeg vs DE_results_sig2 (same filters): only sig "
          f"{raw['nodeg'].join(sig2, on=KEY, how='anti').height:,}, only sig2 "
          f"{sig2.join(raw['nodeg'], on=KEY, how='anti').height:,}")
    del sig, sig2

    edges = {(v, c): condition_edges(raw[v], c, 0.5) for v in MIN_DEGS for c in CONDS}
    counts = pl.DataFrame([{"variant": v, "condition": c, "knockdowns": e["cis_ensg"].n_unique(),
                            "genes_with_incoming": e["trans_ensg"].n_unique(), "edges": e.height}
                           for (v, c), e in edges.items()])
    counts.write_csv(OUT / "edge_counts.csv")
    print(counts)

    # --- 2. matched pseudo-anchor overlap ------------------------------------
    traits, real_anchors, all_coloc = load_traits()
    draws = []
    for (variant, cond), e in edges.items():
        rows, _ = matched_draws(e, cond, traits, real_anchors, all_coloc, args.n_draws,
                                args.knn_k, args.seed, {"variant": variant})
        draws += rows
        print(f"[{variant} {cond}] done", file=sys.stderr)
    draws = pl.DataFrame(draws).select("study", "trait", "group", "condition", "variant", pl.all()
                                       .exclude("study", "trait", "group", "condition", "variant"))
    draws.write_csv(OUT / "overlap_draws.csv")

    # --- 3. summary: one row per trait x condition x set x variant -----------
    ids = ["study", "trait", "group", "condition", "variant", "n_anchors", "n_anchors_drawn"]
    summarise(draws, ids).sort("study", "condition", "set", "variant").write_csv(OUT / "overlap_summary.csv")


if __name__ == "__main__":
    main()
