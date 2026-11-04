#!/usr/bin/env python
"""Hub-contribution curve: how much of the regulator result rides on hub knockdowns?

For variants nodeg and deg75 (overlap_nodeg.py edge tables: |z| >= 2,
|logFC| >= 0.5, no self-edges), per condition, knockdowns are ranked by
out-degree in that variant's graph (ties by Ensembl id) and the top k are
excluded as regulators (k = 0, 10, 25, ..., 1000; all of them if fewer exist).
Exclusion applies to regulator sets only: a removed knockdown cannot be a
regulator, but anchors keep their outgoing edges, so target sets do not change.
Pseudo-anchor sets follow the same rule. The graph itself is unchanged, so
matching (degrees, draws) is done once per variant x condition on the full
graph and every cut reuses the same 500 draws.

Per variant x k x condition x trait: real regulator / target set sizes; anchors
contributing (>= 1 remaining regulator; knockdown with >= 1 target); share of
the k = 0 real regulator set among the removed knockdowns; pseudo-anchor
recovery of the real sets (median, 95% range). Per variant x k x condition
(repeated on each row): median pairwise Jaccard of real regulator sets for
  IMD-IMD pairs; IMD-IMD without same-study pairs (SAME_STUDY); the same pairs
  with anchors shared by the two traits removed from both anchor sets before
  expansion; all other pairs; other pairs with shared anchors removed.
Pairs whose two sets are both empty are skipped.

Writes data/permutations/hub_curve/hub_curve_summary.csv and hub_curve.{pdf,png};
prints the top 10 nodeg knockdowns, the out-degree at each cut, the medians
across IMDs, and (nodeg) the first cut at which median regulator recovery falls
below half its k = 0 value.

  uv run python scripts/hub_curve.py
"""

from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parent), str(HERE.parents[2] / "vis")]
from overlap_nodeg import (CONDS, DATA, MIN_DEGS, condition_edges, expand, frac,  # noqa: E402
                           jaccard, load_inputs, load_traits, neighbours,
                           pseudo_anchor_sets, screen_matcher)
from common import COND_COLORS, COND_LABEL, PRINT_W, PT, print_style, save  # noqa: E402

OUT = DATA / "permutations/hub_curve"
CUTS = [0, 10, 25, 50, 100, 250, 500, 1000]
VARIANT_LS = {"nodeg": "-", "deg75": "--"}
VARIANT_LABEL = {"nodeg": "No DEG filter", "deg75": "≥ 75 DEGs"}
# GWAS from the same study (shared cases or controls), excluded as IMD-IMD pairs.
SAME_STUDY = [{"IBD_DeLange_28067908_1", "CD_DeLange_28067908_1", "UC_DeLange_28067908_1"},
              {"RA_Ishigaki_36333501_1", "RA_Ishigaki_36333501_2"},
              {"ASTCO_Ferreira_30929738_1", "ASTAO_Ferreira_30929738_1"}]
JACCARDS = ["imd_imd", "imd_imd_no_same_study", "imd_imd_no_same_study_no_shared",
            "other", "other_no_shared"]


def same_study(a: str, b: str) -> bool:
    return any(a in fam and b in fam for fam in SAME_STUDY)


def pair_jaccards(anchors: dict, group: dict, regs_of: dict, removed: set) -> dict:
    """Median pairwise Jaccard of real regulator sets (removed knockdowns excluded)."""
    regs = lambda a: expand(a, regs_of) - removed  # noqa: E731
    full = {s: regs(a) for s, a in anchors.items()}
    js = {k: [] for k in JACCARDS}
    for a, b in itertools.combinations(sorted(anchors), 2):
        shared = anchors[a] & anchors[b]
        # Regulator sets from the anchors the two traits do not share.
        own = (regs(anchors[a] - shared), regs(anchors[b] - shared)) if shared else (full[a], full[b])
        pair = (full[a], full[b])
        if group[a] == "IMD" and group[b] == "IMD":
            cases = [("imd_imd", pair)]
            if not same_study(a, b):
                cases += [("imd_imd_no_same_study", pair), ("imd_imd_no_same_study_no_shared", own)]
        else:
            cases = [("other", pair), ("other_no_shared", own)]
        for key, (ra, rb) in cases:
            if ra or rb:
                js[key].append(jaccard(ra, rb))
    return {f"jaccard_{k}": float(np.median(v)) if v else None for k, v in js.items()}


def quantile(v, p):
    return float(np.quantile(v, p)) if len(v) else float("nan")


def run_variant_condition(variant, cond, e, traits, real_anchors, all_coloc, group, args):
    """Summary rows for one variant x condition over all cuts, and the knockdown ranking."""
    ranked = (e.group_by("cis_ensg").len("out_degree")
               .sort(["out_degree", "cis_ensg"], descending=[True, False]))
    regs_of, tgts_of = neighbours(e)
    matcher, _ = screen_matcher(e, cond, all_coloc, args.knn_k)
    # Per trait: full real sets and the draws' full regulator / target sets.
    per_trait = {}
    for study in traits["study"]:
        real = real_anchors[study]
        measured, draws = pseudo_anchor_sets(matcher, real, cond, args.n_draws, args.seed)
        per_trait[study] = (expand(real, regs_of), expand(real, tgts_of), measured,
                            [expand(p, regs_of) for p in draws], [expand(p, tgts_of) for p in draws])
    rows = []
    for k in CUTS:
        n = min(k, ranked.height)
        removed = set(ranked["cis_ensg"].head(n))
        cut = {"variant": variant, "condition": cond, "cut": k, "n_removed": n,
               "out_degree_at_cut": ranked["out_degree"][n - 1] if n else None}
        cut |= pair_jaccards(real_anchors, group, regs_of, removed)
        for study, trait, grp in traits.iter_rows():
            regs0, tgts, measured, pregs, ptgts = per_trait[study]
            regs, a = regs0 - removed, real_anchors[study]
            base = cut | {"study": study, "trait": trait, "group": grp,
                          "n_anchors": len(a), "n_anchors_drawn": len(measured)}
            sets = [("regulators", regs, [p - removed for p in pregs],
                     sum(bool(regs_of.get(g, set()) - removed) for g in a),
                     len(regs0 & removed) / len(regs0) if regs0 else None),
                    ("targets", tgts, ptgts, sum(g in tgts_of for g in a), None)]
            for s, real_set, pseudo, contributing, share in sets:
                rec = [frac(p, real_set) for p in pseudo] if real_set else []
                rows.append(base | {
                    "set": s, "real_size": len(real_set), "anchors_contributing": contributing,
                    "share_k0_regulators_removed": share,
                    "recovery_median": float(np.median(rec)) if rec else float("nan"),
                    "recovery_q025": quantile(rec, 0.025), "recovery_q975": quantile(rec, 0.975),
                    "pseudo_size_median": float(np.median([len(p) for p in pseudo]))})
    return rows, ranked


def imd_medians(summary: pl.DataFrame) -> pl.DataFrame:
    """Median across IMDs per variant x cut x condition (NaN recovery = empty real set, dropped)."""
    per_set = lambda s, c: pl.col(c).filter(pl.col("set") == s).drop_nans().median()  # noqa: E731
    return (summary.filter(pl.col("group") == "IMD")
            .group_by("variant", "condition", "cut", "n_removed", "out_degree_at_cut").agg(
                per_set("regulators", "real_size").alias("reg_size"),
                per_set("targets", "real_size").alias("tgt_size"),
                per_set("regulators", "anchors_contributing").alias("anchors_reg"),
                per_set("targets", "anchors_contributing").alias("anchors_tgt"),
                per_set("regulators", "share_k0_regulators_removed").alias("reg_removed"),
                per_set("regulators", "recovery_median").alias("reg_recovery"),
                per_set("targets", "recovery_median").alias("tgt_recovery"),
                *[pl.col(f"jaccard_{k}").first().alias(f"j_{k}") for k in JACCARDS])
            .sort(pl.col("variant").replace_strict({"nodeg": 0, "deg75": 1}),
                  pl.col("condition").replace_strict({c: i for i, c in enumerate(CONDS)}), "cut"))


def plot(med: pl.DataFrame) -> None:
    """Regulator recovery and three Jaccard panels vs knockdowns excluded."""
    panels = [("reg_recovery", "Regulator recovery"),
              ("j_imd_imd", "IMD–IMD regulator Jaccard"),
              ("j_imd_imd_no_same_study", "IMD–IMD, same-study pairs excluded"),
              ("j_other", "Other pairs")]
    x = np.arange(len(CUTS))  # cuts spaced evenly, labelled with their value
    with print_style():
        fig, axes = plt.subplots(2, 2, figsize=(PRINT_W, 3.9), constrained_layout=True)
        for ax, (col, title) in zip(axes.flat, panels):
            for v, c in itertools.product(VARIANT_LS, CONDS):
                y = med.filter(pl.col("variant") == v, pl.col("condition") == c).sort("cut")[col]
                ax.plot(x, y, color=COND_COLORS[c], ls=VARIANT_LS[v], marker="o", ms=2.2)
            ax.set_xticks(x, [str(k) for k in CUTS])
            ax.set_title(title, fontsize=PT["title"])
            ax.set_ylim(0, None)
        for ax in axes[1]:
            ax.set_xlabel("Top knockdowns excluded as regulators")
        axes[0, 0].set_ylabel("Median across IMDs")
        for ax in (axes[0, 1], axes[1, 0]):
            ax.set_ylabel("Median over trait pairs")
        jac_max = max(axes.flat[i].get_ylim()[1] for i in (1, 2, 3))
        for i in (1, 2, 3):
            axes.flat[i].set_ylim(0, jac_max)  # Jaccard panels on one scale
        fig.legend(handles=[Line2D([], [], color=COND_COLORS[c], label=COND_LABEL[c]) for c in CONDS]
                   + [Line2D([], [], color="#52514e", ls=VARIANT_LS[v], label=VARIANT_LABEL[v])
                      for v in VARIANT_LS], loc="outside right center")
        save(fig, "hub_curve", OUT)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--n-draws", type=int, default=500)
    p.add_argument("--knn-k", type=int, default=20)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    _, sig, raw = load_inputs()
    symbols = dict(sig.select("cis_ensg", "cis_gene").unique("cis_ensg").iter_rows())
    del sig
    traits, real_anchors, all_coloc = load_traits()
    group = dict(traits.select("study", "group").iter_rows())

    rows, top10 = [], []
    for variant, cond in itertools.product(MIN_DEGS, CONDS):
        e = condition_edges(raw[variant], cond, 0.5)
        r, ranked = run_variant_condition(variant, cond, e, traits, real_anchors, all_coloc, group, args)
        rows += r
        if variant == "nodeg":
            top10.append(ranked.head(10).with_columns(
                pl.lit(cond).alias("condition"),
                pl.col("cis_ensg").replace_strict(symbols, default=None).alias("symbol")))
        print(f"[{variant} {cond}] {ranked.height} knockdowns", file=sys.stderr)

    summary = pl.DataFrame(rows).sort("variant", "condition", "cut", "study", "set")
    summary.write_csv(OUT / "hub_curve_summary.csv")
    med = imd_medians(summary)
    plot(med)

    # First cut where nodeg median regulator recovery is below half its k = 0 value.
    half = (med.filter(pl.col("variant") == "nodeg")
               .with_columns((pl.col("reg_recovery").first().over("condition") / 2).alias("half_k0"))
               .filter(pl.col("reg_recovery") < pl.col("half_k0"))
               .group_by("condition", maintain_order=True).first()
               .select("condition", "half_k0", "cut", "out_degree_at_cut", "reg_recovery"))

    with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=250, float_precision=3):
        print("Top 10 knockdowns by out-degree per condition (nodeg)")
        print(pl.concat(top10).select("condition", "symbol", "cis_ensg", "out_degree"))
        print("Knockdowns removed and out-degree of the last one removed, per cut")
        print(med.pivot(on="variant", index=["condition", "cut"], values=["n_removed", "out_degree_at_cut"]))
        print("Medians across the 12 IMDs (regulators; targets are unchanged by exclusion)")
        print(med.drop("n_removed", "out_degree_at_cut", "tgt_size", "anchors_tgt", "tgt_recovery"))
        print("Targets: median size / anchors / recovery (constant across cuts)")
        print(med.filter(pl.col("cut") == 0).select("variant", "condition", "tgt_size",
                                                    "anchors_tgt", "tgt_recovery"))
        print("nodeg: first cut with median regulator recovery below half its k = 0 value")
        print(half)


if __name__ == "__main__":
    main()
