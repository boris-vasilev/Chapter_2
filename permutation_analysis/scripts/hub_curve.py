#!/usr/bin/env python
"""Hub-contribution curve: how much of the coloc result rides on hub knockdowns?

For variants nodeg and deg75 (the overlap_nodeg.py edge tables: |z| >= 2,
|logFC| >= 0.5, no self-edges), per condition, knockdowns are ranked by
out-degree in that variant's graph (ties by Ensembl id) and the top q% are
removed together with all their out-edges (q = 0, 1, 2, 5, 10, 20). On each
pruned graph, per trait:
  - real regulator / target set sizes, and anchors contributing to each
    (anchors with >= 1 regulator; anchors that are knockdowns with >= 1 target);
  - share of the q = 0 real regulator set that is among the removed knockdowns;
  - matched pseudo-anchor recovery of the real sets (500 draws, corrected matching
    with degrees recomputed on the pruned graph; overlap_nodeg.matched_draws);
and, per variant x q x condition, the median pairwise Jaccard of real regulator
sets across traits, IMD-IMD pairs vs all other pairs (pairs with two empty sets
are skipped). Those two Jaccard columns repeat on every row of their group.

Writes data/permutations/hub_curve/hub_curve_summary.csv and hub_curve.{pdf,png};
prints the top 10 nodeg knockdowns by out-degree per condition and the medians
across the 12 IMDs.

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
from overlap_nodeg import (CONDS, DATA, MIN_DEGS, condition_edges, jaccard,  # noqa: E402
                           load_inputs, load_traits, matched_draws, summarise)
from common import COND_COLORS, COND_LABEL, PRINT_W, PT, print_style, save  # noqa: E402

OUT = DATA / "permutations/hub_curve"
QS = [0, 1, 2, 5, 10, 20]
VARIANT_LS = {"nodeg": "-", "deg75": "--"}
VARIANT_LABEL = {"nodeg": "No DEG filter", "deg75": "≥ 75 DEGs"}


def ranked_knockdowns(e: pl.DataFrame) -> pl.DataFrame:
    """Knockdowns by out-degree, highest first (ties by id, so removal is deterministic)."""
    return (e.group_by("cis_ensg").len("out_degree")
             .sort(["out_degree", "cis_ensg"], descending=[True, False]))


def regulator_jaccard(real_sets: dict, group: dict) -> dict:
    """Median pairwise Jaccard of real regulator sets: IMD-IMD pairs vs other pairs."""
    js = {"imd_imd": [], "other": []}
    for a, b in itertools.combinations(sorted(real_sets), 2):
        ra, rb = real_sets[a][0], real_sets[b][0]
        if ra or rb:
            both_imd = group[a] == "IMD" and group[b] == "IMD"
            js["imd_imd" if both_imd else "other"].append(jaccard(ra, rb))
    return {f"reg_jaccard_{k}_median": float(np.median(v)) if v else None for k, v in js.items()}


def plot(med: pl.DataFrame, path_stem: str) -> None:
    """Median across IMDs vs q: regulator recovery, target recovery, IMD-IMD Jaccard."""
    panels = [("reg_recovery", "Regulator recovery"), ("tgt_recovery", "Target recovery"),
              ("reg_jaccard_imd_imd", "IMD–IMD regulator Jaccard")]
    x = np.arange(len(QS))  # q values spaced evenly, labelled with their value
    with print_style():
        fig, axes = plt.subplots(1, 3, figsize=(PRINT_W, 2.1), constrained_layout=True)
        for ax, (col, title) in zip(axes, panels):
            for v, c in itertools.product(MIN_DEGS, CONDS):
                y = med.filter(pl.col("variant") == v, pl.col("condition") == c).sort("q")[col]
                ax.plot(x, y, color=COND_COLORS[c], ls=VARIANT_LS[v], marker="o", ms=2.5)
            ax.set_xticks(x, [str(q) for q in QS])
            ax.set(xlabel="Top knockdowns removed (%)", ylim=(0, None))
            ax.set_title(title, fontsize=PT["title"])
        axes[0].set_ylabel("Median across IMDs")
        fig.legend(handles=[Line2D([], [], color=COND_COLORS[c], label=COND_LABEL[c]) for c in CONDS]
                   + [Line2D([], [], color="#52514e", ls=VARIANT_LS[v], label=VARIANT_LABEL[v])
                      for v in VARIANT_LS], loc="outside right center")
        save(fig, path_stem, OUT)


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

    draws, extra, top10 = [], [], []
    for variant, cond in itertools.product(MIN_DEGS, CONDS):
        e0 = condition_edges(raw[variant], cond, 0.5)
        ranked = ranked_knockdowns(e0)
        if variant == "nodeg":
            top10.append(ranked.head(10).with_columns(
                pl.lit(cond).alias("condition"),
                pl.col("cis_ensg").replace_strict(symbols, default=None).alias("symbol")))
        regs_q0 = None
        for q in QS:
            n_removed = int(ranked.height * q / 100 + 0.5)
            removed = set(ranked["cis_ensg"].head(n_removed))
            e = e0.filter(~pl.col("cis_ensg").is_in(list(removed)))
            tag = {"variant": variant, "q": q, "n_removed": n_removed}
            rows, real_sets = matched_draws(e, cond, traits, real_anchors, all_coloc,
                                            args.n_draws, args.knn_k, args.seed, tag)
            draws += rows
            regs_q0 = regs_q0 or {s: r for s, (r, _) in real_sets.items()}
            have_regs, knocked = set(e["trans_ensg"]), set(e["cis_ensg"])
            jac = regulator_jaccard(real_sets, group)
            for s in real_sets:
                a = real_anchors[s]
                extra.append(tag | jac | {
                    "study": s, "condition": cond,
                    "anchors_with_regulators": len(a & have_regs),
                    "anchors_with_targets": len(a & knocked),
                    "share_q0_regulators_removed":
                        len(regs_q0[s] & removed) / len(regs_q0[s]) if regs_q0[s] else None})
            print(f"[{variant} {cond} q={q}] removed {n_removed} of {ranked.height}", file=sys.stderr)

    ids = ["study", "trait", "group", "condition", "variant", "q", "n_removed",
           "n_anchors", "n_anchors_drawn"]
    summary = (summarise(pl.DataFrame(draws), ids)
               .join(pl.DataFrame(extra), on=["study", "condition", "variant", "q", "n_removed"])
               .with_columns(pl.when(pl.col("set") == "regulators")
                               .then(pl.col("anchors_with_regulators"))
                               .otherwise(pl.col("anchors_with_targets")).alias("anchors_contributing"))
               .drop("anchors_with_regulators", "anchors_with_targets")
               .sort("variant", "condition", "q", "study", "set"))
    summary.write_csv(OUT / "hub_curve_summary.csv")

    # Medians across the 12 IMDs, one row per variant x q x condition.
    imd = summary.filter(pl.col("group") == "IMD")
    # Traits with an empty real set have NaN recovery; they are left out of the median.
    per_set = lambda s, col: pl.col(col).filter(pl.col("set") == s).drop_nans().median()  # noqa: E731
    n_nonempty = lambda s: (pl.col("real_size").filter(pl.col("set") == s) > 0).sum()  # noqa: E731
    med = (imd.group_by("variant", "q", "condition").agg(
               per_set("regulators", "real_size").alias("reg_size"),
               per_set("targets", "real_size").alias("tgt_size"),
               per_set("regulators", "anchors_contributing").alias("anchors_reg"),
               per_set("targets", "anchors_contributing").alias("anchors_tgt"),
               per_set("regulators", "share_q0_regulators_removed").alias("reg_removed"),
               per_set("regulators", "recovery_median").alias("reg_recovery"),
               per_set("targets", "recovery_median").alias("tgt_recovery"),
               n_nonempty("regulators").alias("n_imd_reg"), n_nonempty("targets").alias("n_imd_tgt"),
               pl.col("reg_jaccard_imd_imd_median").first().alias("reg_jaccard_imd_imd"),
               pl.col("reg_jaccard_other_median").first().alias("reg_jaccard_other"))
           .sort(pl.col("variant").replace_strict({"nodeg": 0, "deg75": 1}),
                 pl.col("condition").replace_strict({c: i for i, c in enumerate(CONDS)}), "q"))
    plot(med, "hub_curve")

    with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=220, float_precision=3):
        print("Top 10 knockdowns by out-degree per condition (nodeg)")
        print(pl.concat(top10).select("condition", "symbol", "cis_ensg", "out_degree"))
        print("Medians across the 12 IMDs")
        print(med)


if __name__ == "__main__":
    main()
