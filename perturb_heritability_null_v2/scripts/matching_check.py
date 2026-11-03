#!/usr/bin/env python
"""Stage 1: does the matching actually match? (no S-LDSC)

For each study at one condition, draw N pseudo-anchor sets with
  OLD   - the pilot sampler, imported unchanged from
          perturb_ldsc_permutation/scripts/make_permutation_gene_sets.py with the
          pilot's pool, degree, expression and seeds (checked against the
          pilot's saved perm gene sets);
  FIXED - matching.Matcher (v2),
and score both against the real anchors on the v2 variables in that
condition's edge set (|logFC| >= 0.5): expression (screen baseMean),
in-degree and out-degree, plus the regulator/target sets they expand to and
their overlap with the real sets.

Outputs in --out-dir:
  matching_check_summary.csv   one row per study x method x variable
  matching_check_draws.csv     one row per study x method x draw
  matching_check.png           distributions vs real

  uv run python scripts/matching_check.py --out-dir ../data/permutations/matching_check
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import polars as pl  # noqa: E402
from scipy.stats import ks_2samp  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from matching import (  # noqa: E402
    build_matcher, convergent_regulators, convergent_targets,
    load_coloc_genes, perm_rng,
)

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "data"
PILOT_SCRIPT = REPO / "perturb_ldsc_permutation/scripts/make_permutation_gene_sets.py"
STUDIES = {
    # study: (gene coord, protein-coding coord, pilot out dir for this cond)
    "AITD_Saevarsdottir_32581359_1": ("gene_coord_hg38.txt",
                                      "gene_coord_protein_coding_hg38.txt",
                                      "ldsc_permutation0.5_rest_AITD"),
    "ASTCO_Ferreira_30929738_1": ("gene_coord.txt",
                                  "gene_coord_protein_coding.txt",
                                  "ldsc_permutation0.5_rest"),
}
SHORT = {"AITD_Saevarsdottir_32581359_1": "AITD", "ASTCO_Ferreira_30929738_1": "ASTCO"}
COLORS = {"real": "#222222", "old": "#eb6834", "fixed": "#2a78d6"}


def load_pilot_module():
    spec = importlib.util.spec_from_file_location("pilot_gs", PILOT_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class OldSampler:
    """The pilot's draw, reproduced line for line from its main()."""

    def __init__(self, pilot, de_full, study):
        self.pilot = pilot
        cis = {g for g in de_full["cis_ensg"].unique().to_list() if g}
        trans = {g for g in de_full["trans_ensg"].unique().to_list() if g}
        self.pool = sorted(trans | cis)
        self.real = pilot.load_real_anchors(str(DATA / "coloc_dat.tsv"), study, "CD4", 0.8)
        grp = de_full.group_by("cis_ensg").agg(pl.col("trans_ensg").unique())
        self.c2t = {c: set(t) for c, t in zip(grp["cis_ensg"].to_list(),
                                              grp["trans_ensg"].to_list()) if c}
        agg = de_full.group_by("trans_ensg").agg(pl.col("base_DESeq2").mean().alias("mb"))
        self.mb = dict(zip(agg["trans_ensg"].to_list(), agg["mb"].to_list()))
        self.match_by = ["perturb_degree", "expr_decile"]
        self.strata = pilot.build_match_strata(self.pool, self.c2t, self.mb,
                                               self.match_by, 5)

    def draw(self, perm):
        rng = np.random.default_rng(42 + perm)
        return set(self.pilot.sample_pseudo_anchors(
            sorted(set(self.real)), self.pool, self.strata, self.c2t, self.mb,
            self.match_by, 5, rng))


def read_set(path):
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def jaccard(a, b):
    return len(a & b) / len(a | b) if (a or b) else float("nan")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--condition", default="Rest")
    p.add_argument("--n-draws", type=int, default=500)
    p.add_argument("--min-abs-logfc", type=float, default=0.5)
    p.add_argument("--n-verify-pilot", type=int, default=5,
                   help="Check OLD draws 1..n against the pilot's saved H4 sets.")
    args = p.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    cond = args.condition

    de_full = pl.read_csv(DATA / "DE_results_strong.csv",
                          columns=["cis_ensg", "trans_ensg", "culture_condition",
                                   "logFC", "base_DESeq2"])
    pilot = load_pilot_module()
    all_coloc = load_coloc_genes(str(DATA / "coloc_dat.tsv"), "CD4", 0.8)

    draws, values, summary = [], [], []
    for study, (coord, pc_coord, pilot_dir) in STUDIES.items():
        matcher, edges = build_matcher(
            de_full, cond, args.min_abs_logfc,
            str(DATA / "permutations/inputs/screen_expression.tsv"), all_coloc,
            str(DATA / coord), str(DATA / pc_coord))
        b = matcher.binning
        real = set(load_coloc_genes(str(DATA / "coloc_dat.tsv"), "CD4", 0.8, study))
        real_regs = set(convergent_regulators(edges, real))
        real_tgts = set(convergent_targets(edges, real))
        old = OldSampler(pilot, de_full, study)

        for perm in range(1, args.n_verify_pilot + 1):
            saved = read_set(DATA / pilot_dir / f"perm{perm}/gene_sets"
                             / f"{study}_{cond}_H4_genes.txt")
            if old.draw(perm) != saved:
                raise SystemExit(f"OLD reproduction differs from pilot perm{perm} "
                                 f"for {study}")
        print(f"[{SHORT[study]}] OLD sampler reproduces pilot perms "
              f"1..{args.n_verify_pilot}", file=sys.stderr)

        def record(method, perm, anchors, log=None):
            regs = set(convergent_regulators(edges, anchors))
            tgts = set(convergent_targets(edges, anchors))
            draws.append({
                "study": SHORT[study], "method": method, "perm": perm,
                "n_anchors": len(anchors),
                "n_measured": sum(matcher.is_measured(g) for g in anchors),
                "n_coloc_in_anchors": len(anchors & set(all_coloc)),
                "n_real_in_anchors": len(anchors & real),
                "n_regulators": len(regs), "n_targets": len(tgts),
                "reg_overlap_frac_of_real": len(regs & real_regs) / max(1, len(real_regs)),
                "tgt_overlap_frac_of_real": len(tgts & real_tgts) / max(1, len(real_tgts)),
                "reg_jaccard": jaccard(regs, real_regs),
                "tgt_jaccard": jaccard(tgts, real_tgts),
                "n_fallback": None if log is None else log["n_fallback"],
                "n_unmeasured_draws": None if log is None else log["n_unmeasured"],
            })
            for g in anchors:
                o, i = b.out_deg.get(g, 0), b.in_deg.get(g, 0)
                e = matcher.expr.get(g)
                values.append({"study": SHORT[study], "method": method, "perm": perm,
                               "gene": g, "out_degree": o, "in_degree": i,
                               "expression": e})

        record("real", 0, real)
        for perm in range(1, args.n_draws + 1):
            record("old", perm, old.draw(perm))
            anchors, log = matcher.sample(sorted(real), perm_rng(42, perm, cond))
            record("fixed", perm, anchors, log)
        print(f"[{SHORT[study]}] drew {args.n_draws} OLD + FIXED sets", file=sys.stderr)

    ddf = pl.DataFrame(draws)
    vdf = pl.DataFrame(values)
    ddf.write_csv(out / "matching_check_draws.csv")

    # Per-variable summary: anchor-level variables (expression over measured
    # anchors only; degrees over all anchors, unmeasured = 0) and set sizes.
    for st in ddf["study"].unique(maintain_order=True).to_list():
        rv = vdf.filter((pl.col("study") == st) & (pl.col("method") == "real"))
        rd = ddf.filter((pl.col("study") == st) & (pl.col("method") == "real")).row(0, named=True)
        for method in ("old", "fixed"):
            mv = vdf.filter((pl.col("study") == st) & (pl.col("method") == method))
            md = ddf.filter((pl.col("study") == st) & (pl.col("method") == method))
            for var in ("expression", "in_degree", "out_degree"):
                r = rv[var].drop_nulls().to_numpy().astype(float)
                pooled = mv[var].drop_nulls().to_numpy().astype(float)
                per_draw = (mv.drop_nulls(var).group_by("perm")
                              .agg(pl.col(var).median().alias("m"),
                                   (pl.col(var) == 0).mean().alias("z"))
                            )
                meds = per_draw["m"].to_numpy()
                summary.append({
                    "study": st, "method": method, "variable": var,
                    "n_real": r.size,
                    "real_median": float(np.median(r)),
                    "null_median_of_draw_medians": float(np.median(meds)),
                    "null_draw_median_q025": float(np.quantile(meds, 0.025)),
                    "null_draw_median_q975": float(np.quantile(meds, 0.975)),
                    "real_frac_zero": float(np.mean(r == 0)),
                    "null_frac_zero": float(per_draw["z"].mean()),
                    "real_mean": float(r.mean()), "null_mean": float(pooled.mean()),
                    "ks_stat_pooled": float(ks_2samp(r, pooled).statistic),
                    "real_within_95pct_of_draw_medians":
                        bool(np.quantile(meds, 0.025) <= np.median(r) <= np.quantile(meds, 0.975)),
                })
            for var in ("n_regulators", "n_targets", "reg_overlap_frac_of_real",
                        "tgt_overlap_frac_of_real", "reg_jaccard", "tgt_jaccard",
                        "n_measured", "n_coloc_in_anchors", "n_fallback"):
                x = md[var].drop_nulls().to_numpy().astype(float)
                if x.size == 0:
                    continue
                summary.append({
                    "study": st, "method": method, "variable": var,
                    "real_median": float(rd[var]) if rd[var] is not None else None,
                    "null_median_of_draw_medians": float(np.median(x)),
                    "null_draw_median_q025": float(np.quantile(x, 0.025)),
                    "null_draw_median_q975": float(np.quantile(x, 0.975)),
                    "null_mean": float(x.mean()),
                })
    sdf = pl.DataFrame(summary)
    sdf.write_csv(out / "matching_check_summary.csv")

    plot(vdf, ddf, out / "matching_check.png", cond)
    with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=250, float_precision=3):
        print(sdf.select("study", "method", "variable", "real_median",
                         "null_median_of_draw_medians", "null_draw_median_q025",
                         "null_draw_median_q975", "real_frac_zero", "null_frac_zero",
                         "ks_stat_pooled"))


def ecdf(ax, x, color, label, lw=2.0, ls="-"):
    x = np.sort(np.asarray(x, dtype=float))
    ax.step(x, np.arange(1, x.size + 1) / x.size, where="post",
            color=color, lw=lw, ls=ls, label=label)


def plot(vdf, ddf, path, cond):
    studies = ddf["study"].unique(maintain_order=True).to_list()
    cols = [("expression", "Expression (screen baseMean)", True),
            ("in_degree", "In-degree (knockdowns hitting gene)", True),
            ("out_degree", "Out-degree (genes hit by knockdown)", True),
            ("n_regulators", "Regulator set size", False),
            ("n_targets", "Target set size", False)]
    fig, axes = plt.subplots(len(studies), len(cols),
                             figsize=(3.6 * len(cols), 3.0 * len(studies)),
                             squeeze=False)
    for r, st in enumerate(studies):
        for c, (var, label, anchor_level) in enumerate(cols):
            ax = axes[r][c]
            if anchor_level:
                for m in ("old", "fixed", "real"):
                    x = vdf.filter((pl.col("study") == st) & (pl.col("method") == m))[var]
                    x = np.log1p(x.drop_nulls().to_numpy().astype(float))
                    ecdf(ax, x, COLORS[m], m, lw=2.6 if m == "real" else 1.6,
                         ls="--" if m == "real" else "-")
                ax.set_xlabel(f"log1p({label.split(' (')[0].lower()})", fontsize=8)
                ax.set_ylabel("ECDF of anchors" if c == 0 else "", fontsize=8)
            else:
                sub = ddf.filter(pl.col("study") == st)
                vals = sub.filter(pl.col("method") != "real")[var].to_numpy()
                bins = np.linspace(vals.min(), vals.max(), 31)
                for m in ("old", "fixed"):
                    ax.hist(sub.filter(pl.col("method") == m)[var].to_numpy(),
                            bins=bins, color=COLORS[m], alpha=0.55, label=m,
                            edgecolor="white", linewidth=0.4)
                rv = sub.filter(pl.col("method") == "real")[var].item(0)
                ax.axvline(rv, color=COLORS["real"], lw=2, ls="--", label="real")
                ax.set_xlabel(label, fontsize=8)
                ax.set_ylabel("draws" if c == 3 else "", fontsize=8)
            ax.tick_params(labelsize=7)
            for s in ("top", "right"):
                ax.spines[s].set_visible(False)
            ax.grid(alpha=0.25, lw=0.5)
            if r == 0:
                ax.set_title(label, fontsize=9)
            if c == 0:
                ax.annotate(st, xy=(-0.32, 0.5), xycoords="axes fraction",
                            rotation=90, va="center", ha="center",
                            fontsize=11, fontweight="bold")
    axes[0][0].legend(fontsize=8, frameon=False, loc="lower right")
    axes[0][3].legend(fontsize=8, frameon=False)
    fig.suptitle(f"Pseudo-anchor matching vs real anchors — {cond}, |logFC| ≥ 0.5 edges "
                 f"(500 draws; expression over measured anchors only)", fontsize=11)
    fig.tight_layout(rect=(0.02, 0, 1, 0.95))
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {path}", file=sys.stderr)


if __name__ == "__main__":
    main()
