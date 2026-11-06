#!/usr/bin/env python
"""Collect the v2 S-LDSC permutation null and run the pre-specified tests.

Reads {root}/{study}/perm{P}/ldsc_results/{cond}_combined.results for every
study x condition in --runs (whatever has completed; safe on a partial
sweep) and writes, per study x condition, to {out_dir}/{study}_{cond}/:

  permutation_ldsc_perperm.csv   one row per (perm, annot): prop_snps, prop_h2,
                                 enrichment (+SE, p), coefficient (+SE, z),
                                 gene-set size, matching log
  permutation_ldsc_summary.csv   one row per annot (set): real value, null
                                 mean / SD / median / 95th pct, z, empirical p,
                                 real set size, null set-size median / 95% range
  permutation_null_{cond}.png    null histograms, real value marked

and across runs, to {out_dir}/ (default: --root):

  summary_all.csv                all per-run summaries stacked
  pooled_mean_z.csv              per set type x condition: per-trait z and their
                                 mean (descriptive only; no test)
  matching_log_summary.csv       per study x condition: fallbacks, k-NN
                                 relaxations, unmeasured anchors over null perms

Statistics (for enrichment and coefficient z):
  per set:  z = (obs - null mean) / null SD
            empirical p = (1 + #{null >= obs}) / (N + 1), N = non-missing nulls
  pooled:   mean over traits of the per-set z. With two traits this is
            descriptive, so no null or p-value is attached.

  python bin/collect_permutation_ldsc.py --root data/permutations/<tag> \
      --runs AITD_Saevarsdottir_32581359_1:Rest,Stim8hr ASTCO_Ferreira_30929738_1:Rest
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import polars as pl

ANNOTS = ["H4_genes", "regulators_H4", "targets_H4"]
SET_NAMES = {"H4_genes": "anchors", "regulators_H4": "regulators",
             "targets_H4": "targets"}
STATS = ["enrichment", "coefficient_z"]
_FIELDS = {
    "prop_snps": "Prop._SNPs",
    "prop_h2": "Prop._h2",
    "prop_h2_se": "Prop._h2_std_error",
    "enrichment": "Enrichment",
    "enrichment_se": "Enrichment_std_error",
    "enrichment_p": "Enrichment_p",
    "coefficient": "Coefficient",
    "coefficient_se": "Coefficient_std_error",
    "coefficient_z": "Coefficient_z-score",
}
PERM_RE = re.compile(r"^perm(real|\d+)$")
_LOG_FIELDS = ["n_fallback", "max_fallback_dist", "n_knn", "n_knn_relaxed",
               "max_knn_relax_dist", "n_unmeasured"]


def count_genes(path: Path):
    if not path.exists():
        return None
    with open(path) as f:
        return sum(1 for line in f if line.strip())


def single_kept(results: Path):
    """The one annotation combine_annots kept, read from its .combine.txt, or None."""
    f = results.with_suffix(".combine.txt")
    kept = re.findall(r"\[KEEP\] (\S+):", f.read_text()) if f.exists() else []
    return kept[0] if len(kept) == 1 else None


def parse_results(path: Path):
    df = pl.read_csv(path, separator="\t")
    cats = df["Category"].to_list()
    # ldsc names a lone annotation plain L2_1, so take its name from .combine.txt.
    single = single_kept(path) if "L2_1" in cats else None
    out = {}
    for a in ANNOTS:
        # combine_annots names columns {annot}L2_<k>; match on the prefix.
        match = [c for c in cats if c.startswith(a + "L2")] or (["L2_1"] if a == single else [])
        if match:
            r = df.filter(pl.col("Category") == match[0]).row(0, named=True)
            out[a] = {k: r[col] for k, col in _FIELDS.items()}
    return out


def collect_run(root: Path, study: str, cond: str, pilot: bool = False) -> pl.DataFrame:
    """pilot=True reads the pilot layout {root}/perm*/ldsc_results/{study}/."""
    rows = []
    for perm_dir in sorted((root if pilot else root / study).glob("perm*")):
        m = PERM_RE.match(perm_dir.name)
        res_dir = perm_dir / "ldsc_results" / (study if pilot else "")
        rf = res_dir / f"{cond}_combined.results"
        if not m or not rf.exists():
            continue
        perm = m.group(1)
        gs = perm_dir / "gene_sets"
        mlog = gs / f"{study}_{cond}_matching.tsv"
        log = pl.read_csv(mlog, separator="\t").row(0, named=True) if mlog.exists() else {}
        parsed = parse_results(rf)
        for a in ANNOTS:
            size = count_genes(gs / f"{study}_{cond}_{a}.txt")
            row = {"study": study, "condition": cond, "perm": perm,
                   "is_real": perm == "real",
                   "perm_idx": -1 if perm == "real" else int(perm),
                   "annot": a, "set": SET_NAMES[a],
                   "size": size,
                   "empty": not size,
                   "dropped": a not in parsed}
            row |= {k: log.get(k) for k in _LOG_FIELDS}
            row.update(parsed.get(a, {k: None for k in _FIELDS}))
            rows.append(row)
    schema_over = ({k: pl.Float64 for k in _FIELDS}
                   | {k: pl.Int64 for k in _LOG_FIELDS} | {"size": pl.Int64})
    return pl.DataFrame(rows, schema_overrides=schema_over).sort("annot", "perm_idx")


def emp_p(obs, null):
    null = null[~np.isnan(null)]
    if obs is None or np.isnan(obs) or null.size == 0:
        return float("nan")
    return float((1 + np.sum(null >= obs)) / (null.size + 1))


def summarise(df: pl.DataFrame) -> pl.DataFrame:
    recs = []
    for a in ANNOTS:
        sub = df.filter(pl.col("annot") == a)
        real = sub.filter(pl.col("is_real"))
        null = sub.filter(~pl.col("is_real"))
        sz = null["size"].drop_nulls().to_numpy().astype(float)
        rec = {"study": df["study"][0], "condition": df["condition"][0], "annot": a,
               "set": SET_NAMES[a],
               "n_perm_completed": null.height, "real_present": real.height > 0,
               "real_size": real["size"].item(0) if real.height else None,
               "null_size_mean": float(sz.mean()) if sz.size else float("nan"),
               "null_size_median": float(np.median(sz)) if sz.size else float("nan"),
               "null_size_q025": float(np.quantile(sz, 0.025)) if sz.size else float("nan"),
               "null_size_q975": float(np.quantile(sz, 0.975)) if sz.size else float("nan"),
               "real_prop_snps": real["prop_snps"].item(0) if real.height else None,
               "null_prop_snps_mean": null["prop_snps"].mean(),
               "n_null_dropped": int(null["dropped"].sum()),
               "n_null_empty": int(null["empty"].sum())}
        for s in STATS:
            obs = real[s].item(0) if real.height else None
            nv = null[s].drop_nulls().to_numpy().astype(float)
            mu = float(nv.mean()) if nv.size else float("nan")
            sd = float(nv.std(ddof=1)) if nv.size > 1 else float("nan")
            rec |= {f"real_{s}": obs, f"null_{s}_mean": mu, f"null_{s}_sd": sd,
                    f"null_{s}_median": float(np.median(nv)) if nv.size else float("nan"),
                    f"null_{s}_p95": float(np.quantile(nv, 0.95)) if nv.size else float("nan"),
                    f"z_{s}": (obs - mu) / sd if obs is not None and sd > 0 else float("nan"),
                    f"empirical_p_{s}": emp_p(obs, nv), f"n_null_{s}": int(nv.size)}
        recs.append(rec)
    return pl.DataFrame(recs)


def pooled_mean_z(summary: pl.DataFrame) -> pl.DataFrame:
    """Mean over traits of the per-set z, per set type x condition (descriptive)."""
    recs = []
    for (cond, a), sub in summary.group_by(["condition", "annot"], maintain_order=True):
        for s in STATS:
            ok = sub.filter(pl.col(f"z_{s}").is_not_nan())
            recs.append({"condition": cond, "set": SET_NAMES[a], "statistic": s,
                         "n_traits": ok.height,
                         "traits": ";".join(ok["study"].to_list()),
                         "per_trait_z": ";".join(f"{z:.3f}" for z in ok[f"z_{s}"]),
                         "mean_z": ok[f"z_{s}"].mean() if ok.height else None})
    return pl.DataFrame(recs).sort("statistic", "set", "condition")


def matching_summary(perperm: pl.DataFrame) -> pl.DataFrame:
    """Matching-log counts over the null perms, one row per study x condition."""
    logs = (perperm.filter(~pl.col("is_real"))
            .unique(["study", "condition", "perm_idx"]))
    return (logs.group_by("study", "condition", maintain_order=True)
            .agg(pl.len().alias("n_perm"),
                 (pl.col("n_fallback") > 0).sum().alias("n_perm_with_fallback"),
                 pl.col("n_fallback").sum().alias("n_fallback_total"),
                 pl.col("max_fallback_dist").max().alias("max_fallback_dist"),
                 pl.col("n_knn").first().alias("n_knn_anchors"),
                 (pl.col("n_knn_relaxed") > 0).sum().alias("n_perm_with_knn_relaxed"),
                 pl.col("n_knn_relaxed").sum().alias("n_knn_relaxed_total"),
                 pl.col("max_knn_relax_dist").max().alias("max_knn_relax_dist"),
                 pl.col("n_unmeasured").min().alias("n_unmeasured_min"),
                 pl.col("n_unmeasured").max().alias("n_unmeasured_max"))
            .sort("study", "condition"))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", required=True, help="{study}/perm{P}/... tree to read")
    p.add_argument("--out-dir", help="where summaries go (default: --root)")
    p.add_argument("--runs", nargs="+", required=True, metavar="STUDY:COND,COND")
    p.add_argument("--plot", action="store_true", default=True)
    args = p.parse_args()

    root = Path(args.root)
    out_root = Path(args.out_dir or root)
    runs = {s: c.split(",") for s, c in (r.split(":", 1) for r in args.runs)}

    all_perperm, all_summary = [], []
    for study, conds in runs.items():
        for cond in conds:
            df = collect_run(root, study, cond)
            if df.height == 0:
                print(f"[{study} {cond}] no .results yet")
                continue
            out = out_root / f"{study}_{cond}"
            out.mkdir(parents=True, exist_ok=True)
            df.write_csv(out / "permutation_ldsc_perperm.csv")
            sdf = summarise(df)
            sdf.write_csv(out / "permutation_ldsc_summary.csv")
            all_perperm.append(df)
            all_summary.append(sdf)
            if args.plot:
                try:
                    from plot_permutation_null import plot_condition
                except ImportError:
                    # Moved to Chapter_2/vis/heritability/ (reads the published
                    # results); plots are made there, not in the pipeline.
                    plot_condition = None
                if plot_condition is not None:
                    plot_condition(df, cond, ANNOTS, STATS, 30, study,
                                   out / f"permutation_null_{cond}.png")
            with pl.Config(tbl_cols=-1, tbl_width_chars=200, float_precision=3):
                print(f"[{study} {cond}]")
                print(sdf.select("annot", "n_perm_completed", "real_enrichment",
                                 "null_enrichment_mean", "z_enrichment",
                                 "empirical_p_enrichment", "real_coefficient_z",
                                 "z_coefficient_z", "empirical_p_coefficient_z"))

    if all_perperm:
        summary = pl.concat(all_summary, how="diagonal_relaxed")
        summary.write_csv(out_root / "summary_all.csv")
        pooled = pooled_mean_z(summary)
        pooled.write_csv(out_root / "pooled_mean_z.csv")
        mlog = matching_summary(pl.concat(all_perperm, how="diagonal_relaxed"))
        mlog.write_csv(out_root / "matching_log_summary.csv")
        with pl.Config(tbl_cols=-1, tbl_width_chars=200, float_precision=3):
            print(pooled)
            print(mlog)


if __name__ == "__main__":
    main()
