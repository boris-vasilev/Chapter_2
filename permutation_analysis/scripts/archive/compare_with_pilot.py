#!/usr/bin/env python
"""Stage 2 comparison: pilot (old matching) vs v2 (fixed matching).

For each pilot combination, set (anchors = H4_genes, regulators, targets) and
statistic (enrichment, coefficient z), one row with, for old and new: real
value, null mean / SD, z, empirical p, number of nulls, real set size and null
set-size median / 95% range. The pilot (500 perms, old matching) is re-read from
its .results files with the same collector as v2, because one pilot
permutation_ldsc_perperm.csv (ASTCO Rest) was written before its sweep finished.
The `old_*_first_n` columns restrict the pilot null to perm indices 1..N_new.
`p_crosses_0.05` flags old and new empirical p on opposite sides of 0.05.

Run after collect_permutation_ldsc.py:
  uv run python scripts/archive/compare_with_pilot.py --config old_snakemake/config.yaml
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import polars as pl
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "nf_ldsc_permutation" / "bin"))
from collect_permutation_ldsc import ANNOTS, SET_NAMES, STATS, collect_run, summarise

# pilot data dir per (study, condition)
PILOT = {
    ("AITD_Saevarsdottir_32581359_1", "Rest"):     "ldsc_permutation0.5_rest_AITD",
    ("AITD_Saevarsdottir_32581359_1", "Stim8hr"):  "ldsc_permutation0.5_stim8_AITD",
    ("AITD_Saevarsdottir_32581359_1", "Stim48hr"): "ldsc_permutation0.5_stim48_AITD",
    ("ASTCO_Ferreira_30929738_1", "Rest"):         "ldsc_permutation0.5_rest",
    ("ASTCO_Ferreira_30929738_1", "Stim48hr"):     "ldsc_permutation0.5",
}
SIZE_COLS = ["real_size", "null_size_median", "null_size_q025", "null_size_q975"]


def stat_cols(srow: dict, stat: str) -> dict:
    return {"real": srow[f"real_{stat}"], "null_mean": srow[f"null_{stat}_mean"],
            "null_sd": srow[f"null_{stat}_sd"], "z": srow[f"z_{stat}"],
            "empirical_p": srow[f"empirical_p_{stat}"], "n_null": srow[f"n_null_{stat}"]}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default="old_snakemake/config.yaml")
    args = p.parse_args()
    cfg = yaml.safe_load(open(args.config))
    repo = Path(args.config).resolve().parents[2]
    root = repo / cfg["out_root"]

    rows = []
    for (study, cond), pilot_dir in PILOT.items():
        new_path = root / f"{study}_{cond}" / "permutation_ldsc_summary.csv"
        if not new_path.exists():
            print(f"skip {study} {cond}: {new_path} missing")
            continue
        new = pl.read_csv(new_path, infer_schema_length=None)
        n_new = new["n_perm_completed"].max()
        old_pp = collect_run(repo / "data" / pilot_dir, study, cond, pilot=True)
        old = summarise(old_pp)
        old_n = summarise(old_pp.filter(pl.col("is_real") | (pl.col("perm_idx") <= n_new)))
        for annot in ANNOTS:
            o, on, nw = (d.filter(pl.col("annot") == annot).row(0, named=True)
                         for d in (old, old_n, new))
            for stat in STATS:
                so, son, sn = stat_cols(o, stat), stat_cols(on, stat), stat_cols(nw, stat)
                row = {"study": study, "condition": cond, "set": SET_NAMES[annot],
                       "statistic": stat}
                for k in so:
                    row[f"old_{k}"], row[f"new_{k}"] = so[k], sn[k]
                row["null_mean_change"] = sn["null_mean"] - so["null_mean"]
                row["p_crosses_0.05"] = (so["empirical_p"] < 0.05) != (sn["empirical_p"] < 0.05)
                row["old_null_mean_first_n"] = son["null_mean"]
                row["old_empirical_p_first_n"] = son["empirical_p"]
                for c in SIZE_COLS:
                    row[f"old_{c}"], row[f"new_{c}"] = o[c], nw[c]
                row["old_n_null_dropped"] = o["n_null_dropped"]
                row["new_n_null_dropped"] = nw["n_null_dropped"]
                rows.append(row)
    out = root / "old_vs_new_comparison.csv"
    df = pl.DataFrame(rows)
    df.write_csv(out)
    with pl.Config(tbl_cols=-1, tbl_rows=-1, tbl_width_chars=250, float_precision=3):
        print(df.select("study", "condition", "set", "statistic", "old_real",
                        "new_real", "old_null_mean", "new_null_mean", "old_z", "new_z",
                        "old_empirical_p", "new_empirical_p", "old_n_null",
                        "p_crosses_0.05"))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
