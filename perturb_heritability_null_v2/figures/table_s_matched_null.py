#!/usr/bin/env python
"""Supplementary table: real S-LDSC enrichment and coefficient z against the
matched null for all 5 runs x 3 sets. Null median and 2.5-97.5% range are
computed from the per-perm CSVs; real values and empirical p match
summary_all.csv (asserted in fig_common.run_stats).

  uv run python figures/table_s_matched_null.py
"""

from __future__ import annotations

import polars as pl

from fig_common import COND_NAMES, DISEASE_NAMES, OUT, SET_LABELS, run_stats


def fmt_p(p):
    return f"{p:.3f}"


def main():
    st = run_stats()
    wide = None
    for stat, tag in (("enrichment", "enr"), ("coefficient_z", "coefz")):
        s = (st.filter(pl.col("statistic") == stat)
               .select("trait", "condition", "set", "n_null",
                       pl.col("real").alias(f"{tag}_real"),
                       pl.col("null_median").alias(f"{tag}_null_median"),
                       pl.col("null_q025").alias(f"{tag}_null_q025"),
                       pl.col("null_q975").alias(f"{tag}_null_q975"),
                       pl.col("empirical_p").alias(f"{tag}_empirical_p")))
        wide = s if wide is None else wide.join(s.drop("n_null"),
                                                on=["trait", "condition", "set"],
                                                maintain_order="left")
    wide = wide.with_columns(
        pl.col("trait").replace_strict(DISEASE_NAMES).alias("disease"),
        pl.col("condition").replace_strict(COND_NAMES).alias("condition_label"))
    OUT.mkdir(parents=True, exist_ok=True)
    cols = ["disease", "trait", "condition", "condition_label", "set", "n_null"] + [
        c for c in wide.columns if c.startswith(("enr_", "coefz_"))]
    wide.select(cols).write_csv(OUT / "table_s_matched_null.csv")

    lines = [r"\begin{tabular}{lllrrlrrrlr}", r"\toprule",
             r" & & & \multicolumn{4}{c}{Enrichment} & \multicolumn{4}{c}{Coefficient $z$} \\",
             r"\cmidrule(lr){4-7}\cmidrule(lr){8-11}",
             r"Disease & Condition & Set & Real & Null median & Null 95\% & $p$"
             r" & Real & Null median & Null 95\% & $p$ \\", r"\midrule"]
    prev = None
    for r in wide.iter_rows(named=True):
        dis = r["disease"] if (r["disease"], r["condition"]) != prev else ""
        cond = r["condition_label"] if (r["disease"], r["condition"]) != prev else ""
        if prev is not None and (r["disease"], r["condition"]) != prev:
            lines.append(r"\addlinespace")
        prev = (r["disease"], r["condition"])
        cells = [dis, cond, SET_LABELS[r["set"]].split(" (")[0]]
        for tag in ("enr", "coefz"):
            cells += [f"{r[f'{tag}_real']:.2f}", f"{r[f'{tag}_null_median']:.2f}",
                      f"[{r[f'{tag}_null_q025']:.2f}, {r[f'{tag}_null_q975']:.2f}]",
                      fmt_p(r[f"{tag}_empirical_p"])]
        lines.append(" & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    (OUT / "table_s_matched_null.tex").write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUT / 'table_s_matched_null.csv'} and .tex")


if __name__ == "__main__":
    main()
