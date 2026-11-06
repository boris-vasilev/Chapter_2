#!/usr/bin/env python
"""Supplementary table: matched-null permutation results, deg0 vs deg75.

One row per study x condition x set (regulators, targets): real set size,
one-sided empirical p for S-LDSC enrichment and for the coefficient (tau)
z-score, and the BH q-value of each, for both designs side by side. BH is
applied within each statistic for each design (over all its regulator and
target rows). Reads summary_all.csv from both designs' permutation folders
(common.DESIGNS / perm_dir; not CH2_VARIANT). A design with no permutation
results gets empty columns and a note under the table.

Writes figures2/comparison/table_s_deg_filter.{csv,tex}. As for Table 2.1, the
.tex is a complete table float (caption + \\label{tab:deg_filter_perm}); needs
\\usepackage{booktabs,tabularx}.

  uv run python permutation/table_s_deg_filter.py      (from vis/)
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import COND_LABEL, CONDS, DESIGNS, FIGURES2, load_diseases, perm_dir  # noqa: E402

OUT = FIGURES2 / "comparison"
SETS = ["regulators", "targets"]
SET_LABEL = {"regulators": "Regulators", "targets": "Targets"}
DESIGN_LABEL = {"deg0": "No DEG filter (deg0)", "deg75": r"$\geq$ 75 DEGs (deg75)"}
STATS = {"enrichment": "enrichment", "coefficient_z": "tau_z"}   # summary column -> csv tag


def bh(p):
    """Benjamini-Hochberg q-values (NaN-aware)."""
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        x = p[ok]
        order = np.argsort(x)
        ranked = x[order] * x.size / np.arange(1, x.size + 1)
        q_sorted = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
        q_ok = np.empty_like(x)
        q_ok[order] = q_sorted
        q[ok] = q_ok
    return q


def design_table(design):
    f = perm_dir(DESIGNS[design]) / "summary_all.csv"
    if not f.exists():
        print(f"WARNING: {design}: no permutation results ({f}); columns left empty",
              file=sys.stderr)
        return None
    s = pd.read_csv(f)
    s = s[s["set"].isin(SETS)]
    out = s[["study", "condition", "set"]].copy()
    out[f"{design}_size"] = s["real_size"].astype("Int64")
    for col, tag in STATS.items():
        out[f"{design}_p_{tag}"] = s[f"empirical_p_{col}"]
        out[f"{design}_q_{tag}"] = bh(s[f"empirical_p_{col}"])
    return out


def main():
    diseases = load_diseases()
    keys = ["study", "condition", "set"]
    parts = {d: design_table(d) for d in DESIGNS}
    present = [t for t in parts.values() if t is not None]
    if not present:
        raise SystemExit("No permutation results for any design")
    rows = pd.concat([t[keys] for t in present]).drop_duplicates()
    tab = rows
    for d, t in parts.items():
        if t is None:
            t = pd.DataFrame(columns=keys + [f"{d}_size"] + [
                f"{d}_{k}_{tag}" for tag in STATS.values() for k in ("p", "q")])
        tab = tab.merge(t, on=keys, how="left")
    tab.insert(1, "disease", tab["study"].map(diseases["display_name"]))
    order = {s: i for i, s in enumerate(diseases.index)}
    tab = tab.sort_values(["study", "condition", "set"], key=lambda c: (
        c.map(order) if c.name == "study" else
        c.map(CONDS.index) if c.name == "condition" else c.map(SETS.index))).reset_index(drop=True)
    cols = ["study", "disease", "condition", "set"] + [
        f"{d}_{c}" for d in DESIGNS for c in
        ["size", "p_enrichment", "q_enrichment", "p_tau_z", "q_tau_z"]]
    tab = tab[cols]
    OUT.mkdir(parents=True, exist_ok=True)
    tab.to_csv(OUT / "table_s_deg_filter.csv", index=False)

    def num(v, fmt):
        return "--" if pd.isna(v) else format(v, fmt)

    missing = [d for d, t in parts.items() if t is None]
    caption = (
        r"\caption[Matched-null permutation results with and without the DEG "
        r"filter]{\textbf{Matched-null permutation results with and without the "
        r"DEG filter.} For each trait, condition and expanded set: the size of the "
        r"real set and the one-sided empirical $p$ of the real S-LDSC enrichment "
        r"and coefficient ($\tau$) $z$-score against the matched-anchor null, with "
        r"Benjamini--Hochberg $q$ computed within each statistic and design over "
        r"all regulator and target sets. Left, no minimum number of DEGs per "
        r"knockdown (deg0); right, knockdowns with $\geq$ 75 DEGs (deg75).}")
    ncol = 3 + 5 * len(DESIGNS)
    spec = r"@{}>{\raggedright\arraybackslash}X l l" + " r r r r r" * len(DESIGNS) + "@{}"
    head1 = " & & & " + " & ".join(
        rf"\multicolumn{{5}}{{c}}{{{DESIGN_LABEL[d]}}}" for d in DESIGNS) + r" \\"
    rules1 = "".join(rf"\cmidrule(lr){{{4 + 5 * i}-{8 + 5 * i}}}" for i in range(len(DESIGNS)))
    head2 = " & & & " + " & ".join(
        r" & \multicolumn{2}{c}{Enrichment} & \multicolumn{2}{c}{$\tau$ $z$}"
        for _ in DESIGNS) + r" \\"
    rules2 = "".join(rf"\cmidrule(lr){{{5 + 5 * i}-{6 + 5 * i}}}\cmidrule(lr){{{7 + 5 * i}-{8 + 5 * i}}}"
                     for i in range(len(DESIGNS)))
    head3 = r"Disease & Condition & Set & " + " & ".join(
        r"Size & $p$ & $q$ & $p$ & $q$" for _ in DESIGNS) + r" \\"
    lines = [
        r"\begin{table}[tbp]",
        r"\centering",
        caption,
        r"\label{tab:deg_filter_perm}",
        r"\vspace{4pt}",
        r"\small",
        r"\setlength{\tabcolsep}{3pt}",
        rf"\begin{{tabularx}}{{\linewidth}}{{{spec}}}",
        r"\toprule",
        head1, rules1, head2, rules2, head3,
        r"\midrule",
    ]
    prev = None
    for r in tab.itertuples(index=False):
        key = (r.study, r.condition)
        if prev is not None and key != prev:
            lines.append(r"\addlinespace")
        cells = [r.disease if key != prev else "",
                 COND_LABEL[r.condition] if key != prev else "", SET_LABEL[r.set]]
        prev = key
        for d in DESIGNS:
            cells += [num(getattr(r, f"{d}_size"), "d")] + [
                num(getattr(r, f"{d}_{k}_{tag}"), ".3f")
                for tag in STATS.values() for k in ("p", "q")]
        lines.append(" & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabularx}"]
    if missing:
        lines += [r"\par\smallskip",
                  r"{\footnotesize\raggedright -- : no permutation results for "
                  + ", ".join(missing) + r".\par}"]
    lines.append(r"\end{table}")
    assert all(line.count("&") == ncol - 1 for line in lines
               if line.endswith(r"\\") and "multicolumn" not in line)
    (OUT / "table_s_deg_filter.tex").write_text("\n".join(lines) + "\n")
    print(tab.to_string(index=False))
    print(f"Wrote {OUT}/table_s_deg_filter.{{csv,tex}}")


if __name__ == "__main__":
    main()
