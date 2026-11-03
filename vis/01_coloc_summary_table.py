#!/usr/bin/env python
"""Table 2.1 - GWAS traits, CD4+ eQTL studies and number of coloc genes.

One row per disease (display name): ancestry, GWAS \\cite key (cite_key in
diseases.csv), cases / controls, the number of CD4+ eQTL studies in which the
GWAS signal colocalised with an eQTL (PP.H4 > 0.8; coloc_dat.tsv is already
restricted to CD4, so other cell types never count), and the number of
distinct coloc genes. IMDs first, then controls, each by coloc gene count
(descending). A note under the table lists traits from the same GWAS
(sample_group in diseases.csv).

The eQTL study is the suffix of ``blood_cell`` in coloc_dat.tsv
(e.g. ``Th17_UT`` -> ``UT``), so several CD4 subsets from one study count
once. EQTL_STUDY maps those suffixes to study names (unverified); the
per-disease list goes to table_coloc_eqtl_sources.csv.

The .tex is a complete table float (caption + \\label{tab:traits}); \\input
it where the table belongs and refer to it with \\ref{tab:traits}. Needs
\\usepackage{booktabs,tabularx}. Also writes a CSV.
"""

import argparse
from pathlib import Path

import numpy as np
import polars as pl

from common import (DATA, H4_THRESHOLD, add_variant_args, load_coloc,
                    load_diseases, resolve_variant)

# blood_cell suffix -> eQTL resource. Inferred from the dataset labels in
# coloc_dat.tsv; check against the coloc source before publishing.
EQTL_STUDY = {
    "OneK1K": "OneK1K (Yazar 2022)",
    "Perez": "Perez 2022",
    "Randolph": "Randolph 2021",
    "BP": "BLUEPRINT",
    "Sch": "DICE (Schmiedel 2018)",
    "UT": "ImmuNexUT (Ota 2021)",
    "Bo": "Bossini-Castillo 2022",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coloc", default=DATA / "coloc_dat.tsv")
    ap.add_argument("--h4", type=float, default=H4_THRESHOLD)
    add_variant_args(ap)
    args = resolve_variant(ap.parse_args())

    diseases = load_diseases()
    h4 = (
        load_coloc(args.coloc)
        .filter(pl.col("H4") > args.h4, pl.col("imd").is_in(list(diseases.index)))
        .with_columns(pl.col("blood_cell").str.extract(r"_([^_]+)$").alias("eqtl_key"))
    )

    per_study = h4.group_by("imd").agg(
        pl.col("Gene_id").n_unique().alias("n_coloc_genes"),
        pl.col("eqtl_key").unique().sort().alias("eqtl_keys"),
    ).to_pandas().set_index("imd")

    tab = diseases.join(per_study, how="left")
    tab["n_coloc_genes"] = tab["n_coloc_genes"].fillna(0).astype(int)
    tab["eqtl_source"] = tab["eqtl_keys"].apply(
        lambda ks: "; ".join(EQTL_STUDY.get(k, k) for k in ks)
        if isinstance(ks, (list, np.ndarray)) else "")
    tab["n_eqtl_studies"] = tab["eqtl_keys"].apply(
        lambda ks: len(ks) if isinstance(ks, (list, np.ndarray)) else 0)
    # IMDs first, then controls; each block by coloc genes (descending).
    tab = tab.reset_index().sort_values(["group", "n_coloc_genes", "display_name"],
                                        ascending=[True, False, True])

    cols = {"display_name": "Disease", "group": "Class", "ancestry": "Ancestry",
            "cite_key": "GWAS", "n_cases": "N cases", "n_controls": "N controls",
            "n_eqtl_studies": "CD4+ eQTL studies",
            "n_coloc_genes": f"Coloc genes (PP.H4 > {args.h4})"}
    out = tab[list(cols)].rename(columns=cols)
    out["Class"] = out["Class"].astype(str)
    for c in ["N cases", "N controls", "CD4+ eQTL studies",
              f"Coloc genes (PP.H4 > {args.h4})"]:
        out[c] = out[c].map(lambda v: f"{int(v):,}")

    shared = (tab.dropna(subset=["sample_group"]).groupby("sample_group", sort=False)
              .display_name.apply(list))
    footnote = ("Traits from the same GWAS share samples: " + "; ".join(
        f"{', '.join(names)} (\\cite{{{key}}})" for key, names in shared.items()) + ".")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out.rename(columns={"GWAS": "GWAS cite key"}).to_csv(
        out_dir / "table_coloc_summary.csv", index=False)
    tab[["display_name", "n_eqtl_studies", "eqtl_source"]].rename(
        columns={"display_name": "disease",
                 "eqtl_source": "eqtl_studies_unverified_mapping"}).to_csv(
        out_dir / "table_coloc_eqtl_sources.csv", index=False)

    def tex(s):
        return str(s).replace("&", r"\&").replace("_", r"\_").replace("%", r"\%")

    # LaTeX: a complete table float. The tabularx fits \linewidth by
    # construction (Disease absorbs the slack and wraps if needed); IMD /
    # control is shown as group rows, not a column.
    caption = (
        r"\caption[GWAS traits and colocalised genes]{\textbf{GWAS traits and "
        r"colocalised genes.} For each disease: ancestry, GWAS, numbers of cases "
        r"and controls, the number of CD4$^+$ T-cell eQTL studies in which the "
        r"GWAS signal colocalised with an eQTL (PP.H4 $>$ " + f"{args.h4}" + r"), "
        r"and the number of colocalised genes (anchors) across those studies. "
        r"Immune-mediated diseases are listed first, then non-IMD controls, each "
        r"by number of colocalised genes.}")
    group_title = {"IMD": "Immune-mediated diseases", "non-IMD": "Non-IMD controls"}
    body = out.drop(columns="Class")
    ncol = body.shape[1]
    lines = [
        r"\begin{table}[tbp]",
        r"\centering",
        caption,
        r"\label{tab:traits}",
        r"\vspace{4pt}",
        r"\small",
        r"\setlength{\tabcolsep}{4pt}",
        r"\begin{tabularx}{\linewidth}{@{}>{\raggedright\arraybackslash}X l l r r r r@{}}",
        r"\toprule",
        r"Disease & Ancestry & GWAS & Cases & Controls & CD4$^+$ eQTL & Coloc \\",
        r" & & & & & studies & genes \\",
        r"\midrule",
    ]
    groups = out["Class"].tolist()
    for i, r in enumerate(body.itertuples(index=False)):
        if i == 0 or groups[i] != groups[i - 1]:
            if i:
                lines.append(r"\addlinespace")
            lines.append(rf"\multicolumn{{{ncol}}}{{@{{}}l}}{{\textit{{{group_title[groups[i]]}}}}} \\")
        cells = [tex(v) for v in r]
        cells[2] = f"\\cite{{{r[2]}}}"
        lines.append(" & ".join(cells) + r" \\")
    lines += [
        r"\bottomrule",
        r"\end{tabularx}",
        r"\par\smallskip",
        r"{\footnotesize\raggedright " + footnote + r"\par}",
        r"\end{table}",
    ]
    (out_dir / "table_coloc_summary.tex").write_text("\n".join(lines) + "\n")

    print(out.to_string(index=False))
    print(footnote)
    print(f"Wrote {out_dir}/table_coloc_summary.{{csv,tex}}")


if __name__ == "__main__":
    main()
