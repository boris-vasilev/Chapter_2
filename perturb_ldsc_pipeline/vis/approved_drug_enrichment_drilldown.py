#!/usr/bin/env python
"""Drill down into the genes driving an enrichment result.

For one (disease, gene_set, condition) triple from approved_drug_enrichment*.py,
list every gene in the gene set that has an approved or Phase II+ drug for the
disease, together with:
  - the drug names and clinical stages for that gene-disease combo;
  - the coloc (H4) gene(s) it is connected to in the Perturb-seq DE results
    (for targets_H4: cis_gene whose perturbation produced this gene as a DEG;
     for regulators_H4: trans_gene that responded when this gene was perturbed);
  - the direction (logFC sign) of each connection.
"""

import argparse
from pathlib import Path

import polars as pl


APPROVED_STAGES = [
    "APPROVAL",
    "PREAPPROVAL",
    "PHASE_4",
    "PHASE_3",
    "PHASE_2_3",
    "PHASE_2",
]
# Order from strongest evidence to weakest, used to pick `max_stage`.
STAGE_PRIORITY = [
    "APPROVAL",
    "PHASE_4",
    "PREAPPROVAL",
    "PHASE_3",
    "PHASE_2_3",
    "PHASE_2",
]


def expand_subtree(disease_ids: list[str], ontology: pl.DataFrame) -> set[str]:
    matches = ontology.filter(pl.col("id").is_in(disease_ids))
    subtree = set(disease_ids)
    for d in matches.get_column("descendants").to_list():
        if d:
            subtree.update(d)
    return subtree


def load_gene_set(
    gene_sets_dir: Path, studies: list[str], gene_set: str, condition: str
) -> set[str]:
    """Union of gene IDs across `studies` for a specific gene_set/condition pair.

    Gene-set filenames are STUDY_AUTHOR_PMID_NUM_CONDITION_GENESET.txt; we glob
    `{study}_*_{condition}_{gene_set}.txt` to pull exactly the matching files.
    """
    genes: set[str] = set()
    for study in studies:
        for file in gene_sets_dir.glob(f"{study}_*_{condition}_{gene_set}.txt"):
            with open(file) as f:
                for line in f:
                    g = line.strip()
                    if g:
                        genes.add(g)
    return genes


def pick_max_stage(stages: list[str]) -> str:
    for s in STAGE_PRIORITY:
        if s in stages:
            return s
    return ""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-sets-dir", required=True)
    parser.add_argument("--studies", nargs="+", required=True)
    parser.add_argument("--drug-targets", required=True)
    parser.add_argument("--disease-ontology", required=True)
    parser.add_argument("--disease-id", nargs="+", required=True)
    parser.add_argument("--de-results", required=True)
    parser.add_argument(
        "--gene-set",
        required=True,
        choices=["H4_genes", "regulators_H4", "targets_H4"],
    )
    parser.add_argument("--condition", required=True, help="e.g. Rest, Stim8hr, Stim48hr")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # 1. Disease subtree.
    ontology = (
        pl.scan_parquet(args.disease_ontology)
        .select(["id", "descendants"])
        .collect()
    )
    subtree = expand_subtree(args.disease_id, ontology)

    # 2. Approved drug-evidence rows for this disease (or any descendant).
    drug_evidence = (
        pl.scan_parquet(args.drug_targets)
        .select(["targetId", "diseaseId", "clinicalStage", "drugFromSource"])
        .filter(pl.col("diseaseId").is_in(subtree))
        .filter(pl.col("clinicalStage").is_in(APPROVED_STAGES))
        .collect()
    )
    approved = set(drug_evidence.get_column("targetId").unique().to_list())

    # 3. Gene set we're drilling into + H4 coloc genes (for connection lookup).
    gene_sets_dir = Path(args.gene_sets_dir)
    chosen_set = load_gene_set(
        gene_sets_dir, args.studies, args.gene_set, args.condition
    )
    h4_genes = load_gene_set(gene_sets_dir, args.studies, "H4_genes", args.condition)

    hits = chosen_set & approved
    print(
        f"{args.gene_set} | {args.condition} | studies={args.studies}: "
        f"|set|={len(chosen_set)}, |approved-for-disease|={len(approved)}, "
        f"|hits|={len(hits)}"
    )
    if not hits:
        raise SystemExit("No hits to drill into.")

    # 4. DE results for the matched condition, used for coloc-gene linkage.
    de = pl.read_csv(args.de_results).filter(
        pl.col("culture_condition") == args.condition
    )

    # 5. Per-hit rows.
    rows = []
    for ensg in sorted(hits):
        # Drugs: deduplicate (drug, stage) pairs for this gene.
        drug_rows = drug_evidence.filter(pl.col("targetId") == ensg)
        pairs = (
            drug_rows.select(["drugFromSource", "clinicalStage"])
            .drop_nulls()
            .unique()
            .iter_rows()
        )
        drug_pairs = sorted(set(pairs))
        drugs_str = "; ".join(f"{d} ({s})" for d, s in drug_pairs)
        max_stage = pick_max_stage(
            drug_rows.get_column("clinicalStage").unique().to_list()
        )

        # Coloc-gene linkage via DE results, gated on H4 membership.
        symbol = ""
        coloc_links: list[tuple[str, float, float]] = []  # (coloc_symbol, logFC, padj)
        if args.gene_set == "targets_H4":
            sub = de.filter(
                (pl.col("trans_ensg") == ensg)
                & (pl.col("cis_ensg").is_in(h4_genes))
            )
            if sub.height > 0:
                symbol = sub.get_column("trans_gene").item(0)
                for r in sub.iter_rows(named=True):
                    coloc_links.append((r["cis_gene"], r["logFC"], r["padj"]))
        elif args.gene_set == "regulators_H4":
            sub = de.filter(
                (pl.col("cis_ensg") == ensg)
                & (pl.col("trans_ensg").is_in(h4_genes))
            )
            if sub.height > 0:
                symbol = sub.get_column("cis_gene").item(0)
                for r in sub.iter_rows(named=True):
                    coloc_links.append((r["trans_gene"], r["logFC"], r["padj"]))
        else:  # H4_genes: the gene itself is the coloc gene
            # Try symbol from either cis_gene or trans_gene anywhere in DE results.
            cis_match = de.filter(pl.col("cis_ensg") == ensg).get_column("cis_gene")
            trans_match = de.filter(pl.col("trans_ensg") == ensg).get_column("trans_gene")
            if cis_match.len() > 0:
                symbol = cis_match.item(0)
            elif trans_match.len() > 0:
                symbol = trans_match.item(0)

        # Fall back to symbol from a different DE row if we still don't have one.
        if not symbol:
            for col_ensg, col_gene in [("cis_ensg", "cis_gene"), ("trans_ensg", "trans_gene")]:
                match = de.filter(pl.col(col_ensg) == ensg).get_column(col_gene)
                if match.len() > 0:
                    symbol = match.item(0)
                    break

        # Collapse repeated (coloc_symbol, sign) pairs; keep extreme logFC.
        by_coloc: dict[str, float] = {}
        for sym, lfc, _ in coloc_links:
            if sym not in by_coloc or abs(lfc) > abs(by_coloc[sym]):
                by_coloc[sym] = lfc
        coloc_str = "; ".join(
            f"{sym} ({'+' if lfc > 0 else '-'}{abs(lfc):.2f})"
            for sym, lfc in sorted(by_coloc.items(), key=lambda kv: -abs(kv[1]))
        )

        rows.append(
            {
                "gene_ensg": ensg,
                "gene_symbol": symbol,
                "n_drugs": len(drug_pairs),
                "max_stage": max_stage,
                "drugs": drugs_str,
                "n_coloc_links": len(by_coloc),
                "coloc_genes": coloc_str,
            }
        )

    df = pl.DataFrame(rows).sort(
        ["n_coloc_links", "gene_symbol"], descending=[True, False]
    )
    df.write_csv(args.output)
    print(f"Wrote {args.output}")
    with pl.Config(
        tbl_rows=-1, tbl_cols=-1, tbl_width_chars=240, fmt_str_lengths=120
    ):
        print(df)


if __name__ == "__main__":
    main()
