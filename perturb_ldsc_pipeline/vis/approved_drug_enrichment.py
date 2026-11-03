#!/usr/bin/env python
"""Enrichment of approved drug targets in coloc gene sets."""

import argparse
import polars as pl
from pathlib import Path
from scipy.stats import fisher_exact

APPROVED_STAGES = [
    "APPROVAL",
    "PREAPPROVAL",
    "PHASE_4",
    "PHASE_3",
    "PHASE_2_3",
    "PHASE_2",
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gene-sets-dir",
        required=True,
        help="e.g. data/gene_sets/coloc or data/gene_sets/coloc_strict",
    )

    parser.add_argument(
        "--studies",
        nargs="+",
        required=True,
        help="Study identifiers to match in gene-set file names (e.g. AST, ASTCO, UC, IBD)",
    )

    parser.add_argument(
        "--drug-targets",
        required=True,
        help=(
            "Drug-target-indication data from Open Targets (dataset: evidence_clinical_precedence). "
            "Expected columns: targetId (Ensembl gene ID), diseaseFromSource (disease name), "
            "clinicalStage, drugFromSource (drug name)."
        ),
    )

    parser.add_argument(
        "--disease-id",
        nargs="+",
        required=True,
        help=(
            "One or more Open Targets disease identifiers (EFO/MONDO/HP IDs, e.g. "
            "MONDO_0004979 for asthma). The test is run against the union of each "
            "ID and all its descendants in the disease ontology."
        ),
    )

    parser.add_argument(
        "--disease-ontology",
        required=True,
        help=(
            "Path to the Open Targets `disease` parquet, used to expand each "
            "--disease-id into its full set of ontology descendants."
        ),
    )

    parser.add_argument(
        "--de-results",
        required=True,
        help=(
            "Path to the Perturb-seq DE results CSV (e.g. DE_results_strong.csv). "
            "Used to build per-set background universes: unique cis_ensg defines "
            "the perturbation universe (background for regulator sets); unique "
            "trans_ensg defines the detected-gene universe (background for target "
            "sets and, as an approximation, for coloc gene sets)."
        ),
    )

    args = parser.parse_args()

    # Expand each --disease-id into itself + all ontology descendants.
    ontology = (
        pl.scan_parquet(args.disease_ontology)
        .filter(pl.col("id").is_in(args.disease_id))
        .select(["id", "descendants"])
        .collect()
    )
    missing = set(args.disease_id) - set(ontology.get_column("id").to_list())
    if missing:
        raise ValueError(f"Disease IDs not found in ontology: {sorted(missing)}")

    disease_subtree = set(args.disease_id)
    for descendants in ontology.get_column("descendants").to_list():
        if descendants:
            disease_subtree.update(descendants)

    # All drug-target evidence rows, kept unfiltered so we can derive a background universe.
    drug_targets_all = (
        pl.scan_parquet(args.drug_targets)
        .select(["targetId", "diseaseId", "clinicalStage"])
        .collect()
    )

    # Genes with an approved / phase-4 drug for the disease of interest (or any descendant).
    approved_for_disease = set(
        drug_targets_all.filter(pl.col("diseaseId").is_in(disease_subtree))
        .filter(pl.col("clinicalStage").is_in(APPROVED_STAGES))
        .get_column("targetId")
        .unique()
        .to_list()
    )

    # Per-set background universes derived from the Perturb-seq DE results.
    de = pl.read_csv(args.de_results, columns=["cis_ensg", "trans_ensg"])
    cis_universe = {g for g in de.get_column("cis_ensg").unique().to_list() if g}
    trans_universe = {g for g in de.get_column("trans_ensg").unique().to_list() if g}

    def universe_for(gene_set: str) -> set[str]:
        # Regulators are drawn from the perturbation library; everything else
        # (targets, coloc H4 genes) is drawn from the detected-gene pool.
        if gene_set.startswith("regulators"):
            return cis_universe
        return trans_universe

    # Load gene-set memberships.
    gene_sets_dir = Path(args.gene_sets_dir)
    rows = []
    for study in args.studies:
        matches = list(gene_sets_dir.glob(f"{study}*"))
        if not matches:
            raise FileNotFoundError(
                f"No gene-set file found for study '{study}' in {gene_sets_dir}"
            )
        for file in matches:
            parts = file.stem.split("_")
            gene_set = "_".join(parts[-2:])  # e.g. targets_H4
            condition = parts[-3]  # e.g. Stim8hr
            with open(file) as f:
                for gene in map(str.strip, f):
                    if gene:
                        rows.append(
                            {
                                "gene": gene,
                                "study": study,
                                "gene_set": gene_set,
                                "condition": condition,
                            }
                        )

    gene_sets_df = pl.DataFrame(rows)

    # One Fisher test per (gene_set, condition), pooled across studies.
    # Genes outside the universe cannot be classified and are dropped from each set.
    combos = gene_sets_df.select(["gene_set", "condition"]).unique().iter_rows()

    results = []
    for gs, cond in combos:
        universe = universe_for(gs)
        approved_in_universe = approved_for_disease & universe
        n_universe = len(universe)
        n_approved = len(approved_in_universe)

        genes_in_set = (
            set(
                gene_sets_df.filter(
                    (pl.col("gene_set") == gs) & (pl.col("condition") == cond)
                )
                .get_column("gene")
                .to_list()
            )
            & universe
        )

        a = len(genes_in_set & approved_in_universe)
        b = len(genes_in_set) - a
        c = n_approved - a
        d = n_universe - len(genes_in_set) - c

        odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")

        results.append(
            {
                "gene_set": gs,
                "condition": cond,
                "odds_ratio": odds_ratio,
                "p_value": p_value,
                "a": a,
                "b": b,
                "c": c,
                "d": d,
                "n_universe": n_universe,
                "n_approved": n_approved,
            }
        )

    fisher_df = pl.DataFrame(results).sort(["gene_set", "condition"])
    with pl.Config(tbl_cols=-1, tbl_width_chars=-1):
        print(fisher_df)


if __name__ == "__main__":
    main()
