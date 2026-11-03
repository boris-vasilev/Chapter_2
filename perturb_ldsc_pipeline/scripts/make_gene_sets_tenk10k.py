#!/usr/bin/env python
"""Coloc gene-set module for the tenk10k coloc results.

Mirror of ``make_gene_sets_coloc.py`` but reading the tenk10k coloc
directory tree instead of a single ``coloc_dat.tsv``. The tree is laid
out as::

    {coloc_root}/{disease}/{cell_type}/chr{N}.csv

with one CSV per (disease, cell_type, chromosome). Each CSV has columns
``gene`` (Ensembl ID, no version suffix), ``PP.H3.abf``, ``PP.H4.abf``,
``celltype``, ``chrom``, ``top_snp`` (and the other PP.* columns).

For each (disease, condition) it produces the same gene-set files as the
original module:

  * regulators_H4 — perturbed cis genes that change expression of H4
                    (colocalised) genes
  * targets_H4    — DEGs when H4 genes themselves are perturbed
  * H4_genes      — H4 genes themselves
  * regulators_H3, targets_H3, H3_genes — same for H3 (independent
                    signals); skipped when --h4-only is set

H4/H3 genes are taken as the **union across CD4 subtypes**: a gene is an
H4 gene for a disease if its PP.H4 exceeds the threshold in *any* of the
selected cell types (default the CD4 subtypes), matching the CD4-based
Perturb-seq DE results used to derive regulators/targets.

The disease (study) code in every output filename is the disease
*directory name* (e.g. ``ra_GCST90132223``). Outputs are written as
``{out_dir}/{disease}_{condition}_{annot}.txt``, one Ensembl gene ID per
line, no header — the same convention ``approved_drug_enrichment_all.py``
expects.

A ``disease_groups.csv`` is also written to ``out_dir`` with the inferred
IMD / other classification, display names and OpenTargets disease IDs, so
the enrichment heatmap separation works (pass it via
``--disease-groups``).

NOTE: ``approved_drug_enrichment_all.py`` resolves gene sets with
``glob('{study}_*')``, so a study code must not be a prefix of another
(e.g. ``ibd`` vs ``ibd_liu2023``). Keep DiseaseTraits and NewDiseases in
separate output directories / runs to avoid collisions.
"""

import argparse
import glob
import sys
from pathlib import Path

import polars as pl

# Default CD4 subtypes to union over — chosen to match the CD4-based
# Perturb-seq DE results used for regulators/targets.
CD4_SUBTYPES = [
    "CD4_CTL",
    "CD4_Naive",
    "CD4_Proliferating",
    "CD4_TCM",
    "CD4_TEM",
]

# Inferred metadata for the tenk10k DiseaseTraits + NewDiseases diseases.
# Keyed by directory name -> (display_name, short_label, group, ot_id).
# group is "IMD" or "other". OT IDs are best-guess defaults; the enrichment
# script falls back to dbXRefs resolution, so MONDO_/EFO_ forms both work.
DISEASE_META: dict[str, tuple[str, str, str, str]] = {
    # --- DiseaseTraits (mixed IMD / other) ---
    "NHL_GCST90011819": ("Non-Hodgkin lymphoma", "NHL", "other", "EFO_0005952"),
    "alzheimer_GCST90027158": ("Alzheimer's disease", "Alzheimer's", "other", "MONDO_0004975"),
    "breastca_GCST004988": ("Breast cancer", "Breast ca.", "other", "MONDO_0007254"),
    "covid_GCST011071": ("COVID-19", "COVID-19", "other", "MONDO_0100096"),
    "ibd_liu2023": ("Inflammatory bowel disease", "IBD", "IMD", "EFO_0003767"),
    "kiryluk_IgAN": ("IgA nephropathy", "IgA neph.", "IMD", "EFO_0004194"),
    "lungca_GCST004748": ("Lung cancer", "Lung ca.", "other", "MONDO_0008903"),
    "lymphocytic_leukemia_GCST90011814": (
        "Chronic lymphocytic leukemia", "CLL", "other", "EFO_0000095",
    ),
    "lymphoma_GCST90018878": ("Lymphoma", "Lymphoma", "other", "EFO_0000574"),
    "myeloproliferative_GCST90000032": (
        "Myeloproliferative neoplasm", "Myeloprolif.", "other", "MONDO_0006574",
    ),
    "nephrotic_GCST90258619": ("Nephrotic syndrome", "Nephrotic", "IMD", "EFO_0004255"),
    "parkinson_GCST009325": ("Parkinson's disease", "Parkinson's", "other", "MONDO_0005180"),
    "prostateca_GCST90274713": ("Prostate cancer", "Prostate ca.", "other", "MONDO_0008315"),
    "ra_GCST90132223": ("Rheumatoid arthritis", "Rheum. arthritis", "IMD", "EFO_0000685"),
    "sle_GCST003156": ("Systemic lupus erythematosus", "SLE", "IMD", "MONDO_0007915"),
    # --- NewDiseases (all IMD) ---
    "asthma": ("Asthma", "Asthma", "IMD", "MONDO_0004979"),
    "cd": ("Crohn's disease", "Crohn's", "IMD", "EFO_0000384"),
    "ibd": ("Inflammatory bowel disease", "IBD", "IMD", "EFO_0003767"),
    "ms": ("Multiple sclerosis", "MS", "IMD", "MONDO_0005301"),
    "t1dm": ("Type 1 diabetes", "T1D", "IMD", "MONDO_0005147"),
    "uc": ("Ulcerative colitis", "Ulcer. colitis", "IMD", "EFO_0000729"),
}


def load_disease_coloc(
    coloc_root: Path, disease: str, cell_types: list[str]
) -> pl.DataFrame:
    """Read + concatenate all chr*.csv for a disease across cell types.

    Returns a DataFrame with columns Gene_id, H3, H4, cell_type (gene IDs
    stripped of any version suffix to match the DE results).
    """
    frames = []
    for ct in cell_types:
        pattern = str(coloc_root / disease / ct / "chr*.csv")
        if not glob.glob(pattern):
            continue
        frames.append(
            pl.scan_csv(pattern, infer_schema_length=10000).select(
                pl.col("gene").str.replace(r"\..*", "").alias("Gene_id"),
                pl.col("PP.H3.abf").alias("H3"),
                pl.col("PP.H4.abf").alias("H4"),
                pl.col("celltype").alias("cell_type"),
            )
        )
    if not frames:
        return pl.DataFrame(
            schema={"Gene_id": pl.Utf8, "H3": pl.Float64, "H4": pl.Float64, "cell_type": pl.Utf8}
        )
    return pl.concat(frames).collect()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--de-results", required=True)
    parser.add_argument(
        "--coloc-root",
        required=True,
        help="Directory holding per-disease subdirectories, e.g. "
        "data/tenk10k_coloc/DiseaseTraits or data/tenk10k_coloc/NewDiseases.",
    )
    parser.add_argument(
        "--diseases",
        nargs="+",
        default=None,
        help="Disease directory names to process. Default: every "
        "subdirectory of --coloc-root.",
    )
    parser.add_argument(
        "--conditions", nargs="+", default=["Rest", "Stim8hr", "Stim48hr"]
    )
    parser.add_argument(
        "--cell-types",
        nargs="+",
        default=CD4_SUBTYPES,
        help="tenk10k cell types to union H4/H3 genes over. "
        f"Default: {CD4_SUBTYPES}.",
    )
    parser.add_argument("--h3-threshold", type=float, default=0.5)
    parser.add_argument("--h4-threshold", type=float, default=0.5)
    parser.add_argument("--h4-only", action="store_true", default=True)
    parser.add_argument(
        "--min-anchors-per-regulator",
        type=int,
        default=1,
        help="Keep a regulator (perturbed cis gene) only if its perturbation "
        "changes expression of at least this many distinct coloc "
        "anchors. Default 1 - no filter.",
    )
    parser.add_argument(
        "--min-perturbations-per-target",
        type=int,
        default=1,
        help="Keep a target (trans gene / DEG) only if it is changed by "
        "perturbations of at least this many distinct coloc anchors. "
        "Default 1 - no filter.",
    )
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    coloc_root = Path(args.coloc_root)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.diseases:
        diseases = args.diseases
    else:
        diseases = sorted(p.name for p in coloc_root.iterdir() if p.is_dir())

    de_dat = pl.read_csv(args.de_results)

    def convergent_regulators(perturb, anchors, min_k):
        """Cis genes whose perturbation changes ≥ min_k distinct anchors."""
        return (
            perturb.filter(
                pl.col("trans_ensg").is_in(anchors),
                pl.col("trans_ensg") != pl.col("cis_ensg"),
            )
            .group_by("cis_ensg")
            .agg(pl.col("trans_ensg").n_unique().alias("n_anchors_hit"))
            .filter(pl.col("n_anchors_hit") >= min_k)["cis_ensg"]
            .unique()
        )

    def convergent_targets(perturb, anchors, min_k):
        """Trans genes hit by ≥ min_k distinct anchor perturbations."""
        return (
            perturb.filter(
                pl.col("cis_ensg").is_in(anchors),
                pl.col("trans_ensg") != pl.col("cis_ensg"),
            )
            .group_by("trans_ensg")
            .agg(pl.col("cis_ensg").n_unique().alias("n_perturbations_hitting"))
            .filter(pl.col("n_perturbations_hitting") >= min_k)["trans_ensg"]
            .unique()
        )

    min_reg_k = args.min_anchors_per_regulator
    min_tgt_k = args.min_perturbations_per_target

    sets = {}
    groups_rows = []
    for disease in diseases:
        coloc_dat = load_disease_coloc(coloc_root, disease, args.cell_types)
        if coloc_dat.height == 0:
            print(
                f"WARNING: no CSVs for '{disease}' under {coloc_root} for cell "
                f"types {args.cell_types}; gene sets will be empty.",
                file=sys.stderr,
            )

        # Union across the selected cell types: a gene qualifies if it
        # passes threshold in any one of them.
        h4_genes = coloc_dat.filter(pl.col("H4") > args.h4_threshold)["Gene_id"].unique()
        h3_genes = coloc_dat.filter(pl.col("H3") > args.h3_threshold)["Gene_id"].unique()

        for condition in args.conditions:
            perturb = de_dat.filter(pl.col("culture_condition") == condition)

            sets[(disease, condition, "regulators_H4")] = convergent_regulators(
                perturb, h4_genes, min_reg_k
            )
            sets[(disease, condition, "targets_H4")] = convergent_targets(
                perturb, h4_genes, min_tgt_k
            )
            sets[(disease, condition, "H4_genes")] = h4_genes

            if not args.h4_only:
                sets[(disease, condition, "regulators_H3")] = convergent_regulators(
                    perturb, h3_genes, min_reg_k
                )
                sets[(disease, condition, "targets_H3")] = convergent_targets(
                    perturb, h3_genes, min_tgt_k
                )
                sets[(disease, condition, "H3_genes")] = h3_genes

        display, short, group, ot_id = DISEASE_META.get(
            disease, (disease, disease, "other", "")
        )
        if disease not in DISEASE_META:
            print(
                f"WARNING: '{disease}' not in DISEASE_META; defaulting to "
                f"group='other', no OT ID. Add it to make_gene_sets_tenk10k.py.",
                file=sys.stderr,
            )
        groups_rows.append(
            {
                "study": disease,
                "display_name": display,
                "short_label": short,
                "group": group,
                "ot_disease_id": ot_id,
            }
        )

    print(
        f"Filters in use: min anchors/regulator = {min_reg_k}, "
        f"min perturbations/target = {min_tgt_k}; "
        f"cell types unioned = {args.cell_types}",
        file=sys.stderr,
    )
    for (disease, condition, annot), genes in sets.items():
        fname = f"{disease}_{condition}_{annot}.txt"
        pl.DataFrame({"GENE": genes}).write_csv(
            out_dir / fname, separator="\t", include_header=False
        )
        print(f"  {fname}\t{len(genes)} genes", file=sys.stderr)

    groups_csv = out_dir / "disease_groups.csv"
    pl.DataFrame(groups_rows).write_csv(groups_csv)
    print(f"Wrote {groups_csv} ({len(groups_rows)} diseases)", file=sys.stderr)


if __name__ == "__main__":
    main()
