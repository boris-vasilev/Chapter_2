#!/usr/bin/env python
"""Coloc gene-set module.

For each (study, condition), produces gene-set files:

  * regulators_H4 — perturbed cis genes that change expression of H4
                    (colocalised) genes
  * targets_H4    — DEGs when H4 genes themselves are perturbed
  * H4_genes      — H4 genes themselves
  * regulators_H3, targets_H3, H3_genes — same for H3 (independent
                    signals); skipped when --h4-only is set

The ``--studies`` arguments must match the values in coloc_dat.tsv's
``imd`` column (typically the full study identifier such as
``UC_DeLange_28067908_1``). Outputs are written as
``{out_dir}/{study}_{condition}_{annot}.txt``, one Ensembl gene ID
per line, no header.
"""

import argparse
import sys
from pathlib import Path

import polars as pl


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--de-results", required=True)
    parser.add_argument("--coloc", required=True)
    parser.add_argument("--studies", nargs="+", required=True)
    parser.add_argument("--conditions", nargs="+", required=True)
    parser.add_argument("--cell-type", default="CD4")
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
    parser.add_argument(
        "--min-abs-logfc",
        type=float,
        default=0.0,
        help="Keep only DE edges with |logFC| >= this when building "
        "regulator/target sets (the anchor/H4 set is unaffected). "
        "Default 0.0 - no filter (current behaviour). Used by the "
        "stringency sensitivity sweep.",
    )

    parser.add_argument(
        "--hubs",
        default=None,
        help="Hub list from extract_de_results.py --hub-list-out. If given, also "
        "write a trait-agnostic 'hubs' set per study and condition (identical across studies); "
        "with several quantiles, one 'hubs_q<quantile>' set per quantile.",
    )

    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    de_dat = pl.read_csv(args.de_results)
    # All columns as strings: the quantile names files exactly as given (0.90, not 0.9).
    hubs = pl.read_csv(args.hubs, infer_schema=False) if args.hubs else None
    if hubs is not None and "quantile" not in hubs.columns:  # older single-quantile list
        hubs = hubs.with_columns(pl.lit("").alias("quantile"))
    hub_qs = hubs["quantile"].unique(maintain_order=True).to_list() if hubs is not None else []

    if args.min_abs_logfc > 0:
        n0 = de_dat.height
        de_dat = de_dat.filter(pl.col("logFC").abs() >= args.min_abs_logfc)
        print(
            f"|logFC| >= {args.min_abs_logfc}: kept {de_dat.height}/{n0} DE edges",
            file=sys.stderr,
        )
    coloc_dat = (
        pl.scan_csv(args.coloc, separator="\t", infer_schema_length=10000)
        .filter(pl.col("cell_type") == args.cell_type)
        # Remove version suffix from gene IDs to match DE results
        .with_columns(pl.col("Gene_id").str.replace(r"\..*", "").alias("Gene_id"))
        .collect()
    )

    h3_all = coloc_dat.filter(pl.col("H3") > args.h3_threshold)
    h4_all = coloc_dat.filter(pl.col("H4") > args.h4_threshold)

    available_imds = set(coloc_dat["imd"].unique().to_list())
    missing = [s for s in args.studies if s not in available_imds]
    if missing:
        print(
            f"WARNING: {len(missing)} study identifier(s) not present in "
            f"coloc_dat 'imd' column — gene set files will be empty for: "
            f"{missing}",
            file=sys.stderr,
        )

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
    for study in args.studies:
        h4 = h4_all.filter(pl.col("imd") == study)
        h3 = h3_all.filter(pl.col("imd") == study)
        h4_genes = h4["Gene_id"]
        h3_genes = h3["Gene_id"]

        for condition in args.conditions:
            perturb = de_dat.filter(pl.col("culture_condition") == condition)

            sets[(study, condition, "regulators_H4")] = convergent_regulators(
                perturb, h4_genes, min_reg_k
            )
            sets[(study, condition, "targets_H4")] = convergent_targets(
                perturb, h4_genes, min_tgt_k
            )
            sets[(study, condition, "H4_genes")] = h4_genes.unique()

            for q in hub_qs:
                name = "hubs" if len(hub_qs) == 1 else f"hubs_q{q}"
                sets[(study, condition, name)] = hubs.filter(
                    (pl.col("culture_condition") == condition) & (pl.col("quantile") == q)
                )["cis_ensg"].unique()

            if not args.h4_only:
                sets[(study, condition, "regulators_H3")] = convergent_regulators(
                    perturb, h3_genes, min_reg_k
                )
                sets[(study, condition, "targets_H3")] = convergent_targets(
                    perturb, h3_genes, min_tgt_k
                )
                sets[(study, condition, "H3_genes")] = h3_genes.unique()

    print(
        f"Filters in use: min anchors/regulator = {min_reg_k}, "
        f"min perturbations/target = {min_tgt_k}",
        file=sys.stderr,
    )
    for (study, condition, annot), genes in sets.items():
        fname = f"{study}_{condition}_{annot}.txt"
        pl.DataFrame({"GENE": genes}).sort("GENE").write_csv(
            out_dir / fname, separator="\t", include_header=False
        )
        print(f"  {fname}\t{len(genes)} genes", file=sys.stderr)


if __name__ == "__main__":
    main()
