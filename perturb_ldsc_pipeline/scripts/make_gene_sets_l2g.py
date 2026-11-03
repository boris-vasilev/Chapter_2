#!/usr/bin/env python
"""L2G (Locus-to-Gene) gene-set module — OpenTargets parquet input.

For each (study, condition), produces gene-set files:

  * l2g_genes      — OpenTargets-prioritised genes above score threshold
  * regulators_l2g — perturbed cis genes whose perturbation changes
                     expression of any L2G gene
  * targets_l2g    — trans DEGs when L2G genes themselves are perturbed

OpenTargets parquet layout (under --opentargets-dir):
  l2g/*.parquet           studyLocusId, geneId, score
  credible_set/*.parquet  studyLocusId, studyId        (used for the join)

The `--studies` arguments are expected to be OpenTargets studyIds
directly — matching the `study` column in the GWAS_study_info.csv
that step 2 (compile_studies_from_opentargets.py) wrote.
"""

import argparse
import sys
from pathlib import Path

import polars as pl


def load_l2g(opentargets_dir, score_threshold):
    """Return DataFrame with columns [studyLocusId, geneId, score, studyId]."""
    l2g = (
        pl.scan_parquet(str(opentargets_dir / "l2g" / "*.parquet"))
        .filter(pl.col("score") > score_threshold)
        .select(["studyLocusId", "geneId", "score"])
    )
    cs = (
        pl.scan_parquet(str(opentargets_dir / "credible_set" / "*.parquet"))
        .select(["studyLocusId", "studyId"])
    )
    return l2g.join(cs, on="studyLocusId", how="inner").collect()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--de-results",      required=True)
    parser.add_argument("--opentargets-dir", required=True,
                        help="Root OpenTargets dir containing l2g/ credible_set/ ...")
    parser.add_argument("--studies",         nargs="+", required=True,
                        help="OpenTargets studyIds (= 'study' column in GWAS_study_info.csv)")
    parser.add_argument("--conditions",      nargs="+", required=True)
    parser.add_argument("--score-threshold", type=float, default=0.5)
    parser.add_argument("--out-dir",         required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ot_dir = Path(args.opentargets_dir)
    l2g = load_l2g(ot_dir, args.score_threshold)
    print(
        f"Loaded {l2g.height} L2G rows "
        f"(score > {args.score_threshold}, joined to studyId)",
        flush=True,
    )

    available_studies = set(l2g["studyId"].unique().to_list())
    missing = [s for s in args.studies if s not in available_studies]
    if missing:
        print(
            f"WARNING: {len(missing)} study(ies) absent from L2G data — "
            f"gene set files will be empty for: {missing}",
            file=sys.stderr,
        )

    de_dat = pl.read_csv(args.de_results)

    sets = {}
    for study in args.studies:
        study_genes = (
            l2g.filter(pl.col("studyId") == study)["geneId"].unique()
        )
        print(
            f"  {study}: {study_genes.len()} prioritised genes",
            flush=True,
        )

        for condition in args.conditions:
            perturb = de_dat.filter(pl.col("culture_condition") == condition)

            sets[(study, condition, "l2g_genes")] = study_genes

            sets[(study, condition, "regulators_l2g")] = perturb.filter(
                pl.col("trans_ensg").is_in(study_genes),
                pl.col("trans_ensg") != pl.col("cis_ensg"),
            )["cis_ensg"].unique()

            sets[(study, condition, "targets_l2g")] = perturb.filter(
                pl.col("cis_ensg").is_in(study_genes),
                pl.col("trans_ensg") != pl.col("cis_ensg"),
            )["trans_ensg"].unique()

    for (study, condition, annot), genes in sets.items():
        fname = f"{study}_{condition}_{annot}.txt"
        pl.DataFrame({"GENE": genes}).write_csv(
            out_dir / fname, separator="\t", include_header=False
        )


if __name__ == "__main__":
    main()
