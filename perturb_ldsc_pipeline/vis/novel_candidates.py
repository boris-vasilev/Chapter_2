#!/usr/bin/env python
"""Rank novel candidate targets from Perturb-seq-expanded coloc gene sets.

For each disease in DISEASE_MAPPING, takes the union of all candidate genes
that appear in `regulators_H4` or `targets_H4` (across all three conditions),
excludes those already approved/Phase II+ for the focal disease, and ranks the
rest. Each candidate is annotated with:

  - the coloc (H4) gene(s) it links to in the Perturb-seq DE results, with the
    strongest connection (largest |logFC| × -log10(padj)) used for ranking;
  - existing drug evidence in Open Targets for OTHER diseases (drug name,
    best clinical stage, disease the drug is approved/in trials for);
  - a tier label:
      Tier 1 (Repositioning)   — has an approved / Phase IV / Preapproval drug
                                 for at least one non-focal, non-oncology disease.
      Tier 2 (Tractable novel) — has some drug development activity (any phase,
                                 incl. oncology), but no approval for a
                                 non-oncology indication. Oncology approvals are
                                 excluded from Tier 1 because cancer T-cell
                                 engagers / BiTEs (most of the CD3, CD20, BCMA
                                 etc. "approvals") use the protein with a
                                 fundamentally different mechanism from what
                                 would be needed for an autoimmune indication.
      Tier 3 (Fully novel)     — no drug evidence anywhere in OT.

Outputs:
  <out-dir>/candidates_<DISEASE>.csv  — one ranked CSV per disease.
  <out-dir>/summary.csv               — per-disease totals and tier counts.
"""

import argparse
import warnings
from math import sqrt
from pathlib import Path

import polars as pl

# Stages used to define "already validated for the focal disease".
APPROVED_FOCAL_STAGES = [
    "APPROVAL",
    "PREAPPROVAL",
    "PHASE_4",
    "PHASE_3",
    "PHASE_2_3",
    "PHASE_2",
]
# Strongest stages — used to classify Tier 1 repositioning candidates.
APPROVED_TIER1_STAGES = {"APPROVAL", "PREAPPROVAL", "PHASE_4"}

# Sorted from strongest evidence to weakest. Used to pick max stage per gene.
STAGE_PRIORITY = [
    "APPROVAL",
    "PHASE_4",
    "PREAPPROVAL",
    "PHASE_3",
    "PHASE_2_3",
    "PHASE_2",
    "PHASE_1_2",
    "PHASE_1",
    "EARLY_PHASE_1",
    "IND",
    "PRECLINICAL",
    "WITHDRAWAL",
    "UNKNOWN",
]

CONDITIONS = ["Rest", "Stim8hr", "Stim48hr"]
GENE_SET_TYPES = ["regulators_H4", "targets_H4"]

# Diseases whose subtree is filtered out of "drugs_other_diseases" and the
# Tier 1 stage check. MONDO_0045024
ONCOLOGY_ROOT_IDS = ["MONDO_0045024"]

DISEASE_MAPPING: dict[str, tuple[str, str]] = {
    # IMDs
    "AST": ("Asthma", "MONDO_0004979"),
    "ASTAO": ("Asthma (adult-onset)", "MONDO_0004979"),
    "ASTCO": ("Asthma (childhood-onset)", "MONDO_0004979"),
    "ATD": ("Atopic dermatitis", "EFO_0000274"),
    "AITD": ("Autoimmune thyroid disease", "EFO_1001054"),
    "CD": ("Crohn's disease", "EFO_0000384"),
    "IBD": ("Inflammatory bowel disease", "EFO_0003767"),
    "UC": ("Ulcerative colitis", "EFO_0000729"),
    "RA": ("Rheumatoid arthritis", "EFO_0000685"),
    "SLE": ("Systemic lupus erythematosus", "MONDO_0007915"),
    "T1D": ("Type 1 diabetes", "MONDO_0005147"),
    # Non-IMDs
    "IS": ("Ischemic stroke", "EFO_0009601"),
    "MYOI": ("Myocardial infarction", "EFO_0000612"),
    "OA": ("Osteoarthritis", "EFO_0002506"),
    "T2D": ("Type 2 diabetes", "MONDO_0005148"),
    "ADFH": ("Alzheimer's disease (Jansen)", "EFO_0000249"),
}


def stage_rank(stage: str) -> int:
    try:
        return STAGE_PRIORITY.index(stage)
    except ValueError:
        return len(STAGE_PRIORITY)


def best_stage(stages: list[str]) -> str:
    return min(stages, key=stage_rank) if stages else ""


def classify_tier(n_drugs_any: int, max_stage_anywhere: str) -> str:
    if n_drugs_any == 0:
        return "Tier3_FullyNovel"
    if max_stage_anywhere in APPROVED_TIER1_STAGES:
        return "Tier1_Repositioning"
    return "Tier2_TractableNovel"


def aggregate_zscores(abs_zs: list[float]) -> dict[str, float]:
    """Aggregate per-connection |z| into the three ranking metrics.

    - max_abs_z       : strongest single coloc-gene link (conservative).
    - sum_abs_z       : total |z| evidence — favours multi-coloc hubs but scales
                        linearly with N, so weak hubs can inflate.
    - stouffer_abs_z  : sum_abs_z / sqrt(N). Under independence, Σz / √N has
                        unit variance; this is the principled hub-aware score.
    """
    if not abs_zs:
        return {"max_abs_z": 0.0, "sum_abs_z": 0.0, "stouffer_abs_z": 0.0}
    s = sum(abs_zs)
    return {
        "max_abs_z": max(abs_zs),
        "sum_abs_z": s,
        "stouffer_abs_z": s / sqrt(len(abs_zs)),
    }


def resolve_disease_id(query: str, ontology: pl.DataFrame) -> str | None:
    if ontology.filter(pl.col("id") == query).height > 0:
        return query
    for candidate in (query, query.replace("_", ":")):
        m = ontology.filter(pl.col("dbXRefs").list.contains(candidate))
        if m.height > 0:
            return m.get_column("id").item(0)
    return None


def expand_subtree(disease_id: str, ontology: pl.DataFrame) -> set[str]:
    m = ontology.filter(pl.col("id") == disease_id)
    if m.height == 0:
        return {disease_id}
    desc = m.get_column("descendants").to_list()[0]
    return {disease_id} | set(desc or [])


def load_gene_set(
    gene_sets_dir: Path, study: str, gene_set: str, condition: str
) -> set[str]:
    genes: set[str] = set()
    for f in gene_sets_dir.glob(f"{study}_*_{condition}_{gene_set}.txt"):
        with open(f) as fh:
            for line in fh:
                g = line.strip()
                if g:
                    genes.add(g)
    return genes


def get_symbol(ensg: str, de: pl.DataFrame) -> str:
    for ensg_col, sym_col in [("cis_ensg", "cis_gene"), ("trans_ensg", "trans_gene")]:
        s = de.filter(pl.col(ensg_col) == ensg).get_column(sym_col)
        if s.len() > 0:
            return s.item(0)
    return ""


def get_connections(
    ensg: str,
    gene_set_type: str,
    h4_genes: set[str],
    de_condition: pl.DataFrame,
) -> list[dict]:
    """Return DE-evidence linking this candidate to coloc (H4) genes."""
    if gene_set_type == "targets_H4":
        sub = de_condition.filter(
            (pl.col("trans_ensg") == ensg) & (pl.col("cis_ensg").is_in(h4_genes))
        )
        coloc_col = "cis_gene"
    elif gene_set_type == "regulators_H4":
        sub = de_condition.filter(
            (pl.col("cis_ensg") == ensg) & (pl.col("trans_ensg").is_in(h4_genes))
        )
        coloc_col = "trans_gene"
    else:
        return []
    return [
        {
            "coloc_gene": r[coloc_col],
            "logFC": r["logFC"],
            "padj": r["padj"],
            "zscore": r["zscore"],
            "gene_set_type": gene_set_type,
            "condition": r["culture_condition"],
        }
        for r in sub.iter_rows(named=True)
    ]


def process_disease(
    code: str,
    disease_id: str,
    ontology: pl.DataFrame,
    drug_targets: pl.DataFrame,
    de: pl.DataFrame,
    gene_sets_dir: Path,
    disease_name_lookup: dict[str, str],
    oncology_subtree: set[str],
) -> pl.DataFrame | None:
    subtree = expand_subtree(disease_id, ontology)

    # Genes already validated for the focal disease (Phase II+) — excluded as candidates.
    focal_drugged = set(
        drug_targets.filter(pl.col("diseaseId").is_in(subtree))
        .filter(pl.col("clinicalStage").is_in(APPROVED_FOCAL_STAGES))
        .get_column("targetId")
        .unique()
        .to_list()
    )

    de_per_cond = {
        cond: de.filter(pl.col("culture_condition") == cond) for cond in CONDITIONS
    }
    h4_per_cond = {
        cond: load_gene_set(gene_sets_dir, code, "H4_genes", cond)
        for cond in CONDITIONS
    }

    # Collect candidates across all (gene-set type, condition) combinations.
    candidates: dict[str, dict] = {}
    for gs_type in GENE_SET_TYPES:
        for cond in CONDITIONS:
            members = load_gene_set(gene_sets_dir, code, gs_type, cond)
            for ensg in members:
                if ensg in focal_drugged:
                    continue
                if ensg not in candidates:
                    candidates[ensg] = {"appearances": set(), "connections": []}
                candidates[ensg]["appearances"].add((gs_type, cond))
                candidates[ensg]["connections"].extend(
                    get_connections(ensg, gs_type, h4_per_cond[cond], de_per_cond[cond])
                )

    if not candidates:
        return None

    # Drug evidence rows for all candidates in one pass.
    cand_ensgs = set(candidates.keys())
    drug_rows = drug_targets.filter(pl.col("targetId").is_in(cand_ensgs))

    rows: list[dict] = []
    for ensg, info in candidates.items():
        gene_drugs = drug_rows.filter(pl.col("targetId") == ensg)
        n_drugs_any = gene_drugs.get_column("drugFromSource").drop_nulls().n_unique()
        all_stages = (
            gene_drugs.get_column("clinicalStage").drop_nulls().unique().to_list()
        )
        ms_anywhere = best_stage(all_stages)

        # Drugs for non-focal, non-oncology diseases — collapse to (drug -> best
        # stage, primary disease). Oncology indications are excluded because
        # cancer T-cell engagers/BiTEs that hit immune-axis proteins like CD3,
        # CD20, BCMA etc. don't represent valid repositioning precedent for
        # autoimmune disease (the modality redirects T cells; the autoimmune use
        # would suppress them).
        other = gene_drugs.filter(
            ~pl.col("diseaseId").is_in(subtree)
            & ~pl.col("diseaseId").is_in(oncology_subtree)
        ).select(["drugFromSource", "diseaseId", "clinicalStage"])
        drug_best: dict[str, tuple[str, str]] = {}
        for r in other.drop_nulls().unique().iter_rows(named=True):
            d, dis, s = r["drugFromSource"], r["diseaseId"], r["clinicalStage"]
            if d not in drug_best or stage_rank(s) < stage_rank(drug_best[d][0]):
                drug_best[d] = (s, dis)
        ms_other = best_stage([s for s, _ in drug_best.values()])
        drugs_other_str = "; ".join(
            f"{d} ({s} for {disease_name_lookup.get(dis, dis)})"
            for d, (s, dis) in sorted(
                drug_best.items(), key=lambda kv: stage_rank(kv[1][0])
            )
        )

        # Z-score-based ranking over all (coloc-gene, condition) appearances.
        conns = info["connections"]
        abs_zs = [abs(c["zscore"]) for c in conns if c.get("zscore") is not None]
        scores = aggregate_zscores(abs_zs)
        if conns:
            # Define the "top" connection as the one with the largest |z|.
            best = max(
                conns,
                key=lambda c: abs(c["zscore"]) if c.get("zscore") is not None else 0.0,
            )
            top_coloc = best["coloc_gene"]
            top_logfc = best["logFC"]
            top_padj = best["padj"]
            top_z = best.get("zscore")
            unique_coloc = sorted({c["coloc_gene"] for c in conns if c["coloc_gene"]})
            coloc_genes_str = "; ".join(unique_coloc)
        else:
            top_coloc, top_logfc, top_padj, top_z = "", None, None, None
            coloc_genes_str = ""

        appearances = sorted(info["appearances"])
        appearance_str = "; ".join(f"{gs}/{c}" for gs, c in appearances)

        rows.append(
            {
                "disease": code,
                "gene_ensg": ensg,
                "gene_symbol": get_symbol(ensg, de),
                "appearances": appearance_str,
                "n_coloc_links": len(
                    {c["coloc_gene"] for c in conns if c["coloc_gene"]}
                ),
                "n_connections": len(conns),
                "coloc_genes": coloc_genes_str,
                "top_coloc_gene": top_coloc,
                "top_logFC": top_logfc,
                "top_padj": top_padj,
                "top_zscore": top_z,
                # Three z-score-based ranking metrics; sort defaults to Stouffer.
                "max_abs_z": scores["max_abs_z"],
                "sum_abs_z": scores["sum_abs_z"],
                "stouffer_abs_z": scores["stouffer_abs_z"],
                "n_drugs_any_disease": n_drugs_any,
                "max_stage_anywhere": ms_anywhere,
                "max_stage_other_disease": ms_other,
                "drugs_other_diseases": drugs_other_str,
                # Tier 1 requires a non-oncology approval, so classification
                # uses ms_other (the filtered max stage) rather than ms_anywhere.
                "tier": classify_tier(n_drugs_any, ms_other),
            }
        )

    return pl.DataFrame(rows).sort(["tier", "max_abs_z"], descending=[False, True])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-sets-dir", required=True)
    parser.add_argument("--drug-targets", required=True)
    parser.add_argument("--disease-ontology", required=True)
    parser.add_argument("--de-results", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    gene_sets_dir = Path(args.gene_sets_dir)

    print("Loading ontology, drug-targets, DE results...")
    ontology = (
        pl.scan_parquet(args.disease_ontology)
        .select(["id", "name", "descendants", "dbXRefs"])
        .collect()
    )
    drug_targets = (
        pl.scan_parquet(args.drug_targets)
        .select(["targetId", "diseaseId", "clinicalStage", "drugFromSource"])
        .collect()
    )
    de = pl.read_csv(args.de_results)
    disease_name_lookup = dict(
        zip(
            ontology.get_column("id").to_list(),
            ontology.get_column("name").to_list(),
        )
    )

    # Build the oncology subtree once. Diseases inside it are filtered out of
    # the "drugs_other_diseases" evidence and the Tier 1 stage check.
    oncology_subtree: set[str] = set()
    for onco_id in ONCOLOGY_ROOT_IDS:
        resolved_onco = resolve_disease_id(onco_id, ontology)
        if resolved_onco is None:
            warnings.warn(f"Oncology root '{onco_id}' not found in ontology")
            continue
        oncology_subtree |= expand_subtree(resolved_onco, ontology)
    print(f"  |oncology subtree|={len(oncology_subtree)}")

    summary_rows: list[dict] = []
    for code, (display, raw_id) in DISEASE_MAPPING.items():
        resolved = resolve_disease_id(raw_id, ontology)
        if resolved is None:
            warnings.warn(f"[{code}] disease ID '{raw_id}' not found; skipping")
            continue

        df = process_disease(
            code,
            resolved,
            ontology,
            drug_targets,
            de,
            gene_sets_dir,
            disease_name_lookup,
            oncology_subtree,
        )
        if df is None or df.height == 0:
            warnings.warn(f"[{code}] no candidates; skipping")
            continue

        out_path = out_dir / f"candidates_{code}.csv"
        df.write_csv(out_path)

        tier_counts = {
            row["tier"]: row["count"]
            for row in df.group_by("tier")
            .agg(pl.len().alias("count"))
            .iter_rows(named=True)
        }
        summary_rows.append(
            {
                "disease": code,
                "display": display,
                "disease_id": resolved,
                "total": df.height,
                "tier1_repositioning": tier_counts.get("Tier1_Repositioning", 0),
                "tier2_tractable": tier_counts.get("Tier2_TractableNovel", 0),
                "tier3_novel": tier_counts.get("Tier3_FullyNovel", 0),
            }
        )

        print(
            f"[{code}] {display}: {df.height} candidates "
            f"(T1={tier_counts.get('Tier1_Repositioning', 0)}, "
            f"T2={tier_counts.get('Tier2_TractableNovel', 0)}, "
            f"T3={tier_counts.get('Tier3_FullyNovel', 0)})  -> {out_path}"
        )

    summary_df = pl.DataFrame(summary_rows)
    summary_path = out_dir / "summary.csv"
    summary_df.write_csv(summary_path)
    print(f"\nWrote summary -> {summary_path}")


if __name__ == "__main__":
    main()
