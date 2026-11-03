#!/usr/bin/env python
"""Permutation test for the drug-target enrichment of coloc-Perturb-seq gene sets.

The original analysis (`approved_drug_enrichment_all.py`) shows the expanded
gene sets (regulators / targets of coloc genes) are enriched for approved drug
targets across IMDs. The concern this script addresses: maybe the enrichment is
not specific to disease-relevant coloc biology, but instead reflects that
expanding *any* size-matched gene set through the CD4+ T cell screen captures
the same CD4+ T cell signature — which itself overlaps with IMD drug targets.

What the script does, per disease and per condition:

  1. Load the real coloc gene set (H4_genes) and the real expansions
     (regulators_H4, targets_H4) — these define the *observed* test.
  2. Build a "pseudo-coloc" set by sampling matched genes (see v2 note
     below) from the screen-measured pool minus all coloc genes.
  3. Expand the pseudo-coloc set into "pseudo-regulators" and "pseudo-targets"
     using the same Perturb-seq DE relationships as the real expansion
     (precomputed cis→trans and trans→cis indices from --de-results).
  4. Re-run the Fisher exact test against the same approved-drug-target set.
  5. Repeat N times to get an empirical null distribution of odds ratios.
  6. Emit empirical p-values = fraction of null ORs ≥ real OR.

v2 (perturb_heritability_null_v2): matching is always on and comes from
scripts/matching.py, shared with the S-LDSC null. Per condition, each real
anchor is matched on its exact (out-degree, in-degree, expression) stratum,
with degrees from the same |logFC| >= --min-abs-logfc edges that expand the
anchors; all coloc genes are excluded from the pool; anchors not measured in
the screen are matched to unmeasured genes of the same coding status. See
NOTES_null_v2.md for the list of changes from
perturb_ldsc_pipeline/vis/approved_drug_enrichment_permutation.py.

Outputs:
  <out>/permutation_results.csv  one row per (disease, gene_set, condition):
      real_or, real_p, real_a, null_or_mean, null_or_p50, null_or_p95,
      null_or_max, empirical_p, n_perm, n_fallback_total, n_unmeasured
  <out>/permutation_nulls.parquet (optional, with --save-null-draws)  raw null
      OR for every (disease, gene_set, condition, draw).
"""

import argparse
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import fisher_exact

# matching.py lives with the pipeline that uses it.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "nf_ldsc_permutation" / "bin"))
from matching import build_matcher, load_coloc_genes as load_all_coloc, perm_rng  # noqa: E402

APPROVED_STAGES = [
    "APPROVAL",
    "PREAPPROVAL",
    "PHASE_4",
    "PHASE_3",
    "PHASE_2_3",
    "PHASE_2",
]

CONDITIONS = ["Rest", "Stim8hr", "Stim48hr"]

# Same mapping as approved_drug_enrichment_all.py — kept in sync.
DISEASE_MAPPING: dict[str, tuple[str, str]] = {
    "AST": ("Asthma", "MONDO_0004979"),
    "ASTAO": ("Asthma (adult-onset)", "MONDO_0004979"),
    "ASTCO": ("Asthma (childhood-onset)", "MONDO_0004979"),
    "ATD": ("Atopic dermatitis", "EFO_0000274"),
    "AITD": ("Autoimmune thyroid disease", "MONDO_0000569"),
    "CD": ("Crohn's disease", "EFO_0000384"),
    "IBD": ("Inflammatory bowel disease", "EFO_0003767"),
    "UC": ("Ulcerative colitis", "EFO_0000729"),
    "RA": ("Rheumatoid arthritis", "EFO_0000685"),
    "SLE": ("Systemic lupus erythematosus", "MONDO_0007915"),
    "T1D": ("Type 1 diabetes", "MONDO_0005147"),
    "IS": ("Ischemic stroke", "HP_0002140"),
    "MYOI": ("Myocardial infarction", "EFO_0000612"),
    "OA": ("Osteoarthritis", "MONDO_0005178"),
    "T2D": ("Type 2 diabetes", "MONDO_0005148"),
    "ADFH": ("Alzheimer's disease", "MONDO_0004975"),
    "SCZ": ("Schizophrenia", "MONDO_0005090"),
}


def resolve_disease_id(query_id: str, ontology: pl.DataFrame) -> str | None:
    if ontology.filter(pl.col("id") == query_id).height > 0:
        return query_id
    for cand in (query_id, query_id.replace("_", ":")):
        m = ontology.filter(pl.col("dbXRefs").list.contains(cand))
        if m.height > 0:
            return m.get_column("id").item(0)
    return None


def expand_subtree(disease_id: str, ontology: pl.DataFrame) -> set[str]:
    m = ontology.filter(pl.col("id") == disease_id)
    if m.height == 0:
        return {disease_id}
    desc = m.get_column("descendants").to_list()[0]
    return {disease_id} | set(desc or [])


def approved_targets(drug_targets: pl.DataFrame, subtree: set[str]) -> set[str]:
    return set(
        drug_targets.filter(pl.col("diseaseId").is_in(subtree))
        .filter(pl.col("clinicalStage").is_in(APPROVED_STAGES))
        .get_column("targetId")
        .unique()
        .to_list()
    )


def load_coloc_genes(gene_sets_dir: Path, study: str, condition: str) -> set[str]:
    """Real H4 coloc genes for one (study, condition)."""
    genes: set[str] = set()
    for f in gene_sets_dir.glob(f"{study}_*_{condition}_H4_genes.txt"):
        with open(f) as fh:
            for line in fh:
                g = line.strip()
                if g:
                    genes.add(g)
    return genes


def expand_pseudo(
    pseudo_coloc: set[str],
    cis_to_trans: dict[str, set[str]],
    trans_to_cis: dict[str, set[str]],
) -> tuple[set[str], set[str]]:
    """For one condition's indices, return (pseudo_regulators, pseudo_targets)."""
    pseudo_targets: set[str] = set()
    for g in pseudo_coloc:
        ts = cis_to_trans.get(g)
        if ts:
            pseudo_targets |= ts
    pseudo_regulators: set[str] = set()
    for g in pseudo_coloc:
        cs = trans_to_cis.get(g)
        if cs:
            pseudo_regulators |= cs
    return pseudo_regulators, pseudo_targets


def fisher_or(
    gene_set: set[str], approved_in_universe: set[str], universe: set[str]
) -> tuple[float, float, int, int, int, int]:
    gs = gene_set & universe
    n_universe = len(universe)
    n_approved = len(approved_in_universe)
    a = len(gs & approved_in_universe)
    b = len(gs) - a
    c = n_approved - a
    d = n_universe - len(gs) - c
    odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
    return odds_ratio, p_value, a, b, c, d


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-sets-dir", required=True)
    parser.add_argument("--drug-targets", required=True)
    parser.add_argument("--disease-ontology", required=True)
    parser.add_argument("--de-results", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument(
        "--n-perm", type=int, default=200, help="Permutation iterations per test."
    )
    parser.add_argument("--seed", type=int, default=42)
    # v2 matching inputs (see scripts/matching.py).
    parser.add_argument("--coloc", required=True,
                        help="coloc_dat.tsv; all coloc genes are excluded from the pool.")
    parser.add_argument("--cell-type", default="CD4")
    parser.add_argument("--h4-threshold", type=float, default=0.8)
    parser.add_argument("--expression", required=True,
                        help="screen_expression.tsv from build_screen_expression.py")
    parser.add_argument("--gene-coord", required=True,
                        help="Gene-coordinate file defining the unmeasured-gene pool.")
    parser.add_argument("--protein-coding-coord", required=True)
    parser.add_argument(
        "--min-abs-logfc", type=float, default=0.5,
        help="Edge filter for BOTH the anchor expansion and the degree matching; "
             "must match the variant in --gene-sets-dir (0.5 = coloc_thr0.5).")
    parser.add_argument("--n-degree-bins", type=int, default=4)
    parser.add_argument("--n-expr-bins", type=int, default=5)
    parser.add_argument(
        "--save-null-draws",
        action="store_true",
        help="Write the raw null OR for every permutation as a parquet file.",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    gene_sets_dir = Path(args.gene_sets_dir)

    print("Loading ontology, drug-targets, DE results...")
    ontology = (
        pl.scan_parquet(args.disease_ontology)
        .select(["id", "descendants", "dbXRefs"])
        .collect()
    )
    drug_targets = (
        pl.scan_parquet(args.drug_targets)
        .select(["targetId", "diseaseId", "clinicalStage"])
        .collect()
    )

    de = pl.read_csv(args.de_results,
                     columns=["cis_ensg", "trans_ensg", "culture_condition", "logFC"])

    # Fisher universes are unchanged: every gene in the DE table.
    cis_universe = {g for g in de.get_column("cis_ensg").unique().to_list() if g}
    trans_universe = {g for g in de.get_column("trans_ensg").unique().to_list() if g}
    print(f"  |cis|={len(cis_universe)}  |trans|={len(trans_universe)}")

    # Per condition: the |logFC|-filtered edges used for expansion AND for
    # in/out-degree, and a matcher over the screen-measured pool minus all
    # coloc genes.
    all_coloc = load_all_coloc(args.coloc, args.cell_type, args.h4_threshold)
    matchers, cis_to_trans, trans_to_cis = {}, {}, {}
    for cond in CONDITIONS:
        matcher, edges = build_matcher(
            de, cond, args.min_abs_logfc, args.expression, all_coloc,
            args.gene_coord, args.protein_coding_coord,
            args.n_degree_bins, args.n_expr_bins)
        matchers[cond] = matcher
        c2t, t2c = defaultdict(set), defaultdict(set)
        for cis, trans in edges.iter_rows():
            c2t[cis].add(trans)
            t2c[trans].add(cis)
        cis_to_trans[cond], trans_to_cis[cond] = c2t, t2c
        print(f"  [{cond}] {edges.height} edges, measured pool "
              f"{len(matcher.measured_pool)}, {len(matcher.strata)} strata")

    summary_rows: list[dict] = []
    null_draw_rows: list[dict] = []

    for code, (display, raw_id) in DISEASE_MAPPING.items():
        resolved = resolve_disease_id(raw_id, ontology)
        if resolved is None:
            warnings.warn(f"[{code}] disease ID '{raw_id}' not found; skipping")
            continue

        subtree = expand_subtree(resolved, ontology)
        approved = approved_targets(drug_targets, subtree)
        approved_in_trans = approved & trans_universe
        approved_in_cis = approved & cis_universe

        # Real coloc — same across conditions, take any.
        real_coloc = load_coloc_genes(gene_sets_dir, code, "Rest")
        if not real_coloc:
            warnings.warn(f"[{code}] no H4 coloc genes; skipping")
            continue
        n_coloc = len(real_coloc)
        print(f"[{code}] {display}: {n_coloc} coloc genes  |approved|={len(approved)}")

        for cond in CONDITIONS:
            c2t = cis_to_trans[cond]
            t2c = trans_to_cis[cond]

            # Observed (real) expansion + Fisher.
            real_regs, real_tgts = expand_pseudo(real_coloc, c2t, t2c)
            real_or_t, real_p_t, real_a_t, *_ = fisher_or(
                real_tgts, approved_in_trans, trans_universe
            )
            real_or_r, real_p_r, real_a_r, *_ = fisher_or(
                real_regs, approved_in_cis, cis_universe
            )

            # Null distributions.
            null_or_t = np.empty(args.n_perm, dtype=float)
            null_or_r = np.empty(args.n_perm, dtype=float)
            n_fallback = n_unmeasured = 0
            for k in range(args.n_perm):
                # Per-(perm, condition) stream, as in the S-LDSC null.
                pseudo, log = matchers[cond].sample(
                    real_coloc, perm_rng(args.seed, k + 1, cond))
                n_fallback += log["n_fallback"]
                n_unmeasured = log["n_unmeasured"]
                p_regs, p_tgts = expand_pseudo(pseudo, c2t, t2c)
                null_or_t[k] = fisher_or(p_tgts, approved_in_trans, trans_universe)[0]
                null_or_r[k] = fisher_or(p_regs, approved_in_cis, cis_universe)[0]

            # Empirical p-values (fraction of null ORs ≥ real OR).
            def emp_p(real: float, null: np.ndarray) -> float:
                if np.isnan(real):
                    return float("nan")
                # (n_ge + 1) / (N + 1) as a conservative estimator.
                return float((np.sum(null >= real) + 1) / (len(null) + 1))

            for gs_name, real_or, real_p, real_a, null in [
                ("Targets", real_or_t, real_p_t, real_a_t, null_or_t),
                ("Regulators", real_or_r, real_p_r, real_a_r, null_or_r),
            ]:
                summary_rows.append(
                    {
                        "disease": code,
                        "display": display,
                        "gene_set": gs_name,
                        "condition": cond,
                        "n_coloc": n_coloc,
                        "real_or": real_or,
                        "real_p": real_p,
                        "real_a": real_a,
                        "null_or_mean": float(np.mean(null)),
                        "null_or_p50": float(np.median(null)),
                        "null_or_p95": float(np.quantile(null, 0.95)),
                        "null_or_max": float(np.max(null)),
                        "empirical_p": emp_p(real_or, null),
                        "n_perm": args.n_perm,
                        "n_fallback_total": n_fallback,
                        "n_unmeasured": n_unmeasured,
                    }
                )
                if args.save_null_draws:
                    for draw_idx, val in enumerate(null):
                        null_draw_rows.append(
                            {
                                "disease": code,
                                "gene_set": gs_name,
                                "condition": cond,
                                "draw": draw_idx,
                                "null_or": float(val),
                            }
                        )

    out_csv = out_dir / "permutation_results.csv"
    pl.DataFrame(summary_rows).write_csv(out_csv)
    print(f"\nWrote {out_csv}")

    if args.save_null_draws and null_draw_rows:
        out_parquet = out_dir / "permutation_nulls.parquet"
        pl.DataFrame(null_draw_rows).write_parquet(out_parquet)
        print(f"Wrote {out_parquet}")


if __name__ == "__main__":
    main()
