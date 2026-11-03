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
  2. Build a "pseudo-coloc" set by sampling size-matched genes from a
     configurable background pool (default: union of cis and trans universes).
  3. Expand the pseudo-coloc set into "pseudo-regulators" and "pseudo-targets"
     using the same Perturb-seq DE relationships as the real expansion
     (precomputed cis→trans and trans→cis indices from --de-results).
  4. Re-run the Fisher exact test against the same approved-drug-target set.
  5. Repeat N times to get an empirical null distribution of odds ratios.
  6. Emit empirical p-values = fraction of null ORs ≥ real OR.

Optional matching beyond size:

  --match-by perturb_degree   Stratify the sampling pool by how many trans
                              genes each cis gene DEs (a proxy for "is this a
                              perturbation hub?"). Controls for the possibility
                              that real coloc genes happen to be more
                              perturbation-active than random.
  --match-by expr_decile      Stratify by mean DE base expression in the DE
                              table (a proxy for expression level in the
                              screen). Controls for "real coloc genes are
                              highly expressed in CD4+ T cells".

Outputs:
  <out>/permutation_results.csv  one row per (disease, gene_set, condition):
      real_or, real_p, real_a, null_or_mean, null_or_p50, null_or_p95,
      null_or_max, empirical_p, n_perm
  <out>/permutation_nulls.parquet (optional, with --save-null-draws)  raw null
      OR for every (disease, gene_set, condition, draw).
"""

import argparse
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import fisher_exact

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


def build_de_indices(
    de: pl.DataFrame,
) -> tuple[dict[str, dict[str, set[str]]], dict[str, dict[str, set[str]]]]:
    """Per-condition cis→trans and trans→cis maps from the DE table.

    cis_to_trans[cond][cis_ensg] = set of trans_ensg DE'd when cis is perturbed.
    trans_to_cis[cond][trans_ensg] = set of cis_ensg whose perturbation DE'd trans.
    """
    cis_to_trans: dict[str, dict[str, set[str]]] = {
        c: defaultdict(set) for c in CONDITIONS
    }
    trans_to_cis: dict[str, dict[str, set[str]]] = {
        c: defaultdict(set) for c in CONDITIONS
    }
    for r in de.iter_rows(named=True):
        cond = r.get("culture_condition")
        cis = r.get("cis_ensg")
        trans = r.get("trans_ensg")
        if not (cond and cis and trans):
            continue
        if cond not in cis_to_trans:
            continue
        cis_to_trans[cond][cis].add(trans)
        trans_to_cis[cond][trans].add(cis)
    return cis_to_trans, trans_to_cis


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


def build_match_strata(
    pool: list[str],
    cis_to_trans_any: dict[str, set[str]],
    de_meanbase: dict[str, float],
    match_by: list[str],
    n_bins: int = 5,
) -> dict[tuple, list[str]]:
    """Bin the sampling pool into strata so the pseudo-coloc sample can match
    real-coloc gene properties beyond just size.

    Each gene gets a stratum key built from the selected attributes; we then
    sample within-stratum to match the real coloc's stratum distribution.
    """
    if not match_by:
        return {(): list(pool)}

    # Precompute per-gene stratum keys.
    keys: dict[str, list] = {}
    for g in pool:
        key: list = []
        if "perturb_degree" in match_by:
            key.append(len(cis_to_trans_any.get(g, set())))
        if "expr_decile" in match_by:
            key.append(de_meanbase.get(g, 0.0))
        keys[g] = key

    # Bin each attribute independently using percentile cuts over the pool.
    if not keys:
        return {(): list(pool)}
    n_attrs = len(next(iter(keys.values())))
    bin_edges: list[np.ndarray] = []
    for i in range(n_attrs):
        vals = np.array([keys[g][i] for g in pool], dtype=float)
        if np.allclose(vals.max(), vals.min()):
            bin_edges.append(np.array([vals.min(), vals.max() + 1e-9]))
        else:
            bin_edges.append(np.quantile(vals, np.linspace(0, 1, n_bins + 1)))

    strata: dict[tuple, list[str]] = defaultdict(list)
    for g in pool:
        key = []
        for i, v in enumerate(keys[g]):
            idx = max(0, min(n_bins - 1, int(np.searchsorted(bin_edges[i][1:-1], v))))
            key.append(idx)
        strata[tuple(key)].append(g)
    return strata


def sample_pseudo_coloc(
    real_coloc: set[str],
    pool: list[str],
    strata: dict[tuple, list[str]],
    cis_to_trans_any: dict[str, set[str]],
    de_meanbase: dict[str, float],
    match_by: list[str],
    n_bins: int,
    rng: np.random.Generator,
) -> set[str]:
    """Sample a pseudo-coloc set the same size as `real_coloc`, respecting any
    matching strata. Falls back to overall-pool sampling for strata where the
    pool is too small."""
    if not match_by:
        n = min(len(real_coloc), len(pool))
        return set(rng.choice(pool, size=n, replace=False))

    # For each real coloc gene, find its stratum, then sample one gene from
    # that stratum (without replacement across all draws).
    drawn: set[str] = set()
    for g in real_coloc:
        # Re-derive g's stratum key the same way as the pool.
        # (If g isn't in the pool, fall back to a random sample from any stratum.)
        key_vals = []
        if "perturb_degree" in match_by:
            key_vals.append(len(cis_to_trans_any.get(g, set())))
        if "expr_decile" in match_by:
            key_vals.append(de_meanbase.get(g, 0.0))
        # Find best matching stratum by minimum L1 distance on raw key
        # — simple and robust to genes not in the original pool binning.
        best_key = None
        best_dist = float("inf")
        for k in strata:
            d = sum(abs(a - b) for a, b in zip(k, key_vals))
            if d < best_dist and strata[k]:
                best_dist = d
                best_key = k
        if best_key is None:
            continue
        candidates = [x for x in strata[best_key] if x not in drawn]
        if not candidates:
            # Stratum exhausted; fall back to global pool minus drawn.
            candidates = [x for x in pool if x not in drawn]
            if not candidates:
                break
        drawn.add(rng.choice(candidates))
    return drawn


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
    parser.add_argument(
        "--pool",
        choices=["trans", "cis", "intersection", "union"],
        default="union",
        help=(
            "Background pool to sample pseudo-coloc genes from. "
            "'trans'=detected genes; 'cis'=perturbed genes; "
            "'intersection'=genes both perturbed and detected; "
            "'union'=either."
        ),
    )
    parser.add_argument(
        "--match-by",
        nargs="*",
        default=[],
        choices=["perturb_degree", "expr_decile"],
        help="Stratify sampling so pseudo-coloc matches real coloc on these.",
    )
    parser.add_argument("--n-bins", type=int, default=5)
    parser.add_argument(
        "--min-abs-logfc",
        type=float,
        default=0.0,
        help="Keep only DE edges with |logFC| >= this when expanding coloc "
        "genes into regulators/targets (real AND null), matching "
        "make_gene_sets_coloc.py and the S-LDSC permutation. Self-edges "
        "(cis==trans) are excluded when the filter is active. The cis/trans "
        "universes and sampling pool stay threshold-independent; only the "
        "expansion graph and the perturb-degree used for matching are filtered. "
        "Default 0.0 = no filter (original behaviour).",
    )
    parser.add_argument(
        "--save-null-draws",
        action="store_true",
        help="Write the raw null OR for every permutation as a parquet file.",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    gene_sets_dir = Path(args.gene_sets_dir)
    rng = np.random.default_rng(args.seed)

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

    # We need the DE results both for cis/trans universes AND for the indices.
    de_cols = ["cis_ensg", "trans_ensg", "culture_condition"]
    if "expr_decile" in args.match_by:
        de_cols.append("base_DESeq2")
    if args.min_abs_logfc > 0:
        de_cols.append("logFC")
    de = pl.read_csv(args.de_results, columns=de_cols)

    cis_universe = {g for g in de.get_column("cis_ensg").unique().to_list() if g}
    trans_universe = {g for g in de.get_column("trans_ensg").unique().to_list() if g}
    print(f"  |cis|={len(cis_universe)}  |trans|={len(trans_universe)}")

    if args.pool == "trans":
        pool = list(trans_universe)
    elif args.pool == "cis":
        pool = list(cis_universe)
    elif args.pool == "intersection":
        pool = list(trans_universe & cis_universe)
    else:  # union
        pool = list(trans_universe | cis_universe)
    print(f"  sampling pool ({args.pool}) size={len(pool)}")

    # Expansion edges: optionally filter by |logFC| (and drop self-edges) so the
    # real and null expansions match make_gene_sets_coloc.py / the S-LDSC
    # permutation. The universes and pool above are intentionally left unfiltered
    # (threshold-independent background), exactly as in permutation_sweep.py.
    de_edges = de
    if args.min_abs_logfc > 0:
        de_edges = de.filter(
            (pl.col("logFC").abs() >= args.min_abs_logfc)
            & (pl.col("cis_ensg") != pl.col("trans_ensg"))
        )
        print(
            f"  |logFC| >= {args.min_abs_logfc}: kept {de_edges.height}/"
            f"{de.height} DE edges for expansion (self-edges excluded)"
        )
    cis_to_trans, trans_to_cis = build_de_indices(de_edges)
    # Condition-agnostic perturb_degree (sum over conditions) for stratification.
    cis_to_trans_any: dict[str, set[str]] = defaultdict(set)
    for cond_map in cis_to_trans.values():
        for cis, ts in cond_map.items():
            cis_to_trans_any[cis] |= ts
    # Per-gene mean base expression if available.
    de_meanbase: dict[str, float] = {}
    if "expr_decile" in args.match_by:
        agg = (
            de.group_by("trans_ensg")
            .agg(pl.col("base_DESeq2").mean().alias("mean_base"))
            .rename({"trans_ensg": "ensg"})
        )
        de_meanbase = dict(
            zip(agg.get_column("ensg").to_list(), agg.get_column("mean_base").to_list())
        )

    strata = build_match_strata(
        pool, cis_to_trans_any, de_meanbase, args.match_by, args.n_bins
    )
    print(f"  strata: {len(strata)} bins")

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
            for k in range(args.n_perm):
                pseudo = sample_pseudo_coloc(
                    real_coloc,
                    pool,
                    strata,
                    cis_to_trans_any,
                    de_meanbase,
                    args.match_by,
                    args.n_bins,
                    rng,
                )
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
                        "min_abs_logfc": args.min_abs_logfc,
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

        # --- Coloc anchor (positive control) ---
        # Drug-target enrichment of the coloc gene set ITSELF (no DE expansion,
        # so condition-independent) vs the SAME size/property-matched random
        # anchor null, scored on the sampling-pool universe. Mirrors the H4
        # anchor row in the S-LDSC permutation. Drawn after the condition loop
        # so the regulator/target null draws above are byte-identical to before.
        pool_set = set(pool)
        approved_in_pool = approved & pool_set
        real_or_c, real_p_c, real_a_c, *_ = fisher_or(
            real_coloc, approved_in_pool, pool_set
        )
        null_or_c = np.empty(args.n_perm, dtype=float)
        for k in range(args.n_perm):
            pseudo = sample_pseudo_coloc(
                real_coloc,
                pool,
                strata,
                cis_to_trans_any,
                de_meanbase,
                args.match_by,
                args.n_bins,
                rng,
            )
            null_or_c[k] = fisher_or(pseudo, approved_in_pool, pool_set)[0]
        emp_p_c = (
            float("nan")
            if np.isnan(real_or_c)
            else float((np.sum(null_or_c >= real_or_c) + 1) / (args.n_perm + 1))
        )
        # One row per condition (identical values) to keep the grid complete.
        for cond in CONDITIONS:
            summary_rows.append(
                {
                    "disease": code,
                    "display": display,
                    "gene_set": "Coloc",
                    "condition": cond,
                    "n_coloc": n_coloc,
                    "real_or": real_or_c,
                    "real_p": real_p_c,
                    "real_a": real_a_c,
                    "null_or_mean": float(np.mean(null_or_c)),
                    "null_or_p50": float(np.median(null_or_c)),
                    "null_or_p95": float(np.quantile(null_or_c, 0.95)),
                    "null_or_max": float(np.max(null_or_c)),
                    "empirical_p": emp_p_c,
                    "n_perm": args.n_perm,
                    "min_abs_logfc": args.min_abs_logfc,
                }
            )
        if args.save_null_draws:
            for draw_idx, val in enumerate(null_or_c):
                null_draw_rows.append(
                    {
                        "disease": code,
                        "gene_set": "Coloc",
                        "condition": "all",
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
