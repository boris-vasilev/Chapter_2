#!/usr/bin/env python
"""EXPLORATORY: do quantitative, signed regulator->target edges carry
disease-specific signal that survives the CD4 cell-type confound?

Premise — the membership-based tests showed:
  (a) Fisher enrichment is reproduced by property/size-matched permutation
      (cell-type identity confound).
  (b) Diagonal/specificity test was weak (emp_p ~0.10-0.17) AND mis-specified
      because IMDs genuinely share drugs.

Membership throws away the QUANTITATIVE information. A random CD4 gene set
matches the membership distribution but has no reason to reproduce the SIGNED
LOGFC pattern of disease-specific regulator->coloc-target edges. So edge
weights might separate signal from confound.

This script:
  1. Builds per-disease directed edges (regulator -> coloc-target, signed logFC).
  2. Per disease, summarises each regulator's effect (mean signed logFC across
     that disease's targets downstream of it).
  3. For every disease pair, on SHARED regulators, asks: do the two diseases
     show concordant or divergent signed effects on their own (different)
     causal targets?
  4. Aggregate sign-agreement rate is compared to a label-permutation null
     that holds the DE network + regulator sets fixed and randomises which
     coloc targets each disease "claims".
  5. Checks the abundance/detectability confound (|logFC| vs trans-degree).
  6. Prints a blunt promising / not promising / underpowered verdict.

This is NOT publication-grade. It's a fast exploratory check.
"""

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import polars as pl


# Shared with approved_drug_enrichment_all.py — kept inline for one-file run.
IMD_DISEASES = [
    "AST", "ASTAO", "ASTCO", "ATD", "AITD",
    "CD", "IBD", "UC", "RA", "SLE", "T1D",
]
NON_IMD = ["IS", "MYOI", "OA", "T2D", "ADFH", "SCZ"]

# Diseases that genuinely share biology: we expect concordance WITHIN family
# (so within-family pairs are partly "true positives" for shared mechanism,
# not specificity); cross-family is where disease-specific signal should appear.
FAMILIES = {
    "atopic": {"AST", "ASTAO", "ASTCO", "ATD"},
    "IBD": {"CD", "IBD", "UC"},
    "endocrine_auto": {"T1D", "AITD"},
}


def family_of(d: str) -> tuple:
    return tuple(sorted(f for f, mem in FAMILIES.items() if d in mem))


def pair_relation(d1: str, d2: str) -> str:
    f1, f2 = set(family_of(d1)), set(family_of(d2))
    if f1 and f2 and (f1 & f2):
        return "same_family"
    in1, in2 = d1 in IMD_DISEASES, d2 in IMD_DISEASES
    if in1 and in2:
        return "different_IMD"
    if in1 != in2:
        return "IMD_vs_nonIMD"
    return "nonIMD_vs_nonIMD"


def load_gene_set(gene_sets_dir: Path, study: str, token: str, condition: str) -> set:
    genes: set = set()
    for file in gene_sets_dir.glob(f"{study}_*_{condition}_{token}.txt"):
        with open(file) as fh:
            for line in fh:
                g = line.strip()
                if g:
                    genes.add(g)
    return genes


def regulator_profile_fast(
    regs: set, coloc: set, cis_to_targets: dict
) -> dict:
    """For each regulator R, mean signed logFC across (R -> t in coloc) edges.
    Uses precomputed cis_to_targets[r] = list[(trans, logFC)] for speed."""
    profile = {}
    for r in regs:
        targets = cis_to_targets.get(r)
        if not targets:
            continue
        lfcs = [lfc for t, lfc in targets if t in coloc]
        if lfcs:
            profile[r] = float(np.mean(lfcs))
    return profile


def sign_agree_rate(profiles: dict) -> tuple[float, int]:
    """Across all disease pairs, fraction of shared regulators whose signed
    mean-logFC matches in both diseases. Returns (rate, n_triples)."""
    match, total = 0, 0
    diseases = list(profiles)
    for ii in range(len(diseases)):
        p1 = profiles[diseases[ii]]
        if not p1:
            continue
        for jj in range(ii + 1, len(diseases)):
            p2 = profiles[diseases[jj]]
            if not p2:
                continue
            for r in p1.keys() & p2.keys():
                v1, v2 = p1[r], p2[r]
                if v1 == 0 or v2 == 0:
                    continue
                total += 1
                if (v1 > 0) == (v2 > 0):
                    match += 1
    return (match / total if total > 0 else float("nan")), total


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-sets-dir", required=True)
    parser.add_argument("--de-results", required=True)
    parser.add_argument(
        "--condition", default="Stim48hr",
        help="Default Stim48hr (where ATD/UC signal lives).",
    )
    parser.add_argument(
        "--target-set", default="H4_genes",
        choices=["H4_genes", "targets_H4"],
        help="Which set to use as causal targets. Default = coloc genes.",
    )
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--n-perm", type=int, default=500)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    gene_sets_dir = Path(args.gene_sets_dir)
    rng = np.random.default_rng(args.seed)

    # ---- 1. Load DE for the chosen condition ------------------------------
    print(f"Loading DE for condition={args.condition}...")
    de = pl.read_csv(args.de_results).filter(
        pl.col("culture_condition") == args.condition
    )
    print(
        f"  {de.height} DE rows, "
        f"{de.get_column('cis_ensg').n_unique()} unique cis, "
        f"{de.get_column('trans_ensg').n_unique()} unique trans"
    )

    # Precompute cis -> list of (trans, logFC). Hot path for permutations.
    cis_to_targets: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for r in de.iter_rows(named=True):
        c, t, lfc = r.get("cis_ensg"), r.get("trans_ensg"), r.get("logFC")
        if c and t and lfc is not None:
            cis_to_targets[c].append((t, float(lfc)))
    trans_universe = sorted({g for g in de.get_column("trans_ensg").unique().to_list() if g})

    # trans-degree (detectability proxy)
    trans_degree: dict[str, int] = defaultdict(int)
    for tlist in cis_to_targets.values():
        for t, _ in tlist:
            trans_degree[t] += 1

    # ---- 2. Build per-disease edges ---------------------------------------
    all_diseases = IMD_DISEASES + NON_IMD
    disease_regs: dict[str, set] = {}
    disease_coloc: dict[str, set] = {}
    edge_rows: list[dict] = []

    for d in all_diseases:
        regs = load_gene_set(gene_sets_dir, d, "regulators_H4", args.condition)
        coloc = load_gene_set(gene_sets_dir, d, args.target_set, args.condition)
        if not regs or not coloc:
            continue
        disease_regs[d] = regs
        disease_coloc[d] = coloc
        for r in regs:
            for t, lfc in cis_to_targets.get(r, []):
                if t in coloc:
                    edge_rows.append({
                        "disease": d, "cis_ensg": r, "trans_ensg": t,
                        "logFC": lfc, "sign": int(np.sign(lfc)),
                    })

    if not edge_rows:
        print("\nNo edges across any disease. Stop.")
        return
    edges_df = pl.DataFrame(edge_rows)
    edges_df.write_csv(out_dir / "edges_per_disease.csv")
    print(f"\nWrote {out_dir / 'edges_per_disease.csv'}  ({edges_df.height} edges)")

    edge_counts = {
        d: edges_df.filter(pl.col("disease") == d).height for d in disease_regs
    }
    print("Edge counts per disease (target_set =", args.target_set, "):")
    for d, n in sorted(edge_counts.items(), key=lambda kv: -kv[1]):
        print(f"  {d}: {n:5d} edges  ({len(disease_regs[d])} regs, "
              f"{len(disease_coloc[d])} coloc)")

    # ---- 3. Regulator profiles & observed sign-agree ----------------------
    profiles = {
        d: regulator_profile_fast(disease_regs[d], disease_coloc[d], cis_to_targets)
        for d in disease_regs
    }
    observed_sign_agree, n_triples = sign_agree_rate(profiles)
    print(f"\nObserved sign-agree across all shared-regulator triples: "
          f"{observed_sign_agree:.3f}  (n={n_triples})")

    # ---- 4. Pairwise concordance (per-disease-pair detail) ----------------
    pairwise_rows = []
    diseases_with_profiles = list(profiles)
    for ii in range(len(diseases_with_profiles)):
        d1 = diseases_with_profiles[ii]
        p1 = profiles[d1]
        for jj in range(ii + 1, len(diseases_with_profiles)):
            d2 = diseases_with_profiles[jj]
            p2 = profiles[d2]
            shared = list((p1.keys() & p2.keys()))
            n_shared = len(shared)
            if n_shared == 0:
                pairwise_rows.append({
                    "d1": d1, "d2": d2, "n_shared_reg": 0,
                    "sign_agree": float("nan"), "pearson_r": float("nan"),
                    "relation": pair_relation(d1, d2),
                })
                continue
            v1 = np.array([p1[r] for r in shared])
            v2 = np.array([p2[r] for r in shared])
            sa = float(np.mean(np.sign(v1) == np.sign(v2)))
            if n_shared >= 3 and np.std(v1) > 0 and np.std(v2) > 0:
                pr = float(np.corrcoef(v1, v2)[0, 1])
            else:
                pr = float("nan")
            pairwise_rows.append({
                "d1": d1, "d2": d2, "n_shared_reg": n_shared,
                "sign_agree": sa, "pearson_r": pr,
                "relation": pair_relation(d1, d2),
            })
    pw_df = pl.DataFrame(pairwise_rows)
    pw_df.write_csv(out_dir / "pairwise_concordance.csv")
    print(f"Wrote {out_dir / 'pairwise_concordance.csv'}")

    # Family-grouped summary
    rel_summary = (
        pw_df.filter(pl.col("n_shared_reg") >= 2)
        .group_by("relation")
        .agg(
            pl.len().alias("n_pairs"),
            pl.col("n_shared_reg").mean().alias("mean_n_shared"),
            pl.col("sign_agree").mean().alias("mean_sign_agree"),
            pl.col("pearson_r").mean().alias("mean_pearson_r"),
        )
    )
    rel_summary.write_csv(out_dir / "concordance_by_relation.csv")
    print("\nFamily-grouped concordance summary:")
    print(rel_summary)

    # ---- 5. Label-permutation null ----------------------------------------
    # Hold: DE network, regulator sets, disease structure fixed.
    # Randomise: which coloc targets each disease claims (sample from
    # trans_universe at original sizes).
    print(f"\nRunning {args.n_perm} label permutations (this is the slow step)...")
    trans_arr = np.array(trans_universe)
    null_rates = np.empty(args.n_perm, dtype=float)
    for k in range(args.n_perm):
        null_profiles: dict[str, dict[str, float]] = {}
        for d, regs in disease_regs.items():
            n = len(disease_coloc[d])
            fake = set(rng.choice(trans_arr, size=n, replace=False).tolist())
            null_profiles[d] = regulator_profile_fast(regs, fake, cis_to_targets)
        rate, _ = sign_agree_rate(null_profiles)
        null_rates[k] = rate
        if (k + 1) % max(1, args.n_perm // 10) == 0:
            print(f"  ... {k + 1}/{args.n_perm}")

    finite = null_rates[~np.isnan(null_rates)]
    null_mean = float(np.mean(finite))
    null_p5 = float(np.percentile(finite, 5))
    null_p95 = float(np.percentile(finite, 95))
    # Two-sided: observed could be below (divergent) or above (concordant).
    emp_p_below = float((np.sum(finite <= observed_sign_agree) + 1) / (len(finite) + 1))
    emp_p_above = float((np.sum(finite >= observed_sign_agree) + 1) / (len(finite) + 1))
    print(f"\nNull sign-agree: mean={null_mean:.3f}, "
          f"5%={null_p5:.3f}, 95%={null_p95:.3f}  (n_perm={len(finite)})")
    print(f"Observed         : {observed_sign_agree:.3f}")
    print(f"emp_p (observed <= null, i.e. divergent): {emp_p_below:.4f}")
    print(f"emp_p (observed >= null, i.e. concordant): {emp_p_above:.4f}")

    pl.DataFrame({
        "null_sign_agree": null_rates,
    }).write_csv(out_dir / "label_perm_null.csv")

    # ---- 6. Abundance / detectability check -------------------------------
    print("\n--- Abundance / detectability check ---")
    abun_rows = []
    for d in disease_regs:
        sub = edges_df.filter(pl.col("disease") == d)
        if sub.height < 5:
            continue
        ts = sub.get_column("trans_ensg").to_list()
        lfcs = np.abs(sub.get_column("logFC").to_numpy())
        degs = np.array([trans_degree.get(t, 0) for t in ts], dtype=float)
        if np.std(degs) == 0 or np.std(lfcs) == 0:
            continue
        r = float(np.corrcoef(np.log1p(degs), lfcs)[0, 1])
        abun_rows.append({"disease": d, "n_edges": sub.height, "corr_logdeg_abslfc": r})
    abun_df = pl.DataFrame(abun_rows)
    abun_df.write_csv(out_dir / "abundance_check.csv")
    print(abun_df)
    mean_abun = (
        float(abun_df.get_column("corr_logdeg_abslfc").mean())
        if abun_df.height > 0 else float("nan")
    )
    print(f"Mean per-disease corr(log(trans_degree), |logFC|) = {mean_abun:.3f}")

    # ---- 7. Blunt verdict -------------------------------------------------
    total_edges = edges_df.height
    median_edges = int(np.median(list(edge_counts.values()))) if edge_counts else 0
    n_pairs_useful = pw_df.filter(pl.col("n_shared_reg") >= 2).height

    print("\n" + "=" * 64)
    print("BLUNT SUMMARY")
    print("=" * 64)
    print(f"(i)  Edge counts: total={total_edges}  "
          f"median per disease={median_edges}  diseases with edges={len(disease_regs)}")
    print(f"     Pairs with >=2 shared regulators: {n_pairs_useful} / {pw_df.height}")
    print(f"(ii) Observed sign-agree (across all disease pairs, shared regs): "
          f"{observed_sign_agree:.3f}  (n_triples={n_triples})")
    if rel_summary.height > 0:
        for row in rel_summary.iter_rows(named=True):
            print(f"     {row['relation']:>18s}: n_pairs={row['n_pairs']:3d}  "
                  f"mean_sign_agree={row['mean_sign_agree']:.3f}  "
                  f"mean_pearson_r={row['mean_pearson_r']:.3f}")
    print(f"(iii) Null mean={null_mean:.3f}, observed={observed_sign_agree:.3f}")
    print(f"      Lower-tail emp_p (divergent: specific) = {emp_p_below:.4f}")
    print(f"      Upper-tail emp_p (concordant: shared)  = {emp_p_above:.4f}")
    print(f"(iv) Abundance confound: mean corr(log_degree, |logFC|) = {mean_abun:.3f}")
    print(f"(v)  Power: median edges/disease = {median_edges}; "
          f"useful pairs = {n_pairs_useful}")

    print("\nVERDICT:")
    underpowered = (total_edges < 50) or (n_pairs_useful < 5) or (median_edges < 5)
    if underpowered:
        print("  UNDERPOWERED — sparse edges and/or few shared regulators across diseases.")
        print("  Treat any nominal signal here as suggestive only. A larger Perturb-seq")
        print("  screen or relaxed coloc threshold may be needed to test this properly.")
    elif emp_p_below < 0.05:
        print("  PROMISING — signed regulator->target effects are MORE DIVERGENT across")
        print("  diseases than label-permutation expects. Edge weights carry disease-")
        print("  specific quantitative signal beyond shared regulator membership.")
        print("  -> worth a proper follow-up: identify which regulators drive divergence,")
        print("     re-test on independent perturbation data, condition on expression bins.")
    elif emp_p_above < 0.05:
        print("  NOT PROMISING — shared regulators show MORE CONCORDANT signed effects")
        print("  than null. The quantitative layer mirrors shared structure; same confound")
        print("  as membership. Edge weights do NOT separate signal from cell-type identity.")
    else:
        print("  INCONCLUSIVE — observed sign-agreement is indistinguishable from the")
        print("  label-permutation null. No evidence for disease-specific edge signal at")
        print("  this sample size. Reconsider only if power can be increased.")

    if abs(mean_abun) > 0.3:
        print("  WARNING: |logFC| correlates with detectability proxy "
              f"(|r|={abs(mean_abun):.2f}). Any signal needs to be reproduced within")
        print("  expression-matched bins before being taken seriously.")
    else:
        print("  Abundance confound looks weak — unlikely to be driving any signal.")


if __name__ == "__main__":
    main()
