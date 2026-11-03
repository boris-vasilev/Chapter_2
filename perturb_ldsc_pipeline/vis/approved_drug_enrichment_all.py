#!/usr/bin/env python
"""Per-disease approved-drug-target enrichment of coloc gene sets.

For each study code in DISEASE_MAPPING, runs the same Fisher-exact enrichment as
approved_drug_enrichment.py:
  - gene-set side: H4_genes / regulators_H4 / targets_H4 from data/gene_sets/<dir>
  - background:    cis_ensg (for regulators), trans_ensg (everything else),
                   from the Perturb-seq DE results
  - positive set:  targets with an approved or Phase II+ drug for the disease,
                   propagated through the OT disease ontology (descendants)

Writes:
  <out-dir>/enrichment_all.csv  — one row per (disease, gene_set, condition)
  <out-dir>/enrichment_all.png  — heatmap of log2(OR) with significance asterisks,
                                  IMDs grouped left of non-IMDs for specificity.
"""

import argparse
import warnings
from pathlib import Path

import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from matplotlib.patches import Rectangle
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

APPROVED_STAGES = [
    "APPROVAL",
    "PREAPPROVAL",
    "PHASE_4",
    "PHASE_3",
    "PHASE_2_3",
    "PHASE_2",
]

IMD_DISEASES = [
    "AST",
    "ASTAO",
    "ASTCO",
    "ASPH",
    "ATD",
    "AITD",
    "CD",
    "IBD",
    "UC",
    "RA",
    "SLE",
    "T1D",
]
NON_IMD = ["IS", "MYOI", "OA", "T2D", "ADFH", "SCZ"]

# Short display labels for heatmap axes — compact enough to read at 45° rotation.
SHORT_NAMES: dict[str, str] = {
    "AST": "Asthma",
    "ASTAO": "Asthma (AO)",
    "ASTCO": "Asthma (CO)",
    "ATD": "Atopic derm.",
    "AITD": "Autoimm. thyroid",
    "CD": "Crohn's",
    "IBD": "IBD",
    "UC": "Ulcer. colitis",
    "RA": "Rheum. arthritis",
    "SLE": "SLE",
    "T1D": "T1D",
    "IS": "Ischemic stroke",
    "MYOI": "Myocard. infarct.",
    "OA": "Osteoarthritis",
    "T2D": "T2D",
    "ADFH": "Alzheimer's",
    "SCZ": "Schizophrenia",
}

# Study code -> (display name, OT disease ID).
# IDs are best-guess defaults; the script falls back to dbXRefs resolution so
# either MONDO_xxx or EFO_xxx forms typically work. Anything still unresolvable
# is skipped with a warning so you can fix the mapping and re-run.
DISEASE_MAPPING: dict[str, tuple[str, str]] = {
    # IMDs
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
    # Non-IMDs
    "IS": ("Ischemic stroke", "HP_0002140"),
    "MYOI": ("Myocardial infarction", "EFO_0000612"),
    "OA": ("Osteoarthritis", "MONDO_0005178"),
    "T2D": ("Type 2 diabetes", "MONDO_0005148"),
    "ADFH": ("Alzheimer's disease", "MONDO_0004975"),
    "SCZ": ("Schizophrenia", "MONDO_0005090"),
}


def load_disease_groups(path: str) -> None:
    """Override the hardcoded disease metadata from a disease_groups.csv.

    Such a CSV is emitted by make_gene_sets_tenk10k.py with columns:
    study, display_name, short_label, group ("IMD"/"other"), ot_disease_id.
    Rebuilds DISEASE_MAPPING, IMD_DISEASES, NON_IMD and SHORT_NAMES in place
    so the heatmap separation works for the tenk10k diseases.
    """
    global DISEASE_MAPPING, IMD_DISEASES, NON_IMD, SHORT_NAMES
    df = pl.read_csv(path)
    DISEASE_MAPPING = {
        r["study"]: (r["display_name"], r["ot_disease_id"])
        for r in df.iter_rows(named=True)
        if r["ot_disease_id"]
    }
    IMD_DISEASES = [r["study"] for r in df.iter_rows(named=True) if r["group"] == "IMD"]
    NON_IMD = [r["study"] for r in df.iter_rows(named=True) if r["group"] != "IMD"]
    SHORT_NAMES = {r["study"]: r["short_label"] for r in df.iter_rows(named=True)}


def resolve_disease_id(query_id: str, ontology: pl.DataFrame) -> str | None:
    """Resolve a disease identifier to the OT canonical ID.

    Tries: (1) direct match on the `id` column; (2) cross-references in
    `dbXRefs` in both underscore and colon notation.
    """
    if ontology.filter(pl.col("id") == query_id).height > 0:
        return query_id
    for candidate in (query_id, query_id.replace("_", ":")):
        match = ontology.filter(pl.col("dbXRefs").list.contains(candidate))
        if match.height > 0:
            return match.get_column("id").item(0)
    return None


def expand_subtree(disease_id: str, ontology: pl.DataFrame) -> set[str]:
    matches = ontology.filter(pl.col("id") == disease_id)
    if matches.height == 0:
        return {disease_id}
    # `descendants` is Array<String>; .to_list() yields a Python list-of-lists,
    # so the first row's value is either a list[str] or None.
    descendants = matches.get_column("descendants").to_list()[0]
    return {disease_id} | set(descendants or [])


def approved_targets(drug_targets: pl.DataFrame, subtree: set[str]) -> set[str]:
    return set(
        drug_targets.filter(pl.col("diseaseId").is_in(subtree))
        .filter(pl.col("clinicalStage").is_in(APPROVED_STAGES))
        .get_column("targetId")
        .unique()
        .to_list()
    )


def load_gene_sets(gene_sets_dir: Path, study: str) -> pl.DataFrame | None:
    # `{study}_*` rather than `{study}*` to avoid prefix collisions
    # (e.g. AST matching ASTAO/ASTCO).
    matches = list(gene_sets_dir.glob(f"{study}_*"))
    if not matches:
        return None
    rows = []
    for file in matches:
        parts = file.stem.split("_")
        gene_set = "_".join(parts[-2:])
        condition = parts[-3]
        with open(file) as f:
            for gene in map(str.strip, f):
                if gene:
                    rows.append(
                        {"gene": gene, "gene_set": gene_set, "condition": condition}
                    )
    return pl.DataFrame(rows) if rows else None


def universe_for(
    gene_set: str, cis_universe: set[str], trans_universe: set[str]
) -> set[str]:
    if gene_set.startswith("regulators"):
        return cis_universe
    return trans_universe


def run_tests(
    gene_sets_df: pl.DataFrame,
    approved_for_disease: set[str],
    cis_universe: set[str],
    trans_universe: set[str],
) -> list[dict]:
    out: list[dict] = []
    combos = gene_sets_df.select(["gene_set", "condition"]).unique().iter_rows()
    for gs, cond in combos:
        universe = universe_for(gs, cis_universe, trans_universe)
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
        out.append(
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
    return out


def plot_heatmap(df: pl.DataFrame, out_path: Path) -> None:
    condition_order = ["Rest", "Stim8hr", "Stim48hr"]
    # Coloc genes don't change with culture condition, so we only show one row
    # (Rest is chosen as the canonical key, but any condition would give the
    # same numbers). Regulators and Targets each get all three conditions.
    row_keys = [("coloc genes", "Rest")]
    row_labels = ["coloc genes"]
    for gs in ["Regulators", "Targets"]:
        for cond in condition_order:
            row_keys.append((gs, cond))
            row_labels.append(f"{gs} | {cond}")

    present = df.get_column("disease").unique().to_list()
    disease_order = [d for d in IMD_DISEASES + NON_IMD if d in present]

    or_matrix = np.full((len(row_keys), len(disease_order)), np.nan)
    p_matrix = np.full_like(or_matrix, np.nan)
    fdr_matrix = np.full_like(or_matrix, np.nan)

    for i, (gs, cond) in enumerate(row_keys):
        for j, disease in enumerate(disease_order):
            sub = df.filter(
                (pl.col("gene_set") == gs)
                & (pl.col("condition") == cond)
                & (pl.col("disease") == disease)
            )
            if sub.height > 0:
                or_matrix[i, j] = sub.get_column("odds_ratio").item(0)
                p_matrix[i, j] = sub.get_column("p_value").item(0)
                fdr_matrix[i, j] = sub.get_column("fdr_bh").item(0)

    # Cap OR for colour scaling; treat inf as upper cap.
    capped = np.where(np.isinf(or_matrix), 16.0, or_matrix)
    log2_or = np.log2(np.clip(capped, 0.25, 16))

    fig_w = max(10, 0.85 * len(disease_order) + 4)

    fig, ax = plt.subplots(figsize=(fig_w, 8))
    im = ax.imshow(log2_or, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")

    for i in range(len(row_keys)):
        for j in range(len(disease_order)):
            p = p_matrix[i, j]
            if np.isnan(p):
                ax.text(j, i, "·", ha="center", va="center", color="grey", fontsize=18)
                continue
            stars = (
                "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
            )
            if stars:
                txt = ax.text(
                    j,
                    i,
                    stars,
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=20,
                    fontweight="bold",
                )
                txt.set_path_effects(
                    [
                        path_effects.Stroke(linewidth=2.5, foreground="black"),
                        path_effects.Normal(),
                    ]
                )

    # Outline cells passing BH-FDR < 0.05 with a thick gold border.
    for i in range(len(row_keys)):
        for j in range(len(disease_order)):
            f = fdr_matrix[i, j]
            if not np.isnan(f) and f < 0.05:
                ax.add_patch(
                    Rectangle(
                        (j - 0.5, i - 0.5),
                        1.0,
                        1.0,
                        fill=False,
                        edgecolor="#FFD700",
                        linewidth=3.0,
                    )
                )

    n_imd = sum(1 for d in disease_order if d in IMD_DISEASES)
    if 0 < n_imd < len(disease_order):
        ax.axvline(n_imd - 0.5, color="black", linewidth=2.0)
        ax.text(
            n_imd / 2 - 0.5,
            -0.85,
            "Immune-mediated diseases",
            ha="center",
            fontsize=18,
            fontweight="bold",
        )
        ax.text(
            (n_imd + len(disease_order)) / 2 - 0.5,
            -0.85,
            "Other diseases",
            ha="center",
            fontsize=18,
            fontweight="bold",
        )

    ax.set_xticks(range(len(disease_order)))
    ax.set_xticklabels(
        [SHORT_NAMES.get(d, d) for d in disease_order],
        rotation=45,
        ha="right",
        fontsize=18,
    )
    ax.set_yticks(range(len(row_keys)))
    ax.set_yticklabels(row_labels, fontsize=18)

    cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("log2(OR), capped at ±4", fontsize=18)
    cbar.ax.tick_params(labelsize=12)

    ax.set_title(
        "Phase 2+ Drug target enrichment of coloc-Perturb-seq gene sets across diseases\n"
        "Fisher's exact test: * p<0.05, ** p<0.01, *** p<0.001; "
        "gold border: BH-FDR < 0.05; · = no data",
        pad=40,
        fontsize=18,
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=400, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-sets-dir", required=True)
    parser.add_argument("--drug-targets", required=True)
    parser.add_argument("--disease-ontology", required=True)
    parser.add_argument("--de-results", required=True)
    parser.add_argument(
        "--disease-groups",
        default=None,
        help="Optional disease_groups.csv (from make_gene_sets_tenk10k.py) to "
        "override the built-in disease mapping / IMD split.",
    )
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    if args.disease_groups:
        load_disease_groups(args.disease_groups)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    gene_sets_dir = Path(args.gene_sets_dir)

    print("Loading ontology, drug-targets, and Perturb-seq DE results...")
    ontology = (
        pl.scan_parquet(args.disease_ontology)
        .select(["id", "name", "descendants", "dbXRefs"])
        .collect()
    )
    drug_targets = (
        pl.scan_parquet(args.drug_targets)
        .select(["targetId", "diseaseId", "clinicalStage"])
        .collect()
    )
    de = pl.read_csv(args.de_results, columns=["cis_ensg", "trans_ensg"])
    cis_universe = {g for g in de.get_column("cis_ensg").unique().to_list() if g}
    trans_universe = {g for g in de.get_column("trans_ensg").unique().to_list() if g}
    print(
        f"  |cis_universe|={len(cis_universe)}, |trans_universe|={len(trans_universe)}"
    )

    all_results: list[dict] = []
    for code, (display, raw_id) in DISEASE_MAPPING.items():
        resolved = resolve_disease_id(raw_id, ontology)
        if resolved is None:
            warnings.warn(
                f"[{code}] disease ID '{raw_id}' not found in ontology (direct or "
                f"dbXRefs); skipping. Look up the OT canonical ID on the platform "
                f"or via a name search and update DISEASE_MAPPING."
            )
            continue
        gene_sets = load_gene_sets(gene_sets_dir, code)
        if gene_sets is None:
            warnings.warn(
                f"[{code}] no gene-set files matching '{code}_*' in {gene_sets_dir}; skipping"
            )
            continue

        subtree = expand_subtree(resolved, ontology)
        approved = approved_targets(drug_targets, subtree)
        print(
            f"[{code}] {display}  id={resolved}  "
            f"|subtree|={len(subtree)}  |approved targets|={len(approved)}"
        )

        results = run_tests(gene_sets, approved, cis_universe, trans_universe)
        for r in results:
            r["disease"] = code
            r["display_name"] = display
            r["disease_id"] = resolved
            r["is_imd"] = code in IMD_DISEASES
        all_results.extend(results)

    if not all_results:
        raise SystemExit(
            "No diseases were successfully processed. Fix mappings and re-run."
        )

    df = pl.DataFrame(all_results)
    p_vals = df.get_column("p_value").to_numpy()
    _, fdr, _, _ = multipletests(p_vals, method="fdr_bh")
    gene_set_rename = {
        "H4_genes": "coloc genes",
        "regulators_H4": "Regulators",
        "targets_H4": "Targets",
    }
    df = (
        df.with_columns(pl.Series("fdr_bh", fdr))
        .sort(
            ["is_imd", "disease", "gene_set", "condition"],
            descending=[True, False, False, False],
        )
        .with_columns(pl.col("gene_set").replace(gene_set_rename).alias("gene_set"))
    )

    out_csv = out_dir / "enrichment_all.csv"
    df.write_csv(out_csv)
    print(
        f"\nWrote {out_csv}  ({df.height} rows across {df.get_column('disease').n_unique()} diseases)"
    )

    out_png = out_dir / "enrichment_all.png"
    plot_heatmap(df, out_png)
    print(f"Wrote {out_png}")


if __name__ == "__main__":
    main()
