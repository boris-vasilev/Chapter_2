#!/usr/bin/env python
"""Disease-specificity test via a gene-set × drug-indication matrix.

RATIONALE
---------
A property-matched permutation that controls for CD4+ T-cell expression may be
too aggressive: disease-relevant signal lives in CD4+ identity itself (e.g.
IL5 is both a Th2 identity marker and the target of asthma drugs). A null
result under that design means "not separable from cell-type identity", not
"no signal".

This test uses an *internal* comparator. For each disease i (gene-set side),
we ask: is disease i's gene set enriched for disease i's own approved drug
targets MORE THAN for some other disease j's drug targets? The cell-type-
identity confound predicts uniform enrichment across the row (any Th2-rich
gene set will hit any IMD's drugs roughly equally). Disease-specific signal
predicts an elevated DIAGONAL relative to OFF-DIAGONAL.

Because disease serves as its own comparator, no external property-matched
null is required — the shared cell-type-identity component cancels in the
diagonal-vs-off-diagonal contrast.

WHAT THE SCRIPT PRODUCES
------------------------
Per (gene_set, condition):
  - A square-ish matrix of log2(OR), Fisher p, a/b/c/d cells, one per
    (gene_set_i, drug_indication_j) pair. SAME Fisher/universe logic as
    approved_drug_enrichment_all.py for direct comparability.
  - A diagonal-vs-off-diagonal contrast with a label-shuffle empirical p:
    permute the mapping between gene sets and their drug-indication labels,
    recompute the mean diagonal log2(OR), and ask where the observed value
    falls in that null distribution.
  - A heatmap with diagonal cells outlined and IMD/non-IMD blocks grouped.

Plus once across the run:
  - A Jaccard-overlap heatmap of the disease drug-target sets (so the reader
    can see how separable the columns even are; overlapping drug sets blur
    the diagonal).
  - A combined summary CSV with the empirical p per (gene_set, condition).

NOTE on asthma:
  AST / ASTAO / ASTCO all map to MONDO_0004979 (asthma) on the ontology side,
  so they share an identical drug-indication column. We dedupe columns by
  resolved disease ID (one asthma column, three asthma rows), and treat every
  one of those three rows' asthma-column cells as diagonal. Each row therefore
  has exactly one diagonal column.
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

# --------------------------------------------------------------------------
# Constants shared with approved_drug_enrichment_all.py (kept in sync)
# --------------------------------------------------------------------------
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

CONDITIONS = ["Rest", "Stim8hr", "Stim48hr"]
# Gene-set display name -> file-name token used by load_gene_sets()
GENE_SET_FILES: dict[str, str] = {
    "coloc genes": "H4_genes",
    "Regulators": "regulators_H4",
    "Targets": "targets_H4",
}


# --------------------------------------------------------------------------
# Helpers reused from approved_drug_enrichment_all.py
# --------------------------------------------------------------------------
def resolve_disease_id(query_id: str, ontology: pl.DataFrame) -> str | None:
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


def load_gene_set_for_study(
    gene_sets_dir: Path, study: str, gene_set_token: str, condition: str
) -> set[str]:
    """Pull genes for one (study, gene_set, condition) triple."""
    genes: set[str] = set()
    pattern = f"{study}_*_{condition}_{gene_set_token}.txt"
    for file in gene_sets_dir.glob(pattern):
        with open(file) as fh:
            for line in fh:
                g = line.strip()
                if g:
                    genes.add(g)
    return genes


def universe_for(
    gene_set_token: str, cis_universe: set[str], trans_universe: set[str]
) -> set[str]:
    # SAME universe logic as the main script: regulators are perturbation-side
    # (cis), everything else (targets, coloc genes) is detection-side (trans).
    if gene_set_token.startswith("regulators"):
        return cis_universe
    return trans_universe


def fisher_cell(gene_set: set[str], approved: set[str], universe: set[str]) -> dict:
    """Same Fisher test as approved_drug_enrichment_all.py:run_tests()."""
    gs = gene_set & universe
    approved_in_universe = approved & universe
    n_universe = len(universe)
    n_approved = len(approved_in_universe)
    a = len(gs & approved_in_universe)
    b = len(gs) - a
    c = n_approved - a
    d = n_universe - len(gs) - c
    odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
    # Cap inf for downstream log2; treat inf as 16x (same as main heatmap)
    or_for_log = 16.0 if np.isinf(odds_ratio) else float(odds_ratio)
    log2_or = float(np.log2(np.clip(or_for_log, 0.25, 16)))
    return {
        "odds_ratio": float(odds_ratio),
        "log2_or": log2_or,
        "p_value": float(p_value),
        "a": int(a),
        "b": int(b),
        "c": int(c),
        "d": int(d),
        "n_universe": int(n_universe),
        "n_approved": int(n_approved),
    }


# --------------------------------------------------------------------------
# New: diagonal-shuffle null
# --------------------------------------------------------------------------
def diagonal_shuffle_test(
    matrix: np.ndarray,
    row_disease_ids: list[str],
    col_disease_ids: list[str],
    n_perm: int,
    rng: np.random.Generator,
) -> dict:
    """Null = permute the mapping between rows (gene sets) and their disease
    labels, then recompute the mean log2(OR) across the *new* diagonal cells.

    Under H0 (no disease-specific signal), the mean of the cells lit by any
    random row-to-column relabelling should equal the observed mean.
    """
    col_idx = {d: i for i, d in enumerate(col_disease_ids)}

    def diag_mean(row_labels: list[str]) -> float:
        vals = []
        for i, d in enumerate(row_labels):
            j = col_idx.get(d)
            if j is None:
                continue
            v = matrix[i, j]
            if not np.isnan(v):
                vals.append(v)
        return float(np.mean(vals)) if vals else float("nan")

    observed = diag_mean(row_disease_ids)
    null = np.empty(n_perm, dtype=float)
    row_labels_arr = np.array(row_disease_ids)
    for k in range(n_perm):
        shuffled = list(rng.permutation(row_labels_arr))
        null[k] = diag_mean(shuffled)

    if np.isnan(observed):
        emp_p = float("nan")
    else:
        # (#null >= observed + 1) / (N + 1) — standard conservative estimator
        finite_null = null[~np.isnan(null)]
        emp_p = float((np.sum(finite_null >= observed) + 1) / (len(finite_null) + 1))

    # Off-diagonal mean as a second sanity stat (not used for the null itself).
    off_diag_vals = []
    for i, d in enumerate(row_disease_ids):
        for j, cd in enumerate(col_disease_ids):
            if cd == d:
                continue
            v = matrix[i, j]
            if not np.isnan(v):
                off_diag_vals.append(v)
    off_diag_mean = float(np.mean(off_diag_vals)) if off_diag_vals else float("nan")

    return {
        "observed_diag_mean_log2or": observed,
        "off_diag_mean_log2or": off_diag_mean,
        "null_diag_mean_log2or": null,
        "empirical_p": emp_p,
        "n_perm": n_perm,
    }


# --------------------------------------------------------------------------
# Plots
# --------------------------------------------------------------------------
def plot_matrix_heatmap(
    log2_or: np.ndarray,
    p_values: np.ndarray,
    row_labels: list[str],
    col_labels: list[str],
    row_is_imd: list[bool],
    col_is_imd: list[bool],
    diagonal_mask: np.ndarray,
    title: str,
    out_path: Path,
) -> None:
    """Heatmap: rows = gene-set studies, cols = unique disease drug indications.
    Same visual grammar as enrichment_all.png. Diagonal cells outlined in gold.
    """
    n_rows, n_cols = log2_or.shape
    fig_w = max(12, 0.85 * n_cols + 4)
    fig_h = max(9, 0.55 * n_rows + 4)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(log2_or, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")

    stroke = [
        path_effects.Stroke(linewidth=2.5, foreground="black"),
        path_effects.Normal(),
    ]
    for i in range(n_rows):
        for j in range(n_cols):
            v = log2_or[i, j]
            p = p_values[i, j]
            if np.isnan(v):
                ax.text(j, i, "·", ha="center", va="center", color="grey", fontsize=14)
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
                    fontsize=14,
                    fontweight="bold",
                )
                txt.set_path_effects(stroke)

    # Diagonal cells outlined gold (same convention as FDR-significant cells
    # in the main heatmap — different *meaning*, same visual cue: "this is
    # the cell to look at").
    for i in range(n_rows):
        for j in range(n_cols):
            if diagonal_mask[i, j]:
                ax.add_patch(
                    Rectangle(
                        (j - 0.5, i - 0.5),
                        1.0,
                        1.0,
                        fill=False,
                        edgecolor="#FFD700",
                        linewidth=2.5,
                    )
                )

    # IMD / non-IMD dividers on both axes.
    n_imd_rows = sum(row_is_imd)
    n_imd_cols = sum(col_is_imd)
    if 0 < n_imd_rows < n_rows:
        ax.axhline(n_imd_rows - 0.5, color="black", linewidth=1.5)
    if 0 < n_imd_cols < n_cols:
        ax.axvline(n_imd_cols - 0.5, color="black", linewidth=1.5)

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=12)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_labels, fontsize=12)
    ax.set_xlabel("Drug-indication disease (j)", fontsize=14)
    ax.set_ylabel("Gene-set disease (i)", fontsize=14)

    cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("log2(OR), capped at ±4", fontsize=13)
    cbar.ax.tick_params(labelsize=10)

    ax.set_title(title, fontsize=14, pad=14)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_jaccard_heatmap(
    j_matrix: np.ndarray,
    labels: list[str],
    is_imd: list[bool],
    out_path: Path,
) -> None:
    n = j_matrix.shape[0]
    fig_w = max(9, 0.85 * n + 3)
    fig, ax = plt.subplots(figsize=(fig_w, fig_w))
    im = ax.imshow(j_matrix, cmap="Reds", vmin=0, vmax=1, aspect="auto")
    for i in range(n):
        for j in range(n):
            v = j_matrix[i, j]
            if not np.isnan(v):
                ax.text(
                    j,
                    i,
                    f"{v:.2f}",
                    ha="center",
                    va="center",
                    color="black" if v < 0.5 else "white",
                    fontsize=9,
                )
    n_imd = sum(is_imd)
    if 0 < n_imd < n:
        ax.axhline(n_imd - 0.5, color="black", linewidth=1.5)
        ax.axvline(n_imd - 0.5, color="black", linewidth=1.5)
    ax.set_xticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=11)
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=11)
    cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Jaccard overlap of approved-target sets", fontsize=12)
    ax.set_title(
        "Drug-target set overlap across diseases\n"
        "(high values blur the diagonal contrast)",
        fontsize=13,
        pad=12,
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_diag_vs_offdiag(
    summary_rows: list[dict],
    out_path: Path,
) -> None:
    """One-stop visual: per (gene_set, condition), bars for diag vs off-diag
    means plus the label-shuffle null distribution."""
    n = len(summary_rows)
    fig, ax = plt.subplots(figsize=(max(8, 1.2 * n + 3), 5.5))

    xs = np.arange(n)
    diag = [r["observed_diag_mean_log2or"] for r in summary_rows]
    off = [r["off_diag_mean_log2or"] for r in summary_rows]

    ax.bar(
        xs - 0.18,
        diag,
        width=0.34,
        label="Diagonal mean log2(OR)",
        color="#FFD700",
        edgecolor="black",
        linewidth=1.0,
    )
    ax.bar(
        xs + 0.18,
        off,
        width=0.34,
        label="Off-diagonal mean log2(OR)",
        color="#888888",
        edgecolor="black",
        linewidth=1.0,
    )
    for i, r in enumerate(summary_rows):
        p = r["empirical_p"]
        if np.isnan(p):
            continue
        stars = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        if stars:
            y = max(r["observed_diag_mean_log2or"], r["off_diag_mean_log2or"])
            ax.text(i, y + 0.05, stars, ha="center", fontsize=14, fontweight="bold")

    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_xticks(xs)
    ax.set_xticklabels(
        [f"{r['gene_set']}\n{r['condition']}" for r in summary_rows],
        fontsize=10,
        rotation=0,
    )
    ax.set_ylabel("Mean log2(OR)", fontsize=12)
    ax.set_title(
        "Diagonal (own-disease) vs off-diagonal (other-disease) enrichment\n"
        "label-shuffle empirical p: * < 0.05, ** < 0.01, *** < 0.001",
        fontsize=13,
    )
    ax.legend(fontsize=11, loc="upper left", frameon=True)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-sets-dir", required=True)
    parser.add_argument("--drug-targets", required=True)
    parser.add_argument("--disease-ontology", required=True)
    parser.add_argument("--de-results", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument(
        "--n-perm",
        type=int,
        default=10000,
        help="Label-shuffle iterations for the diagonal-vs-off-diagonal null.",
    )
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    gene_sets_dir = Path(args.gene_sets_dir)
    rng = np.random.default_rng(args.seed)

    print("Loading ontology, drug-targets, Perturb-seq DE results...")
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
    print(f"  |cis|={len(cis_universe)} |trans|={len(trans_universe)}")

    # ----- Resolve disease IDs and pull approved-target sets per study -----
    resolved_per_study: dict[str, str] = {}
    for code, (display, raw_id) in DISEASE_MAPPING.items():
        r = resolve_disease_id(raw_id, ontology)
        if r is None:
            warnings.warn(f"[{code}] cannot resolve {raw_id}; skipping")
            continue
        resolved_per_study[code] = r

    approved_per_study: dict[str, set[str]] = {}
    for code, resolved in resolved_per_study.items():
        subtree = expand_subtree(resolved, ontology)
        approved_per_study[code] = approved_targets(drug_targets, subtree)

    # ----- Row axis: study codes (one row per gene-set construction) -------
    row_codes = [c for c in IMD_DISEASES + NON_IMD if c in resolved_per_study]
    row_disease_ids = [resolved_per_study[c] for c in row_codes]
    row_is_imd = [c in IMD_DISEASES for c in row_codes]
    row_labels = [SHORT_NAMES.get(c, c) for c in row_codes]

    # ----- Column axis: unique disease IDs (drug-indication sets) ----------
    # Asthma collapses to one column; everyone else stays distinct. Use the
    # first study code with that disease ID for the display label.
    col_disease_ids: list[str] = []
    col_codes: list[str] = []
    seen_ids: set[str] = set()
    for c in IMD_DISEASES + NON_IMD:
        if c not in resolved_per_study:
            continue
        did = resolved_per_study[c]
        if did in seen_ids:
            continue
        seen_ids.add(did)
        col_disease_ids.append(did)
        col_codes.append(c)
    col_labels = [SHORT_NAMES.get(c, c) for c in col_codes]
    col_is_imd = [c in IMD_DISEASES for c in col_codes]
    col_approved = [approved_per_study[c] for c in col_codes]
    print(
        f"\nRows (gene-set studies): {len(row_codes)}  "
        f"Cols (drug-indication diseases): {len(col_codes)}"
    )

    # ----- Jaccard overlap of approved-target sets across columns ----------
    n_col = len(col_codes)
    jaccard_matrix = np.full((n_col, n_col), np.nan)
    for i in range(n_col):
        for j in range(n_col):
            a, b = col_approved[i], col_approved[j]
            denom = len(a | b)
            jaccard_matrix[i, j] = (len(a & b) / denom) if denom else float("nan")
    plot_jaccard_heatmap(
        jaccard_matrix, col_labels, col_is_imd, out_dir / "drug_target_jaccard.png"
    )
    pl.DataFrame(
        {
            "disease_i": [col_codes[i] for i in range(n_col) for _ in range(n_col)],
            "disease_j": [col_codes[j] for _ in range(n_col) for j in range(n_col)],
            "jaccard": jaccard_matrix.flatten().tolist(),
        }
    ).write_csv(out_dir / "drug_target_jaccard.csv")

    # ----- Per (gene_set, condition): build the matrix and shuffle test ----
    all_long_rows: list[dict] = []
    summary_rows: list[dict] = []

    for gs_display, gs_token in GENE_SET_FILES.items():
        for cond in CONDITIONS:
            # H4 (coloc) is condition-invariant; we keep all three rows for
            # consistency with the main heatmap, but the gene set will be
            # identical across Rest/Stim8hr/Stim48hr.
            print(f"\n[{gs_display} | {cond}]")
            n_row = len(row_codes)
            log2_or = np.full((n_row, n_col), np.nan)
            p_values = np.full((n_row, n_col), np.nan)
            a_matrix = np.full((n_row, n_col), np.nan)
            universe = universe_for(gs_token, cis_universe, trans_universe)

            for i, code in enumerate(row_codes):
                gs = load_gene_set_for_study(gene_sets_dir, code, gs_token, cond)
                if not gs:
                    # Missing gene-set file; row stays NaN
                    continue
                for j in range(n_col):
                    res = fisher_cell(gs, col_approved[j], universe)
                    log2_or[i, j] = res["log2_or"]
                    p_values[i, j] = res["p_value"]
                    a_matrix[i, j] = res["a"]
                    all_long_rows.append(
                        {
                            "gene_set": gs_display,
                            "condition": cond,
                            "row_code": code,
                            "row_disease_id": resolved_per_study[code],
                            "col_code": col_codes[j],
                            "col_disease_id": col_disease_ids[j],
                            "is_diagonal": resolved_per_study[code]
                            == col_disease_ids[j],
                            "row_is_imd": code in IMD_DISEASES,
                            "col_is_imd": col_codes[j] in IMD_DISEASES,
                            "odds_ratio": res["odds_ratio"],
                            "log2_or": res["log2_or"],
                            "p_value": res["p_value"],
                            "a": res["a"],
                            "b": res["b"],
                            "c": res["c"],
                            "d": res["d"],
                            "n_universe": res["n_universe"],
                            "n_approved": res["n_approved"],
                        }
                    )

            # Diagonal mask for the heatmap (gold-bordered cells).
            diagonal_mask = np.zeros((n_row, n_col), dtype=bool)
            for i, did in enumerate(row_disease_ids):
                if did in col_disease_ids:
                    diagonal_mask[i, col_disease_ids.index(did)] = True

            # Run the label-shuffle null on this (gene_set, condition).
            shuffle_res = diagonal_shuffle_test(
                log2_or, row_disease_ids, col_disease_ids, args.n_perm, rng
            )
            summary_rows.append(
                {
                    "gene_set": gs_display,
                    "condition": cond,
                    "n_rows": n_row,
                    "n_cols": n_col,
                    "n_perm": args.n_perm,
                    "observed_diag_mean_log2or": shuffle_res[
                        "observed_diag_mean_log2or"
                    ],
                    "off_diag_mean_log2or": shuffle_res["off_diag_mean_log2or"],
                    "null_diag_mean_log2or_median": float(
                        np.nanmedian(shuffle_res["null_diag_mean_log2or"])
                    ),
                    "null_diag_mean_log2or_p95": float(
                        np.nanpercentile(shuffle_res["null_diag_mean_log2or"], 95)
                    ),
                    "empirical_p": shuffle_res["empirical_p"],
                }
            )
            print(
                f"  diag mean log2(OR) = {shuffle_res['observed_diag_mean_log2or']:.3f}  "
                f"off-diag = {shuffle_res['off_diag_mean_log2or']:.3f}  "
                f"emp_p = {shuffle_res['empirical_p']:.4g}"
            )

            # Heatmap per (gene_set, condition).
            title = (
                f"{gs_display} | {cond}  —  disease × drug-indication enrichment\n"
                f"diag mean log2(OR) = {shuffle_res['observed_diag_mean_log2or']:.3f}, "
                f"off-diag = {shuffle_res['off_diag_mean_log2or']:.3f}, "
                f"label-shuffle emp_p = {shuffle_res['empirical_p']:.4g} "
                f"({args.n_perm} perms)"
            )
            safe_gs = gs_display.replace(" ", "_")
            plot_matrix_heatmap(
                log2_or,
                p_values,
                row_labels,
                col_labels,
                row_is_imd,
                col_is_imd,
                diagonal_mask,
                title,
                out_dir / f"matrix_{safe_gs}_{cond}.png",
            )

    # ----- Write long + summary CSVs --------------------------------------
    long_df = pl.DataFrame(all_long_rows)
    long_df.write_csv(out_dir / "specificity_matrix_long.csv")
    print(
        f"\nWrote {out_dir / 'specificity_matrix_long.csv'} "
        f"({long_df.height} cells across {len(summary_rows)} matrices)"
    )

    summary_df = pl.DataFrame(summary_rows)
    summary_df.write_csv(out_dir / "specificity_matrix_summary.csv")
    print(f"Wrote {out_dir / 'specificity_matrix_summary.csv'}")

    plot_diag_vs_offdiag(summary_rows, out_dir / "diag_vs_offdiag.png")
    print(f"Wrote {out_dir / 'diag_vs_offdiag.png'}")


if __name__ == "__main__":
    main()
