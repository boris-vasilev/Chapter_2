#!/usr/bin/env python
"""Rebuild figures/hub_exclusion/ from scratch.

    cd vis/hub_exclusion && uv run 00_build.py

1. Empties figures/hub_exclusion/.
2. Runs 01_collect.py (tables + Check 1, fails on any gene-set mismatch),
   then the figure scripts and 07_numbers_for_text.py.
3. Writes figures/hub_exclusion/INDEX.md.
4. Fails if any file in the folder is not in the index, or vice versa.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import hxcommon  # noqa: F401
from common import CONDS
from hxcommon import BASE, HERE, OUT

SCRIPTS = ["01_collect.py", "01b_degree_null.py", "02_fig_h1_overlap.py",
           "03_fig_s_overlap_heatmaps.py", "03b_fig_s_h1_targets_summary.py",
           "04_fig_h2_heritability.py", "05_fig_s_h2_forest.py", "06_fig_s_h2_delta.py",
           "06b_fig_s_h2_targets_summary.py", "07_numbers_for_text.py"]
VARIANTS = f"{BASE} and hub0.99 / hub0.95 / hub0.90_{BASE}"

FIGS = [  # (label, stem, message, script)
    ("Fig. H1", "main/fig_h1_overlap",
     "Regulator-set overlap vs hub exclusion. IMD–IMD excess stays above both the "
     "size-only and the degree-aware expectation at every threshold; IMD–control and "
     "control–control fall to zero. A: Stim 8 h observed and size-only excess Jaccard per "
     "variant; B: median (IQR) size-only excess by pair class; E: same against the "
     "degree-aware expectation (filled = empirical FDR < 0.05); C: Spearman ρ with 0 %; "
     "D: median regulator-set size.", "02_fig_h1_overlap.py"),
    ("Fig. H2", "main/fig_h2_heritability",
     "IMD regulator-set enrichment and τ* fall with hub exclusion; control sets carried "
     "almost no signal at 0 %. No τ* change reaches FDR < 0.05 under the conservative "
     "test. A/C: Stim 8 h regulator estimates and change vs 0 %; B: IMD vs control "
     "summaries (pooled); D: enrichment SE and log2-enrichment CI width vs annotation "
     "size.", "04_fig_h2_heritability.py"),
]
for annot in ["regulators", "targets"]:
    for c in CONDS:
        if (annot, c) == ("regulators", "Stim8hr"):
            continue
        FIGS.append((f"Fig. S-H1 ({annot}, {c})", f"supp/fig_s_h1_overlap_{annot}_{c}",
                     f"As Fig. H1 A for {annot}, {c}. Open circle on the diagonal = target "
                     "set lost ≥ 1 seeding anchor as a hub.", "03_fig_s_overlap_heatmaps.py"))
FIGS.append(("Fig. S-H1 targets", "supp/fig_s_h1_targets_summary",
             "Target-set excess Jaccard by pair class (all sets; anchors unchanged) and "
             "Spearman ρ vs 0 %. Targets change only through anchor loss.",
             "03b_fig_s_h1_targets_summary.py"))
FIGS.append(("Fig. S-H2 targets", "supp/fig_s_h2_targets_summary",
             "Target-set IMD vs control heritability summaries (all sets; anchors "
             "unchanged) and precision vs annotation size.",
             "06b_fig_s_h2_targets_summary.py"))
for c in CONDS:
    FIGS.append((f"Fig. S-H2 forest ({c})", f"supp/fig_s_h2_forest_{c}",
                 f"Log2 enrichment and τ* with 95 % CI for anchors, regulators and targets, "
                 f"{c}, per variant. Filled = FDR < 0.05; diamond = target set lost an "
                 "anchor; grey = no estimate.", "05_fig_s_h2_forest.py"))
for c in CONDS:
    FIGS.append((f"Fig. S-H2 change ({c})", f"supp/fig_s_h2_delta_{c}",
                 f"Δ log2 enrichment and Δτ* vs 0 % for regulators and targets, {c}. CI "
                 "conservative: ignores positive covariance between nested annotations.",
                 "06_fig_s_h2_delta.py"))

DATA = [  # (file, message, script)
    ("hub_regulators.csv", "Hub regulators per threshold x condition (n_degs, anchor status).", "01_collect.py"),
    ("gene_set_rebuild.csv", "Check 1: every gene set rebuilt from the hub-filtered graph vs gene_sets/.", "01_collect.py"),
    ("set_sizes.csv", "Set sizes, universes, anchors lost, unchanged-anchor flag.", "01_collect.py"),
    ("jaccard_pairs.csv", "Pairwise observed / expected / excess Jaccard (Fig. H1 A, S-H1).", "01_collect.py"),
    ("heritability_long.csv", "S-LDSC estimates for all variants (01_collect_results parsers).", "01_collect.py"),
    ("anchor_shifts.csv", "Anchor estimates vs 0 % (tolerance check).", "01_collect.py"),
    ("top_shared_regulators.csv", "Top 20 regulators by number of IMD sets, per variant x condition.", "01_collect.py"),
    ("hub_anchor_fraction.csv", "Fraction of perturbed anchors that are hubs vs all perturbations.", "01_collect.py"),
    ("jaccard_degree_null_pairs.csv", "Pairwise degree-aware expected and excess Jaccard, regulators.", "01b_degree_null.py"),
    ("h1e_degree_null_summary.csv", "Fig. H1 E: degree-aware median excess, empirical p and FDR by pair class.", "01b_degree_null.py"),
    ("h1b_excess_summary.csv", "Fig. H1 B, S-H1 targets: median / IQR size-only excess Jaccard by pair class.", "02_fig_h1_overlap.py"),
    ("h1c_spearman.csv", "Fig. H1 C: Spearman ρ vs 0 %.", "02_fig_h1_overlap.py"),
    ("h1d_set_size.csv", "Fig. H1 D: median regulator-set size.", "02_fig_h1_overlap.py"),
    ("h2a_estimates.csv", "Fig. H2 A, S-H2 forest: estimates and FDR calls.", "04_fig_h2_heritability.py"),
    ("h2b_summary.csv", "Fig. H2 B: IMD vs control medians and significant counts.", "04_fig_h2_heritability.py"),
    ("h2c_delta.csv", "Fig. H2 C, S-H2 change: Δ vs 0 % with conservative z test.", "04_fig_h2_heritability.py"),
    ("h2d_precision.csv", "Fig. H2 D, S-H2 targets: enrichment SE and CI widths vs proportion of SNPs.", "04_fig_h2_heritability.py"),
]


def rows():
    out = []
    for label, stem, msg, script in FIGS:
        for ext in ("pdf", "png"):
            out.append((label, f"{stem}.{ext}", msg, script))
    out += [("Data", f"data/{f}", m, s) for f, m, s in DATA]
    out += [("Text", "numbers_for_text.md", "Every number quoted in the text, with the "
             "closing statement on Q1 / Q2.", "07_numbers_for_text.py"),
            ("—", "INDEX.md", "This index.", "00_build.py")]
    return out


def write_index(rs):
    lines = ["# Hub exclusion — index", "",
             f"Variants: **{VARIANTS}** (0, 1, 5, 10 % of regulators excluded per condition "
             "by out-degree). Rebuild with `uv run 00_build.py` in `vis/hub_exclusion/`.", "",
             "| Figure | File | Message | Command (from vis/hub_exclusion) |",
             "|---|---|---|---|"]
    lines += [f"| {lab} | `{f}` | {m.replace('|', '/')} | `uv run {s}` |" for lab, f, m, s in rs]
    (OUT / "INDEX.md").write_text("\n".join(lines) + "\n")


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    for s in SCRIPTS:
        print(f"== {s}", flush=True)
        subprocess.run([sys.executable, str(HERE / s)], check=True, cwd=HERE)
    rs = rows()
    write_index(rs)
    indexed = {f for _, f, _, _ in rs}
    present = {str(p.relative_to(OUT)) for p in OUT.rglob("*") if p.is_file()}
    if indexed != present:
        raise SystemExit(f"Index mismatch. Not indexed: {sorted(present - indexed)}; "
                         f"missing: {sorted(indexed - present)}")
    print(f"Built {len(present)} files in {OUT}; INDEX.md complete.")


if __name__ == "__main__":
    main()
