#!/usr/bin/env python
"""Rebuild every Chapter 2 figure, table and data file from scratch.

    uv run 00_build_all.py [--schematic PATH]   # default: assets/ch2_schematic.png
    CH2_VARIANT=<results/ldsc folder name> uv run 00_build_all.py

The design is the single setting VARIANT in common.py (env var CH2_VARIANT);
inputs come from LDSC_DIR / PERM_DIR and outputs go to FIG_DIR.

1. Empties FIG_DIR and the staging folder build/<VARIANT>/.
2. Runs the figure scripts (each writes to staging; the permutation-null
   plots and the matched-null scripts in permutation/ write straight to
   FIG_DIR/permutation/, skipped if PERM_DIR is absent; fig_s_matching_check
   is skipped if its inputs outside results/ are missing). The deg0 vs deg75
   comparison (permutation/*deg_filter*) always runs and writes to
   figures2/comparison/, outside FIG_DIR.
3. Copies the chapter outputs into FIG_DIR/{main,supp,data}/ under their
   figure numbers (MANIFEST below), plus NOTES_overview.md from notes/.
4. Writes FIG_DIR/INDEX.md and FIG_DIR/contact_sheet.pdf.
5. Fails if any file in FIG_DIR is not in the index, or vice versa.
"""

import argparse
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

from common import BUILD_DIR, FIG_DIR, HERE, PERM_DIR, VARIANT

STAGE = BUILD_DIR
PERM = "plot_permutation_null.py"   # in heritability/; writes to FIG_DIR/permutation
P = "permutation"   # matched-null scripts; write to FIG_DIR/permutation (no --out-dir)
# script -> output file stems in FIG_DIR/permutation (for the index)
P_OUTPUTS = {
    "fig_5_matched_null.py": ["fig_5_matched_null"],
    "fig_s_null_histograms.py": ["fig_s_null_histograms"],
    "table_s_matched_null.py": ["table_s_matched_null"],
    "fig_s_matching_check.py": ["fig_s_matching_check"],
}
P_COMPARISON = ["fig_s_deg_filter.py", "table_s_deg_filter.py"]   # -> figures2/comparison
sys.path.insert(0, str(HERE / P))
from pcommon import MATCHING_INPUTS  # noqa: E402
NOTES_SRC = HERE.parent / "notes" / "NOTES_overview.md"
NOTES_H_SRC = HERE.parent / "notes" / "NOTES_heritability.md"
H = "heritability"   # subfolder for the heritability scripts, staging and outputs

RUNS = [
    [f"{H}/01_collect_results.py"],   # first: Fig. 2.2C reads its table
    ["01_coloc_summary_table.py"],
    ["02_gene_set_sizes_barplot.py"],
    ["03_condition_size_change.py"],
    ["04_jaccard_heatmaps.py"],
    ["04_jaccard_heatmaps.py", "--excess"],
    ["04_jaccard_heatmaps.py", "--condition", "Stim8hr", "--pair"],
    ["05_perturbation_graph.py"],
    ["06_anchor_coverage.py"],
    ["07_fig_2_1.py"],
    [f"{H}/02_fig_4_heritability.py"],
    [f"{H}/03_fig_s_forest.py"],
    [f"{H}/04_fig_s_scatter.py"],
    [f"{H}/05_numbers_for_text.py"],
    [f"{H}/{PERM}"],
    *[[f"{P}/{f}"] for f in P_OUTPUTS],
    *[[f"{P}/{f}"] for f in P_COMPARISON],
]

SEC_TRAITS = "Traits and anchors"
SEC_DESIGN = "Study design and anchors"
SEC_GRAPH = "Perturbation-response graph and expanded sets"
SEC_OVERLAP = "Overlap of expanded sets between diseases"
SEC_COND = "Expanded sets across stimulation conditions"
SEC_HERIT = "Heritability enrichment of expanded sets"

# One entry per chapter item. kind: table | figure | data.
# src is the staging file name; for figures the PNG next to it feeds the
# contact sheet. cmd is the command that writes src (from chapter2/vis).
MANIFEST = [
    dict(num="Table 2.1", kind="table", section=SEC_TRAITS,
         files={"main/table_2_1_traits.tex": "table_coloc_summary.tex",
                "main/table_2_1_traits.csv": "table_coloc_summary.csv"},
         cmd="uv run 01_coloc_summary_table.py",
         message="Eighteen GWAS traits (12 immune-mediated, 6 non-IMD controls) "
                 "colocalise with 2–7 CD4+ eQTL studies and contribute 3–63 coloc "
                 "genes (anchors).",
         caption=("Complete LaTeX table float (caption + \\label{tab:traits}); \\input it "
                 "and cite it with Table~\\ref{tab:traits}. Needs booktabs and tabularx. "
                 "Columns: disease, ancestry, GWAS (\\cite), cases, controls, number of "
                 "CD4+ T-cell eQTL studies in which the GWAS signal colocalised with an "
                 "eQTL (PP.H4 > 0.8), number of coloc genes. Immune-mediated diseases "
                 "first, then non-IMD controls, each by coloc genes. Note under the "
                 "table: traits from the same GWAS share samples (de Lange 2017; "
                 "Ferreira 2019; Ishigaki 2022). Per-disease eQTL study names are in "
                 "data/table_2_1_eqtl_sources.csv (study mapping unverified).")),
    dict(num="Fig. 2.1", kind="figure", section=SEC_DESIGN,
         files={"main/fig_2_1_overview.pdf": "fig_2_1_overview.pdf"},
         cmd="uv run 07_fig_2_1.py",
         message="Most anchors respond to some knockdown but few are strong "
                 "knockdowns themselves, and anchor sets overlap mainly between "
                 "related immune-mediated diseases.",
         caption=('(A) Study design: GWAS and cis-eQTL colocalisation gives coloc genes '
                 '(anchors, PP.H4 > 0.8); Perturb-seq expands them to regulators (knockdowns '
                 'that change an anchor) and targets (genes changed when an anchor is '
                 'knocked down); each set is tested for heritability enrichment (S-LDSC) and '
                 'drug-target enrichment. (B) Anchors per disease and condition, classified '
                 'by their edges in the perturbation graph (|logFC| ≥ 0.5): both incoming '
                 'and outgoing (yellow), incoming only (green), outgoing only (orange), or '
                 'none (grey). * no targets in any condition (no anchor has an outgoing '
                 'edge). (C) Jaccard index of anchor sets between diseases, ordered by '
                 'average-linkage clustering; columns are numbered as the rows; diagonal '
                 'masked. Disease names are coloured by class (purple, immune-mediated; '
                 'grey, non-IMD control).')),
    dict(num="Fig. 2.2", kind="figure", section=SEC_GRAPH,
         files={"main/fig_2_2_graph_and_sets.pdf": "perturbation_graph_and_sets.pdf"},
         cmd="uv run 05_perturbation_graph.py",
         message="A few knockdowns affect many genes and some genes respond to many "
                 "knockdowns, so expanded sets are large and often come from a "
                 "single anchor.",
         caption=('(A, B) Complementary CDFs of out-degree (genes changed per knockdown; n = '
                 '894, 970 and 934 knockdowns at Rest, 8 h and 48 h) and in-degree '
                 '(knockdowns changing a gene; n = 9,544, 9,775 and 9,665 genes); the three '
                 'largest hubs are labelled with the condition in which they occur. (C) '
                 'Percentage of reference SNPs in the anchor, regulator and target '
                 'annotations: violins (density estimated on the log scale), one point '
                 'per disease, black bar at the median '
                 "(fitted S-LDSC sets, incl. the single-set runs for osteoarthritis "
                 "regulators and Alzheimer's proxy 48 h targets; no results for "
                 "schizophrenia, the Alzheimer's proxy at Rest and 8 h, or sets under "
                 '1,000 SNPs). (D) Regulator and target set sizes per disease and '
                 'condition; blue ticks, anchor set size; dashed lines, all knockdowns with '
                 "an effect. Empty target sets are shown in the '0' strip left of the axis "
                 'break. (E) Share of each expanded set contributed by its largest anchor; '
                 'dominant anchors are named where the target share exceeds 0.5. Empty sets '
                 'are omitted: Asthma (EAS), Ulcerative colitis, Osteoarthritis, Ischaemic '
                 "stroke (all conditions), Alzheimer's disease (proxy) (Rest, 8 h), Type 2 "
                 'diabetes (Rest).')),
    dict(num="Fig. 2.6A", kind="figure", section=SEC_OVERLAP,
         files={"main/fig_2_6a_regulator_overlap_8h.pdf":
                "jaccard_regulators_Stim8hr_pair.pdf"},
         cmd="uv run 04_jaccard_heatmaps.py --condition Stim8hr --pair",
         message="At 8 h, regulator sets overlap more than expected by chance for "
                 "every pair of diseases, most strongly among immune-mediated "
                 "diseases.",
         caption=('Observed Jaccard index (left) and observed minus expected Jaccard (right) '
                 'between regulator sets at 8 h of stimulation. The expectation is for '
                 'random sets of the same sizes drawn from all knockdowns with an effect; no '
                 'pair falls below expectation at 8 h. Rows and columns share one order '
                 '(clustering on mean observed Jaccard across conditions); columns are '
                 'numbered as the rows; diagonal masked. Black outlines mark pairs of traits '
                 'from the same GWAS, which share samples (de Lange 2017; Ferreira 2019; '
                 'Ishigaki 2022).')),
    dict(num="Fig. 4", kind="figure", section=SEC_HERIT,
         files={f"main/{H}/fig_4_heritability.pdf": f"{H}/fig_4_heritability.pdf"},
         cmd=f"uv run {H}/02_fig_4_heritability.py",
         message="Expanded sets are marginally enriched for immune-mediated disease "
                 "heritability, but not conditionally on the baseline and anchor "
                 "annotations.",
         caption=("S-LDSC for anchor, regulator and target sets (baselineLD v2.2 plus the "
                  "three sets fitted jointly per trait and condition). (A) Heritability "
                  "enrichment: colour, log2 enrichment (clipped at ±3); dot area, −log10 p "
                  "(two-sided). (B) τ*, the standardised per-annotation effect conditional "
                  "on all other annotations in the model (clipped at ±0.5); dot area, "
                  "−log10 one-sided p (H1: τ* > 0). Black outline, BH FDR < 0.05, corrected "
                  "within regulator + target sets (84 tests) and within anchors (45 tests) "
                  "separately, for A and B. Lupus regulators at Rest have negative "
                  "Prop._h2 (enrichment −11.8, drawn at the blue end) and are not counted "
                  "as significant. Grey square, empty set; grey square with ×, no estimate "
                  "(schizophrenia: S-LDSC failed; Alzheimer's proxy at Rest and 8 h: not run; "
                  "sets under 1,000 SNPs dropped before fitting). Rows in the order of "
                  "Fig. 2.2D. * multi-ancestry GWAS analysed with the European LD "
                  "reference; East Asian GWAS used the East Asian reference.")),
    dict(num="Fig. S2A", kind="figure", section=SEC_OVERLAP,
         files={"supp/fig_s2a_regulator_overlap_observed.pdf": "jaccard_regulators_all.pdf"},
         cmd="uv run 04_jaccard_heatmaps.py",
         message="Regulator-set overlap between diseases is similar in all three "
                 "conditions.",
         caption="Observed Jaccard index between regulator sets at Rest, 8 h and "
                 "48 h, shared ordering."),
    dict(num="Fig. S2B", kind="figure", section=SEC_OVERLAP,
         files={"supp/fig_s2b_regulator_overlap_excess.pdf":
                "jaccard_regulators_all_excess.pdf"},
         cmd="uv run 04_jaccard_heatmaps.py --excess",
         message="In every condition, regulator-set overlap exceeds the random "
                 "expectation for almost all disease pairs, most for pairs of "
                 "immune-mediated diseases.",
         caption="Observed minus expected Jaccard between regulator sets per "
                 "condition (universe: knockdowns with ≥1 edge in that condition)."),
    dict(num="Fig. S2C", kind="figure", section=SEC_OVERLAP,
         files={"supp/fig_s2c_target_overlap_observed.pdf": "jaccard_targets_all.pdf"},
         cmd="uv run 04_jaccard_heatmaps.py",
         message="Target sets overlap only between diseases that share a dominant "
                 "anchor.",
         caption="Observed Jaccard index between target sets at Rest, 8 h and 48 h, "
                 "shared ordering."),
    dict(num="Fig. S2D", kind="figure", section=SEC_OVERLAP,
         files={"supp/fig_s2d_target_overlap_excess.pdf": "jaccard_targets_all_excess.pdf"},
         cmd="uv run 04_jaccard_heatmaps.py --excess",
         message="Excess target-set overlap is near zero except for diseases that "
                 "share a dominant anchor.",
         caption="Observed minus expected Jaccard between target sets per condition "
                 "(universe: genes with ≥1 incoming edge). Grey, undefined (a "
                 "disease with no targets)."),
    dict(num="Fig. S3", kind="figure", section=SEC_COND,
         files={"supp/fig_s3_set_size_by_condition.pdf": "supp_set_size_by_condition.pdf"},
         cmd="uv run 03_condition_size_change.py",
         message="Regulator and target sets are largest at 8 h of stimulation.",
         caption="Set size per disease across Rest, 8 h and 48 h (top, log scale) "
                 "and log2 fold change versus Rest (bottom); black, median across "
                 "diseases. Diseases with an empty set in any condition are "
                 "excluded (listed under the figure). Friedman / Wilcoxon results "
                 "are in data/set_size_by_condition_tests.csv, not on the figure "
                 "(diseases are not independent)."),
    dict(num="Fig. S4", kind="figure", section=SEC_GRAPH,
         files={"supp/fig_s4_stacked_set_sizes.pdf": "supp_gene_set_sizes_by_condition.pdf"},
         cmd="uv run 02_gene_set_sizes_barplot.py",
         message="Expanded sets are one to two orders of magnitude larger than "
                 "anchor sets.",
         caption="Coloc genes (own axis, condition-independent) and stacked target "
                 "+ regulator set sizes per condition. Stacked length is not the "
                 "union size: regulators and targets can share genes."),
    dict(num="Fig. S5", kind="figure", section=SEC_HERIT,
         files={f"supp/{H}/fig_s_heritability_forest.pdf":
                f"{H}/fig_s_heritability_forest.pdf"},
         cmd=f"uv run {H}/03_fig_s_forest.py",
         message="Expanded-set enrichment intervals mostly exclude one, while their τ* "
                 "intervals mostly include zero.",
         caption="Enrichment (A) and τ* (B) with 95% CIs per trait, annotation and "
                 "condition (colour); filled points significant as in Fig. 4. Each panel's "
                 "x range covers the bulk of the estimates; CIs past it end in a small "
                 "arrowhead at the panel edge, and estimates past it are drawn as a larger "
                 "arrowhead. * as in Fig. 4."),
    dict(num="Fig. S6", kind="figure", section=SEC_HERIT,
         files={f"supp/{H}/fig_s_enrichment_vs_coverage.pdf":
                f"{H}/fig_s_enrichment_vs_coverage.pdf"},
         cmd=f"uv run {H}/04_fig_s_scatter.py",
         message="Enrichment falls as annotations get larger: small anchor sets "
                 "have the largest enrichment.",
         caption="log2 enrichment against % of reference SNPs in each fitted set; colour, "
                 "annotation; shape, disease class. Lupus regulators at Rest (enrichment "
                 "−11.8) are not shown (no log2)."),
    dict(num="Fig. S7", kind="figure", section=SEC_HERIT,
         files={f"supp/{H}/fig_s_enrichment_vs_tau.pdf": f"{H}/fig_s_enrichment_vs_tau.pdf"},
         cmd=f"uv run {H}/04_fig_s_scatter.py",
         message="No enriched regulator or target set passes the τ* threshold; "
                 "anchors pass both or neither.",
         caption="Signed enrichment z (from the S-LDSC enrichment p) against τ* z per "
                 "fitted set. Dashed lines, FDR 5% thresholds for regulator + target sets "
                 "(enrichment two-sided, |z| = 2.19; τ* one-sided, z = 3.24, the BH hurdle "
                 "for the smallest p since no set passes); dotted, nominal p = 0.05 "
                 "(|z| = 1.96; one-sided z = 1.64). Anchors are corrected separately "
                 "(lines not drawn). Bottom-left point: lupus regulators at Rest "
                 "(negative Prop._h2)."),
]

DATA = [
    ("Table 2.1", "table_coloc_eqtl_sources.csv", "table_2_1_eqtl_sources.csv",
     "uv run 01_coloc_summary_table.py",
     "Per-disease eQTL study list for the Table 2.1 footnote (unverified mapping)."),
    ("Fig. 2.1B", "anchor_coverage.csv", "anchor_coverage.csv",
     "uv run 06_anchor_coverage.py",
     "One row per (disease, condition, anchor) with its edge category."),
    ("Fig. 2.1C", "jaccard_anchors_observed.csv", "jaccard_anchors_observed.csv",
     "uv run 04_jaccard_heatmaps.py", "Anchor Jaccard matrix (clustered order)."),
    ("Fig. 2.2A, B", "degrees.csv", "degrees.csv",
     "uv run 05_perturbation_graph.py", "Out- and in-degree per gene and condition."),
    ("Fig. 2.2D, S3, S4", "set_sizes.csv", "set_sizes.csv",
     "uv run 05_perturbation_graph.py",
     "Anchor, regulator and target set sizes per disease and condition."),
    ("Fig. 2.2C", "genome_coverage.csv", "genome_coverage.csv",
     "uv run 05_perturbation_graph.py",
     "% reference SNPs per annotation (S-LDSC Prop._SNPs × 100)."),
    ("Fig. 2.2E", "largest_anchor_share.csv", "largest_anchor_share.csv",
     "uv run 05_perturbation_graph.py",
     "Largest single-anchor share of each non-empty expanded set."),
    ("Fig. 2.2 (text)", "overview_numbers.csv", "overview_numbers.csv",
     "uv run 05_perturbation_graph.py",
     "Per-condition universe sizes, degree summaries and median coverage."),
    ("Fig. S3", "supp_set_size_by_condition_tests.csv", "set_size_by_condition_tests.csv",
     "uv run 03_condition_size_change.py",
     "Friedman / Wilcoxon tests (not on the figure) and excluded diseases."),
]
for annot in ["regulators", "targets"]:
    for cond in ["Rest", "Stim8hr", "Stim48hr"]:
        for kind in ["observed", "expected"]:
            fig = ("Fig. 2.6A, S2A, S2B" if annot == "regulators" and cond == "Stim8hr"
                   else "Fig. S2A, S2B" if annot == "regulators" else "Fig. S2C, S2D")
            f = f"jaccard_{annot}_{cond}_{kind}.csv"
            DATA.append((fig, f, f, "uv run 04_jaccard_heatmaps.py",
                         f"{kind.capitalize()} Jaccard matrix, {annot}, {cond}."))

H_DATA = [
    ("Fig. 4, S5–S7", "heritability_long.csv",
     "One row per trait × condition × annotation: status, set size, Prop._SNPs, "
     "Prop._h2 (+ negative flag), FDR family, enrichment (SE, p, z, FDR, sig), "
     "τ (SE, z), τ* (SE, z, one-sided p, FDR, sig; two-sided p/FDR/sig kept)."),
    ("Fig. 4 (methods)", "ldsc_runs.csv",
     "One row per trait × condition: LD reference, files present, h2g, implied M, "
     "run status."),
]
for fig, f, desc in H_DATA:
    DATA.append((fig, f"{H}/{f}", f"{H}/{f}", f"uv run {H}/01_collect_results.py", desc))
DATA.append(("Text", f"{H}/numbers_for_text.md", "numbers_for_text.md",
             f"uv run {H}/05_numbers_for_text.py",
             "Numbers for the chapter text (Heritability section)."))

DOCS = [
    ("—", "INDEX.md", "This index."),
    ("—", "contact_sheet.pdf", "One page per figure and table with number, message "
                               "and caption, then the full file list."),
    ("—", "NOTES_overview.md", "Definitions, checks, measured results, "
                               "interpretation and open questions."),
    ("—", "NOTES_heritability.md", "S-LDSC inventory (Step 0), τ* computation and "
                                   "check, measured results, interpretation, not done."),
]


def run_scripts(schematic):
    for argv in RUNS:
        if argv[0] == "07_fig_2_1.py" and schematic:
            argv = argv + ["--schematic", schematic]
        script = Path(argv[0]).name
        needs_perm = argv[0] == f"{H}/{PERM}" or script in P_OUTPUTS
        if needs_perm and not PERM_DIR.is_dir():
            print(f"WARNING: skipping {argv[0]}: no permutation results at "
                  f"{PERM_DIR}", flush=True)
            continue
        missing = [str(f) for f in MATCHING_INPUTS if not f.exists()]
        if script == "fig_s_matching_check.py" and missing:
            print(f"WARNING: skipping {argv[0]}: inputs outside results/ missing: "
                  f"{', '.join(missing)}", flush=True)
            continue
        if argv[0] == f"{H}/{PERM}":
            out = ["--out-dir", str(FIG_DIR / "permutation")]
        elif argv[0].startswith(f"{P}/"):
            out = []   # paths fixed by pcommon / FIGURES2
        else:
            out = ["--out-dir", str(STAGE / H if argv[0].startswith(f"{H}/") else STAGE)]
        print(f"$ uv run {' '.join(argv)}", flush=True)
        subprocess.run([sys.executable, *argv, *out],
                       cwd=HERE, check=True, stdout=subprocess.DEVNULL)


def copy_outputs():
    for item in MANIFEST:
        for dst, src in item["files"].items():
            out = FIG_DIR / dst
            out.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(STAGE / src, out)
    for _, src, dst, _, _ in DATA:
        out = FIG_DIR / "data" / dst
        out.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(STAGE / src, out)
    shutil.copy2(NOTES_SRC, FIG_DIR / "NOTES_overview.md")
    shutil.copy2(NOTES_H_SRC, FIG_DIR / "NOTES_heritability.md")


def index_rows():
    rows = []
    for item in MANIFEST:
        for dst in item["files"]:
            rows.append((item["num"], dst, item["section"], item["message"],
                         item["cmd"], VARIANT))
    for fig, _, dst, cmd, desc in DATA:
        rows.append((fig, f"data/{dst}", "Data", desc, cmd, VARIANT))
    for num, f, desc in DOCS:
        rows.append((num, f, "—", desc, "uv run 00_build_all.py", VARIANT))
    for png in sorted((FIG_DIR / "permutation").glob("*/*.png")):
        rows.append(("Permutation null", str(png.relative_to(FIG_DIR)), SEC_HERIT,
                     "Matched-anchor S-LDSC permutation null with the real value marked.",
                     f"uv run {H}/{PERM}", VARIANT))
    for script, stems in P_OUTPUTS.items():
        for f in sorted((FIG_DIR / "permutation").iterdir() if (FIG_DIR / "permutation").is_dir()
                        else []):
            if f.is_file() and any(f.stem.startswith(s) for s in stems):
                rows.append(("Matched null", str(f.relative_to(FIG_DIR)), SEC_HERIT,
                             "Matched-anchor permutation null (Fig. 5 and supplementary).",
                             f"uv run {P}/{script}", VARIANT))
    return rows


def write_index(rows):
    def cell(s):
        return str(s).replace("|", "\\|")

    lines = [
        "# Chapter 2 figures — index",
        "",
        f"Variant: **{VARIANT}** (anchors CD4 coloc PP.H4 > 0.8; perturbation edges "
        "padj < 0.05, |z| ≥ 2, |logFC| ≥ 0.5, strong perturbations only).",
        "Rebuild everything from an empty folder with `uv run 00_build_all.py` in "
        "`chapter2/vis/`; each command below writes to `build/` and the build "
        "script copies it here under the figure number.",
        "",
        "| Figure | File | Chapter section | Message | Command (from chapter2/vis) | Variant |",
        "|---|---|---|---|---|---|",
    ]
    lines += ["| " + " | ".join(cell(c) for c in (num, f"`{f}`", sec, msg, f"`{cmd}`", var))
              + " |" for num, f, sec, msg, cmd, var in rows]
    lines += [
        "",
        "## Reserved",
        "",
        "| Figure | File | Chapter section | Message | Command (from chapter2/vis) | Variant |",
        "|---|---|---|---|---|---|",
        f"| Fig. 4C | `main/{H}/fig_4c_…` (not yet made) | {SEC_HERIT} | Enrichment with an "
        "all-measured-genes annotation added to the model (S-LDSC not yet run). | — | "
        f"{VARIANT} |",
    ]
    (FIG_DIR / "INDEX.md").write_text("\n".join(lines) + "\n")


def text_block(fig, item, y0=0.27):
    fig.text(0.06, y0, item["num"], fontsize=14, fontweight="bold", va="top")
    fig.text(0.06, y0 - 0.03, textwrap.fill(item["message"], 95), fontsize=10.5,
             va="top", style="italic")
    n_msg = len(textwrap.wrap(item["message"], 95))
    files = ", ".join(item["files"])
    fig.text(0.06, y0 - 0.035 - 0.022 * n_msg,
             textwrap.fill(item["caption"], 115) + f"\n\nFile: {files}\nSection: "
             f"{item['section']}\nCommand: {item['cmd']}",
             fontsize=8.5, va="top", color="#2b2b2b")


def figure_page(pdf, item):
    src = next(iter(item["files"].values()))
    img = plt.imread(STAGE / Path(src).with_suffix(".png"))
    fig = plt.figure(figsize=(8.27, 11.69))
    ax = fig.add_axes([0.05, 0.30, 0.90, 0.66])
    ax.imshow(img)
    ax.set_axis_off()
    text_block(fig, item)
    pdf.savefig(fig)
    plt.close(fig)


def table_page(pdf, item):
    csv = next(f for f in item["files"].values() if f.endswith(".csv"))
    df = pd.read_csv(STAGE / csv)
    df = df.apply(lambda c: c.map(lambda v: f"{v:,}" if isinstance(v, int) and v >= 10000
                                  else v))
    fig = plt.figure(figsize=(8.27, 11.69))
    ax = fig.add_axes([0.03, 0.32, 0.94, 0.64])
    ax.set_axis_off()
    cols = [textwrap.fill(c, 14) for c in df.columns]
    tab = ax.table(cellText=df.values, colLabels=cols, loc="upper center",
                   cellLoc="left", colLoc="left")
    tab.auto_set_font_size(False)
    tab.set_fontsize(7)
    tab.auto_set_column_width(range(len(cols)))
    tab.scale(1, 1.5)
    for (r, _), c in tab.get_celld().items():
        c.set_edgecolor("#d6d5d0")
        if r == 0:
            c.set_height(c.get_height() * 2.2)
            c.set_text_props(fontweight="bold")
    text_block(fig, item)
    pdf.savefig(fig)
    plt.close(fig)


def file_list_pages(pdf, rows):
    per_page = 48
    for start in range(0, len(rows), per_page):
        fig = plt.figure(figsize=(8.27, 11.69))
        fig.text(0.06, 0.96, "All files in figures/" +
                 (" (continued)" if start else ""), fontsize=13, fontweight="bold",
                 va="top")
        y = 0.92
        for num, f, *_ in rows[start:start + per_page]:
            fig.text(0.06, y, num, fontsize=8, va="top")
            fig.text(0.30, y, f, fontsize=8, va="top", family="monospace")
            y -= 0.018
        pdf.savefig(fig)
        plt.close(fig)


def write_contact_sheet(rows):
    with PdfPages(FIG_DIR / "contact_sheet.pdf") as pdf:
        for item in MANIFEST:  # main first, then supplementary (manifest order)
            (table_page if item["kind"] == "table" else figure_page)(pdf, item)
        file_list_pages(pdf, rows)


def verify(rows):
    listed = {r[1] for r in rows}
    present = {str(p.relative_to(FIG_DIR)) for p in FIG_DIR.rglob("*")
               if p.is_file() and not p.name.startswith(".")}  # skip .DS_Store etc.
    missing, extra = listed - present, present - listed
    if missing or extra:
        raise SystemExit(f"Index mismatch. Missing: {sorted(missing)}; "
                         f"not indexed: {sorted(extra)}")
    print(f"OK: {len(present)} files in {FIG_DIR}, all indexed and on the contact sheet.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--schematic", default=None,
                    help="Raster image for Fig. 2.1A (default: 07_fig_2_1.py's, "
                         "assets/ch2_schematic.png; 'none' for the placeholder).")
    args = ap.parse_args()
    schematic = (args.schematic if args.schematic in (None, "none")
                 else str(Path(args.schematic).resolve()))

    for d in (FIG_DIR, STAGE):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)

    run_scripts(schematic)
    copy_outputs()
    rows = index_rows()
    write_index(rows)
    write_contact_sheet(rows)
    verify(rows)


if __name__ == "__main__":
    main()
