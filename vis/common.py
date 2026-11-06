"""Shared paths, disease metadata, palette, axis formatting and loaders for
the chapter 2 figures.

The analysis design is one setting, VARIANT: the folder name under
results/ldsc/ (env var CH2_VARIANT; default the deg75 design). Gene sets are
read from ``LDSC_DIR/gene_sets`` with the pipeline's
``{study}_{condition}_{annot}.txt`` naming and S-LDSC output from
``LDSC_DIR/ldsc_results``; disease display names and the IMD / non-IMD split
come from ``diseases.csv`` next to this file.

    CH2_VARIANT=deg0_cells50_... uv run 00_build_all.py

Scripts write to ``build/<VARIANT>/`` (staging); ``00_build_all.py`` copies the
chapter outputs from there into ``FIG_DIR`` under their figure numbers.
"""

import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
import polars as pl

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
# Design folder names under results/ldsc/. VARIANT is the one in use; the
# design-comparison scripts (vis/permutation/*deg_filter*) read both.
DESIGNS = {"deg0": "deg0_cells50_ot1_p0.05_z2.0_lfc0.5_h40.8_k1-1_w10000_snp1000",
           "deg75": "deg75_cells50_ot1_p0.05_z2.0_lfc0.5_h40.8_k1-1_w10000_snp1000"}
VARIANT = os.environ.get("CH2_VARIANT", DESIGNS["deg75"])


def ldsc_dir(variant):
    return RESULTS / "ldsc" / variant


def perm_dir(variant):
    return RESULTS / "ldsc_permutation" / (variant + "_bins4-5_knn20_seed42")


LDSC_DIR = ldsc_dir(VARIANT)
PERM_DIR = perm_dir(VARIANT)
FIGURES2 = ROOT / "figures2"
FIG_DIR = FIGURES2 / VARIANT.split("_")[0]   # temporary, e.g. figures2/deg75/
BUILD_DIR = ROOT / "build" / VARIANT                  # staging, copied into FIG_DIR

CONDS = ["Rest", "Stim8hr", "Stim48hr"]
COND_LABEL = {"Rest": "Rest", "Stim8hr": "Stim 8 h", "Stim48hr": "Stim 48 h"}
COND_SHORT = {"Rest": "Rest", "Stim8hr": "8 h", "Stim48hr": "48 h"}
# Ordinal ink ramp (light -> dark) for the stimulation time course.
COND_COLORS = {"Rest": "#b5b4ae", "Stim8hr": "#6f6e69", "Stim48hr": "#1f1f1d"}

ANNOTS = ["H4_genes", "targets_H4", "regulators_H4"]
ANNOT_LABEL = {"H4_genes": "Coloc genes (anchors)",
               "targets_H4": "Targets",
               "regulators_H4": "Regulators"}
ANNOT_COLORS = {"H4_genes": "#2a78d6",
                "targets_H4": "#eb6834",
                "regulators_H4": "#1baf7a"}

GROUPS = ["IMD", "non-IMD"]
GROUP_LABEL = {"IMD": "Immune-mediated diseases",
               "non-IMD": "Non-IMD controls"}
GROUP_SHORT = {"IMD": "Immune-mediated", "non-IMD": "Non-IMD control"}
GROUP_COLORS = {"IMD": "#6250d6", "non-IMD": "#8a8985"}
GROUP_LINESTYLE = {"IMD": "-", "non-IMD": "--"}
GROUP_SHADE = {"IMD": "#efedfb", "non-IMD": "#f1f1ef"}

H4_THRESHOLD = 0.8
CELL_TYPE = "CD4"

plt.rcParams.update({
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": "#52514e",
    "axes.labelcolor": "#0b0b0b",
    "xtick.color": "#52514e",
    "ytick.color": "#52514e",
    "axes.grid": False,
    "grid.color": "#e4e3df",
    "grid.linewidth": 0.6,
    "legend.frameon": False,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
})


# Print-size figures (main chapter figures): drawn at final size, 6.3 in wide,
# text 7-9 pt, lines >= 0.5 pt. Use ``with print_style():`` around drawing and
# ``save(..., exact=True)`` so the PDF keeps the exact page width.
PRINT_W = 6.3
PT = {"tick": 7, "label": 7.5, "title": 8, "panel": 9}
PRINT_RC = {
    "font.size": PT["label"],
    "axes.titlesize": PT["title"],
    "axes.labelsize": PT["label"],
    "xtick.labelsize": PT["tick"],
    "ytick.labelsize": PT["tick"],
    "legend.fontsize": PT["tick"],
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.5,
    "ytick.minor.width": 0.5,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "xtick.minor.size": 1.5,
    "ytick.minor.size": 1.5,
    "lines.linewidth": 1.0,
    "grid.linewidth": 0.5,
    "patch.linewidth": 0.5,
    "legend.handlelength": 1.6,
    "legend.borderaxespad": 0.3,
}


def print_style():
    return plt.rc_context(PRINT_RC)


def axes_in(fig, x, y, w, h, **kw):
    """Add axes at (x, y) with size (w, h), all in inches from bottom-left."""
    W, H = fig.get_size_inches()
    return fig.add_axes([x / W, y / H, w / W, h / H], **kw)


def panel_label(fig, x, y, letter, title):
    """Bold panel letter plus title at (x, y) inches (top-left anchor)."""
    W, H = fig.get_size_inches()
    fig.text(x / W, y / H, f"{letter}  ", fontsize=PT["panel"], fontweight="bold",
             ha="left", va="top")
    fig.text((x + 0.17) / W, y / H, title, fontsize=PT["panel"], ha="left", va="top")


def load_diseases():
    """Disease metadata indexed by study id, IMDs first."""
    df = pd.read_csv(HERE / "diseases.csv")
    df["group"] = pd.Categorical(df["group"], GROUPS, ordered=True)
    return df.set_index("study")


def name(diseases, study):
    """Display name used on every figure and table."""
    return diseases.loc[study, "display_name"]


def _plain_number(v, _pos=None):
    if v <= 0:
        return ""
    if v >= 1:
        return f"{v:,.0f}"
    return f"{v:.10f}".rstrip("0").rstrip(".")


def plain_log_axis(ax, axis="x"):
    """Log-axis ticks as plain numbers: 0.001, 0.1, 1, 10, 1,000, 10,000.

    Call after the axis limits are final. Every decade gets a labelled tick;
    2 and 5 are added only when there are >= 0.9 inches per decade, so ticks
    never crowd at print size.
    """
    a = ax.xaxis if axis == "x" else ax.yaxis
    lo, hi = ax.get_xlim() if axis == "x" else ax.get_ylim()
    lo, hi = min(lo, hi), max(lo, hi)
    bbox = ax.get_window_extent().transformed(ax.figure.dpi_scale_trans.inverted())
    length = bbox.width if axis == "x" else bbox.height
    decades = np.log10(hi / lo) if lo > 0 else np.inf
    subs = (1.0, 2.0, 5.0) if length / decades >= 0.9 else (1.0,)
    a.set_major_locator(LogLocator(base=10, subs=subs, numticks=50))
    a.set_major_formatter(FuncFormatter(_plain_number))
    a.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10), numticks=50))
    a.set_minor_formatter(NullFormatter())


def plain_log10_ticks(ax, axis="y"):
    """For data already plotted as log10 values on a linear axis: label each
    decade with its plain number (0.01, 0.1, 1, 10) and add 2-9 minor ticks.
    Call after the axis limits are final."""
    a = ax.yaxis if axis == "y" else ax.xaxis
    lo, hi = ax.get_ylim() if axis == "y" else ax.get_xlim()
    lo, hi = min(lo, hi), max(lo, hi)
    decades = np.arange(np.ceil(lo), np.floor(hi) + 1)
    a.set_ticks(decades)
    a.set_ticklabels([_plain_number(10 ** d) for d in decades])
    minor = np.log10(np.outer(10.0 ** np.arange(np.floor(lo), np.ceil(hi)),
                              np.arange(2, 10)).ravel())
    a.set_ticks(minor[(minor > lo) & (minor < hi)], minor=True)


def read_gene_set(path):
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def load_gene_sets(gene_sets_dir=None):
    """dict[(study, cond, annot)] = set of Ensembl ids."""
    gene_sets_dir = Path(gene_sets_dir or LDSC_DIR / "gene_sets")
    annot_re = "|".join(sorted(ANNOTS, key=len, reverse=True))
    pat = re.compile(rf"^(?P<study>.+)_(?P<cond>[^_]+)_(?P<annot>{annot_re})\.txt$")
    sets = {}
    for p in sorted(gene_sets_dir.glob("*.txt")):
        m = pat.match(p.name)
        if m and m.group("cond") in CONDS:
            sets[(m.group("study"), m.group("cond"), m.group("annot"))] = read_gene_set(p)
    if not sets:
        raise SystemExit(f"No gene-set files found in {gene_sets_dir}")
    return sets


def size_table(sets, diseases):
    """Long table: study, disease, group, condition, annot, n_genes."""
    rows = [
        {"study": s, "condition": c, "annot": a, "n_genes": len(g)}
        for (s, c, a), g in sets.items() if s in diseases.index
    ]
    df = pd.DataFrame(rows)
    df["disease"] = df["study"].map(diseases["display_name"])
    df["group"] = df["study"].map(diseases["group"]).astype(diseases["group"].dtype)
    return df


def disease_order(diseases, anchor_counts):
    """IMDs then controls, each block sorted by anchor count (descending)."""
    d = diseases.assign(n=pd.Series(anchor_counts)).fillna({"n": 0})
    return list(d.sort_values(["group", "n"], ascending=[True, False]).index)


def load_coloc(coloc_path=None):
    """CD4 coloc rows with versionless Gene_id, as in make_gene_sets_coloc.py."""
    return (
        pl.scan_csv(coloc_path or DATA / "coloc_dat.tsv", separator="\t",
                    infer_schema_length=100000)
        .filter(pl.col("cell_type") == CELL_TYPE)
        .with_columns(pl.col("Gene_id").str.replace(r"\..*", ""))
        .collect()
    )


# Edge definition used to build data/gene_sets/coloc (make_gene_sets_coloc.py
# run by the coloc module with default filters). DE_results_strong.csv was
# written by extract_de_results.py with: on-target KD, ontarget_significant,
# n_total_de_genes >= 75, n_cells_target >= 50, padj < 0.05, |zscore| >= 2,
# cis != trans. The gene-set step adds no |logFC| filter (min_abs_logfc = 0)
# and no convergence filter (min_reg_k = min_tgt_k = 1). The padj / zscore
# filters below are re-applied as no-ops to make the definition explicit.
EDGE_PADJ = 0.05
EDGE_ABS_Z = 2.0


def load_edges(de_path=None, min_abs_logfc=0.0):
    """Perturbation -> response edges exactly as the gene sets use them.

    ``min_abs_logfc`` mirrors make_gene_sets_coloc.py --min-abs-logfc
    (the ``lfc`` field of VARIANT).
    """
    return (
        pl.read_csv(de_path or DATA / "DE_results_strong.csv",
                    columns=["culture_condition", "cis_ensg", "cis_gene",
                             "trans_ensg", "trans_gene", "padj", "zscore", "logFC"])
        .filter(pl.col("padj") < EDGE_PADJ,
                pl.col("zscore").abs() >= EDGE_ABS_Z,
                pl.col("logFC").abs() >= min_abs_logfc,
                pl.col("cis_ensg") != pl.col("trans_ensg"))
        .drop("padj", "zscore", "logFC")
    )


# The |logFC| edge filter is the "lfc<x>" field of VARIANT (lfc0.5 matches
# the old coloc_thr0.5 design: H4 > 0.8 and |logFC| >= 0.5).
def variant_config():
    m = re.search(r"(?:^|_)lfc(?P<lfc>[\d.]+)(?:_|$)", VARIANT)
    if not m:
        raise SystemExit(f"CH2_VARIANT={VARIANT!r} has no lfc<x> field")
    return {
        "gene_sets_dir": LDSC_DIR / "gene_sets",
        "ldsc_dir": LDSC_DIR / "ldsc_results",
        "min_abs_logfc": float(m.group("lfc")),
        "out_dir": BUILD_DIR,
    }


def add_variant_args(ap):
    ap.add_argument("--out-dir", default=None,
                    help=f"Override the output directory (default: {BUILD_DIR}). "
                         "The design itself is set with the CH2_VARIANT env var.")


def resolve_variant(args):
    """Fill args.gene_sets_dir / ldsc_dir / min_abs_logfc / out_dir in place."""
    cfg = variant_config()
    args.gene_sets_dir = cfg["gene_sets_dir"]
    args.ldsc_dir = cfg["ldsc_dir"]
    args.min_abs_logfc = cfg["min_abs_logfc"]
    args.out_dir = Path(args.out_dir or cfg["out_dir"])
    args.out_dir.mkdir(parents=True, exist_ok=True)
    return args


def universe(edges, cond):
    """(perturbations with >= 1 edge, genes with >= 1 incoming edge) in cond.

    These are the largest possible regulator and target sets in that
    condition; no list of all measured / perturbed genes is available
    locally (the DE h5ad lives on the cluster).
    """
    e = edges.filter(pl.col("culture_condition") == cond)
    return set(e["cis_ensg"].unique()), set(e["trans_ensg"].unique())


def gene_symbols(edges):
    """Ensembl id -> symbol from both ends of the edge list."""
    a = edges.select(pl.col("cis_ensg").alias("id"), pl.col("cis_gene").alias("sym"))
    b = edges.select(pl.col("trans_ensg").alias("id"), pl.col("trans_gene").alias("sym"))
    return dict(pl.concat([a, b]).unique("id").iter_rows())


def shade_groups(ax, order, diseases, axis="y", label=True):
    """Shade the IMD / non-IMD blocks of a categorical axis and label them.

    ``order`` is the plotted order of study ids at positions 0..n-1.
    """
    groups = [diseases.loc[s, "group"] for s in order]
    for g in GROUPS:
        idx = [i for i, x in enumerate(groups) if x == g]
        if not idx:
            continue
        lo, hi = min(idx) - 0.5, max(idx) + 0.5
        span = ax.axhspan if axis == "y" else ax.axvspan
        span(lo, hi, color=GROUP_SHADE[g], zorder=0, lw=0)
        if label:
            mid = (lo + hi) / 2
            if axis == "y":
                ax.text(1.0, mid, GROUP_LABEL[g], transform=ax.get_yaxis_transform(),
                        rotation=270, ha="left", va="center", fontsize=9,
                        color=GROUP_COLORS[g], fontweight="bold")
            else:
                ax.text(mid, 1.0, GROUP_LABEL[g], transform=ax.get_xaxis_transform(),
                        ha="center", va="bottom", fontsize=9,
                        color=GROUP_COLORS[g], fontweight="bold")


def color_ticklabels(labels, order, diseases):
    for t, s in zip(labels, order):
        t.set_color(GROUP_COLORS[diseases.loc[s, "group"]])


def save(fig, stem, out_dir=None, exact=False):
    """Write PDF + PNG. ``exact`` keeps the figure size (no tight cropping)."""
    out_dir = Path(out_dir or BUILD_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    with plt.rc_context({"savefig.bbox": "standard" if exact else "tight"}):
        for ext in ("pdf", "png"):
            fig.savefig(out_dir / f"{stem}.{ext}")
    plt.close(fig)
    print(f"Wrote {out_dir / stem}.{{pdf,png}}")
