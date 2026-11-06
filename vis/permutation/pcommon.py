"""Shared conventions and data loading for the matched-null figures (Chapter 2).

Paths and the print width come from vis/common.py (design setting CH2_VARIANT):
per-perm results are read from PERM_DIR and figures written to
FIG_DIR/permutation/. The rest is specific to these figures: fonts 7-9 pt,
disease display names, plain-number log axes, one colour per gene set
(anchors / regulators / targets; Okabe-Ito).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import polars as pl  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator  # noqa: E402

from common import DATA, FIG_DIR, PERM_DIR, PRINT_W  # noqa: E402

HERIT = PERM_DIR
OUT = FIG_DIR / "permutation"
# Inputs outside results/ (left in place; fig_s_matching_check.py only).
MATCH = DATA / "permutations/matching_check"
SOURCE = DATA / "permutations/figures/source_data"
MATCHING_INPUTS = [SOURCE / "matching_check_knn_anchor_values.csv",
                   SOURCE / "matching_check_knn_draws.csv",
                   SOURCE / "knn_neighbours.csv",
                   MATCH / "matching_check_draws.csv"]

FIG_WIDTH = PRINT_W  # inches, final print size

SETS = ["anchors", "regulators", "targets"]
ANNOT_OF = {"anchors": "H4_genes", "regulators": "regulators_H4",
            "targets": "targets_H4"}
SET_COLOURS = {"anchors": "#D55E00", "regulators": "#0072B2",
               "targets": "#009E73"}
SET_LABELS = {"anchors": "Anchors (coloc genes)", "regulators": "Regulators",
              "targets": "Targets"}
NULL_GREY = "#9a9a9a"

DISEASE_NAMES = {"AITD": "Autoimmune thyroid disease",
                 "ASTCO": "Asthma (childhood-onset)"}
STUDY_SHORT = {"AITD_Saevarsdottir_32581359_1": "AITD",
               "ASTCO_Ferreira_30929738_1": "ASTCO"}
COND_NAMES = {"Rest": "rest", "Stim8hr": "8 h", "Stim48hr": "48 h"}
# Run order in every figure / table.
RUNS = [("AITD_Saevarsdottir_32581359_1", "Rest"),
        ("AITD_Saevarsdottir_32581359_1", "Stim8hr"),
        ("AITD_Saevarsdottir_32581359_1", "Stim48hr"),
        ("ASTCO_Ferreira_30929738_1", "Rest"),
        ("ASTCO_Ferreira_30929738_1", "Stim48hr")]
STATS = {"enrichment": "Enrichment", "coefficient_z": "Coefficient z-score"}

# Enrichment is a ratio and can be <= 0 in the null; on log axes values below
# ENR_FLOOR are drawn at the floor and labelled "≤0.1".
ENR_FLOOR = 0.1
ENR_TICKS = [0.1, 0.3, 1, 3, 10, 30, 100]


def run_label(study: str, cond: str) -> str:
    return f"{DISEASE_NAMES[STUDY_SHORT[study]]}, {COND_NAMES[cond]}"


def setup_style():
    # Importing common applies the chapter rcParams; these figures start from
    # matplotlib's defaults, as they did before the move.
    matplotlib.rcdefaults()
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300,
    })


def _plain(x, _pos=None):
    if x >= 1:
        return f"{x:,.0f}"
    return f"{x:g}"


def plain_log_axis(ax, axis="x", ticks=None, floor_label=None):
    """Log scale with plain-number tick labels; optional '≤floor' first label."""
    a = ax.xaxis if axis == "x" else ax.yaxis
    (ax.set_xscale if axis == "x" else ax.set_yscale)("log")
    if ticks is not None:
        a.set_major_locator(FixedLocator(ticks))
    fmt = (lambda x, p=None: f"≤{x:g}" if floor_label is not None
           and np.isclose(x, floor_label) else _plain(x))
    a.set_major_formatter(FuncFormatter(fmt))
    a.set_minor_locator(NullLocator())


def panel_letter(ax, letter, x=-0.02, y=1.02, fig=None):
    (fig or ax.figure).text(x, y, letter, transform=ax.transAxes, fontsize=9,
                            fontweight="bold", ha="right", va="bottom")


def fmt_p(p: float) -> str:
    return f"p = {p:.3f}" if p >= 0.001 else f"p = {p:.1e}"


def emp_p(real: float, null: np.ndarray) -> float:
    null = null[~np.isnan(null)]
    return float((1 + np.sum(null >= real)) / (null.size + 1))


def load_perperm(study: str, cond: str) -> pl.DataFrame:
    return pl.read_csv(HERIT / f"{study}_{cond}" / "permutation_ldsc_perperm.csv",
                       infer_schema_length=None)


def null_and_real(df: pl.DataFrame, set_name: str, stat: str):
    sub = df.filter(pl.col("annot") == ANNOT_OF[set_name])
    real = sub.filter(pl.col("is_real"))[stat].item(0)
    null = sub.filter(~pl.col("is_real"))[stat].drop_nulls().to_numpy().astype(float)
    return real, null


def read_set(path: Path) -> set:
    return {l.strip() for l in path.read_text().splitlines() if l.strip()} \
        if path.exists() else set()


def recovery(cond: str = "Rest", perm_root: Path | None = None) -> pl.DataFrame:
    """Fraction of the real regulator / target set recovered by each matched
    permutation: |perm set & real set| / |real set|, from the gene sets in
    PERM_DIR/<study>/perm*/gene_sets (real = permreal). Same columns as the
    old source_data/matching_check_knn_draws.csv (method "knn"). perm_root
    overrides PERM_DIR (fig_s_deg_filter reads both designs)."""
    rows = []
    for study in dict.fromkeys(s for s, c in RUNS if c == cond):
        root = Path(perm_root or HERIT) / study
        def sets(perm):
            g = root / perm / "gene_sets"
            return (read_set(g / f"{study}_{cond}_regulators_H4.txt"),
                    read_set(g / f"{study}_{cond}_targets_H4.txt"))
        real_reg, real_tgt = sets("permreal")
        perms = sorted((p.name for p in root.glob("perm[0-9]*")), key=lambda n: int(n[4:]))
        for perm in perms:
            reg, tgt = sets(perm)
            rows.append({"study": STUDY_SHORT[study], "method": "knn", "draw": int(perm[4:]),
                         "reg_overlap_frac_of_real": len(reg & real_reg) / len(real_reg),
                         "tgt_overlap_frac_of_real": len(tgt & real_tgt) / len(real_tgt)})
    return pl.DataFrame(rows)


def run_stats() -> pl.DataFrame:
    """Per run x set x statistic: real, null median / 2.5% / 97.5%, empirical p.
    Empirical p and real values are checked against summary_all.csv."""
    summ = pl.read_csv(HERIT / "summary_all.csv")
    rows = []
    for study, cond in RUNS:
        df = load_perperm(study, cond)
        for s in SETS:
            for stat in STATS:
                real, null = null_and_real(df, s, stat)
                p = emp_p(real, null)
                ref = summ.filter((pl.col("study") == study) & (pl.col("condition") == cond)
                                  & (pl.col("set") == s)).row(0, named=True)
                assert np.isclose(ref[f"empirical_p_{stat}"], p)
                assert np.isclose(ref[f"real_{stat}"], real)
                rows.append({"study": study, "trait": STUDY_SHORT[study],
                             "condition": cond, "set": s, "statistic": stat,
                             "real": real, "null_median": float(np.median(null)),
                             "null_q025": float(np.quantile(null, 0.025)),
                             "null_q975": float(np.quantile(null, 0.975)),
                             "empirical_p": p, "n_null": int(null.size)})
    return pl.DataFrame(rows)


def save(fig, stem: str):
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{stem}.{ext}", bbox_inches="tight")
    print(f"Wrote {OUT / stem}.{{pdf,png}}")
