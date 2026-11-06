#!/usr/bin/env python
"""Figure and printed summary for overlap_nodeg.py.

Reads data/permutations/overlap_nodeg/overlap_summary.csv and writes
overlap_deg75_vs_nodeg.png there: one panel per set (regulators, targets),
recovery of the real set by matched pseudo-anchors with the >= 75-DEG filter
(x) against without it (y), one point per trait x condition (median, 95%
range), AITD and ASTCO at Rest highlighted.

  uv run python scripts/overlap_nodeg_report.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import polars as pl  # noqa: E402

# fig_common.py moved to vis/permutation/pcommon.py; it imports vis/common.py,
# which sets the chapter rcParams, so reset to matplotlib's defaults as before.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "vis" / "permutation"))
from pcommon import SET_COLOURS  # noqa: E402
matplotlib.rcdefaults()

OUT = Path(__file__).resolve().parents[2] / "data/permutations/overlap_nodeg"
HIGHLIGHT = {"AITD", "ASTCO"}
# Label offsets in points, per (trait, set), so labels clear the other points.
LABEL_OFFSET = {("AITD", "regulators"): (9, -3), ("ASTCO", "regulators"): (9, -10),
                ("AITD", "targets"): (9, -3), ("ASTCO", "targets"): (2, 10)}


def wide(summary: pl.DataFrame) -> pl.DataFrame:
    """One row per trait x condition x set, deg75 and nodeg columns side by side."""
    keys = ["study", "trait", "group", "condition", "set"]
    cols = ["real_size", "recovery_median", "recovery_q025", "recovery_q975", "pseudo_size_median"]
    d, n = (summary.filter(pl.col("variant") == v).select(keys + cols) for v in ("deg75", "nodeg"))
    return d.join(n, on=keys, suffix="_nodeg")


def plot(w: pl.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(6.3, 3.3), constrained_layout=True)
    for ax, s in zip(axes, ["regulators", "targets"]):
        sub = w.filter(pl.col("set") == s).drop_nulls("recovery_median")
        hl = pl.col("trait").is_in(list(HIGHLIGHT)) & (pl.col("condition") == "Rest")
        for part, alpha, size, edge in [(sub.filter(~hl), 0.35, 16, "none"),
                                        (sub.filter(hl), 1.0, 40, "black")]:
            x, y = part["recovery_median"], part["recovery_median_nodeg"]
            ax.errorbar(x, y, xerr=[x - part["recovery_q025"], part["recovery_q975"] - x],
                        yerr=[y - part["recovery_q025_nodeg"], part["recovery_q975_nodeg"] - y],
                        fmt="none", ecolor=SET_COLOURS[s], alpha=alpha * 0.6, lw=0.8)
            ax.scatter(x, y, s=size, color=SET_COLOURS[s], alpha=alpha, edgecolor=edge,
                       linewidth=0.8, zorder=3)
        for r in sub.filter(hl).iter_rows(named=True):
            ax.annotate(f"{r['trait']} rest", (r["recovery_median"], r["recovery_median_nodeg"]),
                        xytext=LABEL_OFFSET[r["trait"], s], textcoords="offset points", fontsize=8)
        ax.plot([0, 1], [0, 1], color="#9a9a9a", lw=0.8, ls="--", zorder=0)
        ax.set(xlim=(-0.02, 1.02), ylim=(-0.02, 1.02), aspect="equal",
               title=s.capitalize(), xlabel="Recovery, ≥ 75 DEGs")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Recovery, no DEG filter")
    fig.savefig(path, dpi=300)


def main() -> None:
    summary = pl.read_csv(OUT / "overlap_summary.csv")
    w = wide(summary)
    plot(w, OUT / "overlap_deg75_vs_nodeg.png")
    pct = lambda c: (pl.col(c) * 100).round(1)  # noqa: E731

    with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=200):
        print("deg75 / nodeg recovery at Rest, AITD and ASTCO (Fig. 5C: regulators 66% and 49%, targets ~8%)")
        print(w.filter(pl.col("trait").is_in(list(HIGHLIGHT)), pl.col("condition") == "Rest")
               .select("trait", "set", "real_size", pct("recovery_median").alias("deg75_%"),
                       pct("recovery_q025").alias("lo"), pct("recovery_q975").alias("hi"),
                       "real_size_nodeg", pct("recovery_median_nodeg").alias("nodeg_%"),
                       pct("recovery_q025_nodeg").alias("lo_n"), pct("recovery_q975_nodeg").alias("hi_n"))
               .sort("set", "trait"))
        print("\nMedian across IMDs (non-empty real sets) of per-trait median recovery, and of real set size")
        print(w.filter(pl.col("group") == "IMD").group_by("set", "condition").agg(
                  # Empty real sets give NaN recovery, which polars sorts above every
                  # number, so drop NaNs before the median.
                  pl.col("recovery_median").is_not_nan().sum().alias("n_traits_deg75"),
                  pl.col("recovery_median_nodeg").is_not_nan().sum().alias("n_traits_nodeg"),
                  pct("recovery_median").drop_nans().median().alias("deg75_%"),
                  pct("recovery_median_nodeg").drop_nans().median().alias("nodeg_%"),
                  pl.col("real_size").median().alias("real_size_deg75"),
                  pl.col("real_size_nodeg").median().alias("real_size_nodeg"),
                  pl.col("pseudo_size_median").median().alias("pseudo_deg75"),
                  pl.col("pseudo_size_median_nodeg").median().alias("pseudo_nodeg"))
              .sort("set", pl.col("condition").replace_strict({"Rest": 0, "Stim8hr": 1, "Stim48hr": 2})))


if __name__ == "__main__":
    main()
