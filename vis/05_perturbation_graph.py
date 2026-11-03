#!/usr/bin/env python
"""Figure 2.2 - the perturbation-response graph and expanded sets.

Message: a few knockdowns affect many genes and some genes respond to many
knockdowns, so expanded sets are large.

  A. Out-degree (DE genes per knockdown) per condition, as a complementary
     CDF P(X >= k) on log-log axes; the largest hubs are labelled.
  B. In-degree (knockdowns that change a gene), same format.
  C. Genome coverage: % of reference SNPs in the anchor / regulator /
     target annotations, from Prop._SNPs in the S-LDSC .results files;
     violins (density on log10 scale), one point per disease, median bar.
  D. Regulator and target set sizes per disease x condition (log axis), with
     the anchor (coloc) set size as a reference mark. Dashed lines in the
     regulator panel mark the regulator universe (all knockdowns with an
     effect) per condition. Empty target sets sit in a separate "0" strip
     left of a broken axis.
  E. Share of each expanded set contributed by its single largest anchor:
     max_a |set(a) & S| / |S|, where set(a) is the regulators of anchor a
     (knockdowns that change a) or the targets of anchor a (genes changed
     by knocking a down). Empty sets are skipped.

Drawn at print size (6.3 in wide, text 7-9 pt): A, B, C on the first row;
D and E on the second row, sharing the disease labels, with a wide gap
between them. One condition legend at the top serves all panels; the
legend for D's reference marks sits inside panel D.
Edges come from common.load_edges(), the exact gene-set edge definition.
Also writes the data behind each panel: degrees.csv (A, B),
genome_coverage.csv (C), set_sizes.csv (D), largest_anchor_share.csv (E) and
overview_numbers.csv (numbers quoted in the text).
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
from matplotlib.lines import Line2D

from common import (ANNOT_COLORS, ANNOT_LABEL, ANNOTS, COND_COLORS,
                    COND_LABEL, COND_SHORT, CONDS, DATA, GROUP_COLORS,
                    GROUP_LABEL, GROUPS, H4_THRESHOLD, PRINT_W, PT, axes_in,
                    color_ticklabels, disease_order, load_coloc, load_diseases,
                    load_edges, load_gene_sets, name, panel_label,
                    plain_log_axis, plain_log10_ticks, print_style, save, shade_groups,
                    size_table, universe, add_variant_args, resolve_variant)

N_HUB_LABELS = 3
SHARE_LABEL_MIN = 0.5
FIG_H = 7.1
DOT = dict(s=9, edgecolor="white", linewidth=0.5, zorder=3)
OFFS = {c: (i - 1) * 0.24 for i, c in enumerate(CONDS)}


def degrees(edges, node_col, other_col):
    return (edges.group_by("culture_condition", node_col)
            .agg(pl.col(other_col).n_unique().alias("degree"))
            .to_pandas())


def ccdf(ax, deg, name_col, xlabel):
    for cond in CONDS:
        k = np.sort(deg.loc[deg.culture_condition == cond, "degree"].values)
        p = 1 - np.arange(len(k)) / len(k)          # P(X >= k)
        ax.step(k, p, where="post", color=COND_COLORS[cond], lw=1.2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    plain_log_axis(ax, "x")
    plain_log_axis(ax, "y")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("P(degree ≥ k)")
    ax.grid(which="major")
    top = deg.sort_values("degree", ascending=False).drop_duplicates(name_col)
    n = deg.groupby("culture_condition").size()
    for rank, (_, r) in enumerate(top.head(N_HUB_LABELS).iterrows()):
        # Point on its own condition's curve; labels in a fixed stack with
        # leader lines so hubs with near-equal degree never overprint.
        d = deg.loc[deg.culture_condition == r.culture_condition, "degree"]
        p = (d >= r.degree).sum() / n[r.culture_condition]
        ax.annotate(f"{r[name_col]} ({COND_SHORT[r.culture_condition]})",
                    (r.degree, p), xytext=(0.04, 0.36 - 0.1 * rank),
                    textcoords="axes fraction", ha="left", va="center",
                    fontsize=PT["tick"], color="#52514e",
                    arrowprops=dict(arrowstyle="-", color="#a3a29d", lw=0.5))


def read_coverage(long_csv, diseases):
    """% reference SNPs per fitted set, from heritability/01_collect_results.py.

    Reading the tidy table (not the .results files) keeps the four runs where
    ldsc wrote the single surviving set as an unnamed "L2_1" category
    (osteoarthritis regulators, Alzheimer's proxy 48 h targets).
    """
    if not Path(long_csv).exists():
        raise SystemExit(f"{long_csv} not found; run heritability/01_collect_results.py first")
    df = pd.read_csv(long_csv)
    df = df[(df.status == "fitted") & df.study.isin(diseases.index) & df.annot.isin(ANNOTS)]
    return (df.assign(pct_snps=100 * df.prop_snps)
            [["study", "condition", "annot", "pct_snps"]].reset_index(drop=True))


def largest_anchor_share(edges, sets, diseases, symbols):
    """One row per non-empty (study, condition, annot) expanded set."""
    rows = []
    for cond in CONDS:
        e = edges.filter(pl.col("culture_condition") == cond)
        regs_of = dict(e.group_by("trans_ensg").agg("cis_ensg").iter_rows())
        tgts_of = dict(e.group_by("cis_ensg").agg("trans_ensg").iter_rows())
        for study in diseases.index:
            anchors = sets[(study, cond, "H4_genes")]
            for annot, per_anchor in [("regulators_H4", regs_of), ("targets_H4", tgts_of)]:
                s = sets[(study, cond, annot)]
                if not s:
                    continue
                share = {a: len(set(per_anchor.get(a, [])) & s) / len(s) for a in anchors}
                top = max(share, key=share.get)
                rows.append({"study": study, "condition": cond, "annot": annot,
                             "n_genes": len(s), "top_anchor": symbols.get(top, top),
                             "top_anchor_ensg": top, "top_share": share[top]})
    return pd.DataFrame(rows)


def disease_rows(ax, order, diseases, labels=False):
    shade_groups(ax, order, diseases, label=False)
    ax.grid(axis="x")
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.5, len(order) - 0.5)
    ax.set_yticks(range(len(order)))
    if labels:
        ax.set_yticklabels([name(diseases, s) for s in order])
        color_ticklabels(ax.get_yticklabels(), order, diseases)
    else:
        ax.tick_params(axis="y", labelleft=False)


def break_marks(ax_left, ax_right):
    """Diagonal // marks on the bottom spine between a zero strip and a log axis."""
    fig = ax_left.figure
    W, H = fig.get_size_inches()
    dx, dy = 0.03 / W, 0.05 / H        # mark half-size in figure fraction
    y0 = ax_left.get_position().y0
    for x in (ax_left.get_position().x1, ax_right.get_position().x0):
        fig.add_artist(plt.Line2D([x - dx, x + dx], [y0 - dy, y0 + dy],
                                  color="#52514e", lw=0.6))


def size_dots(ax, values, y):
    for cond in CONDS:
        ax.scatter(values[cond], y + OFFS[cond], color=COND_COLORS[cond], **DOT)


def coverage_panel(ax, cov, diseases):
    """Violins of % reference SNPs per annotation (and condition).

    Densities are estimated on log10(% SNPs) and drawn on a log10 axis, so
    the violin shapes are not distorted by the log scale. Points: one per
    disease, coloured by class; black bar: median.
    """
    rng = np.random.default_rng(0)
    annots = ["H4_genes", "regulators_H4", "targets_H4"]
    for i, annot in enumerate(annots):
        for j, cond in enumerate(CONDS):
            d = cov[(cov.annot == annot) & (cov.condition == cond)]
            if annot == "H4_genes" and j > 0:
                continue  # anchors do not depend on condition
            x0 = i + (0 if annot == "H4_genes" else (j - 1) * 0.27)
            ly = np.log10(d.pct_snps.values)
            vp = ax.violinplot(ly, positions=[x0], widths=0.24, showextrema=False)
            for body in vp["bodies"]:
                body.set_facecolor("#e4e3df")
                body.set_edgecolor("#a3a29d")
                body.set_linewidth(0.5)
                body.set_alpha(1)
            med = np.median(ly)
            ax.plot([x0 - 0.07, x0 + 0.07], [med, med], color="#0b0b0b", lw=1.0,
                    zorder=4)
            g = d.study.map(diseases["group"])
            ax.scatter(x0 + rng.uniform(-0.05, 0.05, len(d)), ly, s=6,
                       c=[GROUP_COLORS[x] for x in g], edgecolor="white",
                       linewidth=0.5, zorder=3)
    lo, hi = np.log10(cov.pct_snps.min()), np.log10(cov.pct_snps.max())
    ax.set_ylim(lo - 0.25, hi + 0.55)          # headroom for the legend
    plain_log10_ticks(ax, "y")
    ax.set_xticks(range(len(annots)))
    ax.set_xticklabels(["Anchors", "Regulators\nR  8h  48h", "Targets\nR  8h  48h"])
    ax.set_xlim(-0.45, 2.45)
    ax.set_ylabel("Reference SNPs (%)")
    ax.grid(axis="y")
    ax.set_axisbelow(True)
    short = {"IMD": "IMD", "non-IMD": "Non-IMD control"}
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=GROUP_COLORS[g],
                              ms=3.5, label=short[g]) for g in GROUPS],
              loc="upper left", handletextpad=0.2, labelspacing=0.3)


def draw_figure(order, diseases, wide, anchors, reg_universe, out_deg, in_deg,
                cov, share):
    W, H = PRINT_W, FIG_H
    fig = plt.figure(figsize=(W, H))
    n = len(order)
    y = np.arange(n)

    # Shared condition legend (A, B, D, E).
    fig.legend(handles=[Line2D([], [], color=COND_COLORS[c], marker="o", ms=3.5,
                               lw=1.2, label=COND_LABEL[c]) for c in CONDS],
               loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.0),
               columnspacing=1.5)

    # Row 1: A, B (degree CCDFs), C (genome coverage).
    row_y, row_h = 4.85, 1.65
    ax_a = axes_in(fig, 0.58, row_y, 1.45, row_h)
    ax_b = axes_in(fig, 2.62, row_y, 1.45, row_h)
    ax_c = axes_in(fig, 4.68, row_y, 1.56, row_h)
    ccdf(ax_a, out_deg, "cis_gene", "Genes changed per knockdown, k")
    ccdf(ax_b, in_deg, "trans_gene", "Knockdowns changing the gene, k")
    coverage_panel(ax_c, cov, diseases)
    top = row_y + row_h + 0.3
    panel_label(fig, 0.02, top, "A", "Out-degree")
    panel_label(fig, 2.12, top, "B", "In-degree")
    panel_label(fig, 4.2, top, "C", "Genome coverage")

    # Row 2: D (set sizes) and E (largest-anchor share) share the disease rows.
    r2_y, r2_h = 0.42, 3.0
    ax_d1 = axes_in(fig, 1.6, r2_y, 0.95, r2_h)
    ax_z = axes_in(fig, 2.71, r2_y, 0.09, r2_h, sharey=ax_d1)
    ax_d2 = axes_in(fig, 2.84, r2_y, 0.95, r2_h, sharey=ax_d1)
    ax_e1 = axes_in(fig, 4.08, r2_y, 0.58, r2_h, sharey=ax_d1)
    ax_e2 = axes_in(fig, 4.71, r2_y, 0.58, r2_h, sharey=ax_d1)

    for ax, annot in [(ax_d1, "regulators_H4"), (ax_d2, "targets_H4")]:
        disease_rows(ax, order, diseases, labels=(ax is ax_d1))
        vals = {c: np.array([wide.loc[(s, c), annot] for s in order], float) for c in CONDS}
        size_dots(ax, {c: np.where(v > 0, v, np.nan) for c, v in vals.items()}, y)
        ax.scatter([anchors[s] for s in order], y, marker="|", s=30, lw=1.2,
                   color=ANNOT_COLORS["H4_genes"], zorder=4)
        ax.set_xscale("log")
        ax.set_xlabel("Genes in set")
        ax.set_title(ANNOT_LABEL[annot])
        if annot == "targets_H4":
            zero = {c: np.where(v == 0, 0.0, np.nan) for c, v in vals.items()}
    for cond in CONDS:
        ax_d1.axvline(reg_universe[cond], color=COND_COLORS[cond], ls="--",
                      lw=0.8, zorder=2)
    ax_d1.set_xlim(1, max(reg_universe.values()) * 1.4)
    ax_d2.set_xlim(1, None)
    plain_log_axis(ax_d1, "x")
    plain_log_axis(ax_d2, "x")

    # Zero strip for empty target sets, separated by a broken axis.
    disease_rows(ax_z, order, diseases)
    size_dots(ax_z, zero, y)
    ax_z.set_xlim(-1, 1)
    ax_z.set_xticks([0])
    ax_z.set_xticklabels(["0"])
    ax_z.grid(False)
    for sp in ("left", "right"):
        ax_z.spines[sp].set_visible(False)
    ax_d2.spines["left"].set_visible(False)
    break_marks(ax_z, ax_d2)

    for ax, annot in [(ax_e1, "regulators_H4"), (ax_e2, "targets_H4")]:
        disease_rows(ax, order, diseases)
        sh = share[share.annot == annot].set_index(["study", "condition"])
        for cond in CONDS:
            pts = [(i + OFFS[cond], sh.loc[(s, cond), "top_share"])
                   for i, s in enumerate(order) if (s, cond) in sh.index]
            if pts:
                yy, xx = zip(*pts)
                ax.scatter(xx, yy, color=COND_COLORS[cond], **DOT)
        if annot == "targets_H4":
            for i, s in enumerate(order):
                d = sh[(sh.index.get_level_values(0) == s) & (sh.top_share > SHARE_LABEL_MIN)]
                if len(d):
                    ax.text(1.06, i, "/".join(dict.fromkeys(d.top_anchor)),
                            va="center", fontsize=PT["tick"], color="#52514e",
                            transform=ax.get_yaxis_transform())
        ax.set_xlim(-0.05, 1.07)
        ax.set_xticks([0, 0.5, 1])
        ax.set_xticklabels(["0", "0.5", "1"])
        ax.set_xlabel("Share")
        ax.set_title(ANNOT_LABEL[annot])

    top2 = r2_y + r2_h + 0.62
    panel_label(fig, 0.02, top2, "D", "Expanded set sizes")
    panel_label(fig, 3.95, top2, "E", "Share of set from largest anchor")
    # Legend for D's reference marks, inside panel D under its title.
    fig.legend(handles=[Line2D([], [], marker="|", ls="", color=ANNOT_COLORS["H4_genes"],
                               ms=6, mew=1.2, label="Anchor set size"),
                        Line2D([], [], ls="--", color="#6f6e69", lw=0.8,
                               label="All knockdowns with an effect")],
               loc="upper left", bbox_to_anchor=(1.57 / W, (top2 - 0.2) / H),
               ncol=2, columnspacing=1.0, handlelength=1.4, borderaxespad=0)
    return fig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--de-results", default=DATA / "DE_results_strong.csv")
    add_variant_args(ap)
    args = resolve_variant(ap.parse_args())
    out_dir = Path(args.out_dir)

    diseases = load_diseases()
    edges = load_edges(args.de_results, args.min_abs_logfc)
    out_deg = degrees(edges, "cis_gene", "trans_ensg")
    in_deg = degrees(edges, "trans_gene", "cis_ensg")
    reg_universe = {c: len(universe(edges, c)[0]) for c in CONDS}
    pd.concat([out_deg.rename(columns={"cis_gene": "gene"}).assign(direction="out"),
               in_deg.rename(columns={"trans_gene": "gene"}).assign(direction="in")]
              )[["direction", "culture_condition", "gene", "degree"]].to_csv(
        out_dir / "degrees.csv", index=False)

    sets = load_gene_sets(args.gene_sets_dir)
    sizes = size_table(sets, diseases)
    sizes.sort_values(["study", "condition", "annot"]).to_csv(
        out_dir / "set_sizes.csv", index=False)
    wide = sizes.pivot_table(index=["study", "condition"], columns="annot",
                             values="n_genes", fill_value=0)
    anchors = wide.xs("Rest", level="condition")["H4_genes"]
    order = disease_order(diseases, anchors)[::-1]
    cov = read_coverage(args.out_dir / "heritability" / "heritability_long.csv", diseases)
    cov.to_csv(out_dir / "genome_coverage.csv", index=False)

    # Anchor symbols from coloc_dat (anchors need not appear in the edge
    # list); gn2 always holds the symbol, Gene_name is sometimes the ENSG.
    symbols = dict(load_coloc().filter(pl.col("H4") > H4_THRESHOLD)
                   .select("Gene_id", "gn2").unique("Gene_id").iter_rows())
    share = largest_anchor_share(edges, sets, diseases, symbols)
    share.to_csv(out_dir / "largest_anchor_share.csv", index=False)

    with print_style():
        fig = draw_figure(order, diseases, wide, anchors, reg_universe,
                          out_deg, in_deg, cov, share)
    save(fig, "perturbation_graph_and_sets", out_dir, exact=True)

    # Numbers quoted in the text.
    rows = []
    for cond in CONDS:
        regs = wide.xs(cond, level="condition")["regulators_H4"]
        pct = 100 * regs / reg_universe[cond]
        od = out_deg.loc[out_deg.culture_condition == cond, "degree"]
        idg = in_deg.loc[in_deg.culture_condition == cond, "degree"]
        c = cov[cov.condition == cond].groupby("annot").pct_snps.median()
        rows.append({
            "condition": cond,
            "regulator_universe": reg_universe[cond],
            "regulator_set_pct_universe_median": pct.median(),
            "regulator_set_pct_universe_max": pct.max(),
            "out_degree_median": od.median(),
            "out_degree_p99": od.quantile(0.99),
            "out_degree_max": od.max(),
            "in_degree_median": idg.median(),
            "in_degree_p99": idg.quantile(0.99),
            "in_degree_max": idg.max(),
            "pct_snps_median_anchors": c.get("H4_genes"),
            "pct_snps_median_regulators": c.get("regulators_H4"),
            "pct_snps_median_targets": c.get("targets_H4"),
        })
    nums = pd.DataFrame(rows)
    nums.to_csv(out_dir / "overview_numbers.csv", index=False)
    print(nums.round(2).T.to_string(header=False))

    empty = sizes[(sizes.n_genes == 0) & (sizes.annot != "H4_genes")]
    print("Empty expanded sets (zero strip in D, skipped in E):")
    print(empty.groupby(["annot", "disease"]).condition.apply(", ".join).to_string())


if __name__ == "__main__":
    main()
