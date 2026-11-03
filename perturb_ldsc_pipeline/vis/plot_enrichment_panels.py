#!/usr/bin/env python
"""Four-panel S-LDSC heatmap figure (supervisor-requested layout).

  Panel 1 (single row):       anchor annotation (e.g. "coloc genes"),
                              condition-collapsed by mean. One cell per
                              disease.
  Panels 2-4 (2 rows each):   one per condition (Rest, Stim8hr, Stim48hr).
                              Rows = (Regulators, Targets), columns = diseases.

The plotted metric is selectable with --metric:

  enrichment     (default) colour = log2(Enrichment); stars = enrichment_p FDR.
                 This is the MARGINAL metric: a star means enrichment != 1, NOT
                 that the annotation adds heritability beyond baseline.
  coefficient_z  colour = Coefficient z-score (raw, diverging scale centred 0);
                 stars = one-sided coefficient-z FDR (p = P(Z >= z)). This is
                 the CONDITIONAL metric: a star means the annotation adds
                 heritability beyond the baseline model.

Generating both and showing them side by side makes the point that most
enrichment stars do not survive the conditional test.

All cells in a panel set share one colour scale. BH-FDR is computed across
every visible cell; significant cells get a star.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe
import numpy as np
import polars as pl
from scipy.stats import norm

COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]

# Immune-mediated disease short codes (same set as collect_sensitivity_ldsc.py /
# the drug-target heatmaps) — used to group IMDs left of non-IMDs.
IMD_DISEASES = {
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
}

# Readable axis labels keyed by short code (same as the drug-target heatmaps).
SHORT_NAMES = {
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


def pretty_label(label):
    """Readable disease label: map a study ID / short code through SHORT_NAMES
    (falling back to the short-code prefix, then the original label)."""
    if label in SHORT_NAMES:
        return SHORT_NAMES[label]
    return SHORT_NAMES.get(str(label).split("_")[0], label)


def order_diseases_imd(df, label_col, diseases):
    """Return (ordered_diseases, n_imd) with IMDs (sorted) left of non-IMDs
    (sorted), matching the drug-target heatmap layout. Uses an `is_imd` column
    if present, else classifies via trait_short / study-prefix short code."""
    if "is_imd" in df.columns:
        truthy = (
            pl.col("is_imd").cast(pl.Utf8).str.to_lowercase().is_in(["true", "1", "t"])
        )
        imd_labels = set(df.filter(truthy)[label_col].unique().to_list())
        member = {d: d in imd_labels for d in diseases}
    else:
        code_col = "trait_short" if "trait_short" in df.columns else label_col
        codes = {
            r[label_col]: r["code"]
            for r in df.group_by(label_col)
            .agg(pl.col(code_col).first().alias("code"))
            .iter_rows(named=True)
        }

        def is_imd(d):
            c = codes.get(d, d)
            return c in IMD_DISEASES or str(c).split("_")[0] in IMD_DISEASES

        member = {d: is_imd(d) for d in diseases}
    imd = sorted(d for d in diseases if member[d])
    non = sorted(d for d in diseases if not member[d])
    return imd + non, len(imd)


# Per-metric configuration. `p_col=None` means derive a one-sided p from the
# value itself (used for coefficient_z: p = P(Z >= z)).
METRICS = {
    "enrichment": {
        "value_col": "enrichment",
        "p_col": "enrichment_p",
        "log2": True,
        "cbar": "log2(Enrichment)",
        "vmin": 0.0,
        "vmax": 6.0,
        "cmap": "Reds",
    },
    "coefficient_z": {
        "value_col": "coefficient_z",
        "p_col": None,  # derive one-sided p from the z-score
        "log2": False,
        "cbar": "Coefficient z-score",
        "vmin": -4.0,
        "vmax": 4.0,
        "cmap": "RdBu_r",
    },
}


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan)
    mask = ~np.isnan(p)
    n = mask.sum()
    if n == 0:
        return out
    idx = np.where(mask)[0]
    order = np.argsort(p[idx])
    ranked = p[idx][order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out[idx[order]] = adj
    return out


def pick_label_col(df, requested):
    if requested:
        return requested
    if (
        "traitFromSource" in df.columns
        and df["traitFromSource"].null_count() < df.height
    ):
        return "traitFromSource"
    if "trait_short" in df.columns:
        return "trait_short"
    return "study"


def aggregate(df, label_col, value_col, p_col):
    """Mean value (+ mean p if available), per (label, condition, annot)."""
    exprs = [pl.col(value_col).mean().alias("value")]
    if p_col is not None:
        exprs.append(pl.col(p_col).mean().alias("pval"))
    return df.group_by([label_col, "condition", "annot"]).agg(exprs)


def build_panel_matrices(agg, label_col, diseases, anchor, regulators, targets, has_p):
    """Returns value matrices and p matrices:
    anchor_v (1, N), anchor_p (1, N), cond_v {cond -> (2, N)}, cond_p {cond -> (2, N)}
    Rows of the condition matrices are [regulators, targets].
    If has_p is False the p matrices are derived later from the values.
    """
    lookup = {
        (r[label_col], r["condition"], r["annot"]): (
            r["value"],
            r["pval"] if has_p else None,
        )
        for r in agg.iter_rows(named=True)
    }

    N = len(diseases)
    anchor_v = np.full(N, np.nan)
    anchor_p = np.full(N, np.nan)
    for j, d in enumerate(diseases):
        hits = [lookup.get((d, c, anchor)) for c in COND_ORDER]
        hits = [h for h in hits if h is not None]
        if not hits:
            continue
        anchor_v[j] = np.mean([h[0] for h in hits])
        if has_p:
            anchor_p[j] = np.mean([h[1] for h in hits])

    cond_v = {c: np.full((2, N), np.nan) for c in COND_ORDER}
    cond_p = {c: np.full((2, N), np.nan) for c in COND_ORDER}
    for c in COND_ORDER:
        for j, d in enumerate(diseases):
            for i, annot in enumerate([regulators, targets]):
                hit = lookup.get((d, c, annot))
                if hit is None:
                    continue
                cond_v[c][i, j] = hit[0]
                if has_p:
                    cond_p[c][i, j] = hit[1]

    if not has_p:
        # One-sided significance from the z-score: p = P(Z >= z).
        anchor_p = norm.sf(anchor_v)
        cond_p = {c: norm.sf(cond_v[c]) for c in COND_ORDER}

    return anchor_v[None, :], anchor_p[None, :], cond_v, cond_p


def joint_fdr(anchor_p, cond_p):
    """BH-FDR across all visible cells (anchor row + 3×2×N condition cells)."""
    flat = [anchor_p.ravel()] + [m.ravel() for m in cond_p.values()]
    pooled = np.concatenate(flat)
    adj = bh_fdr(pooled)
    sizes = [m.size for m in flat]
    chunks = np.split(adj, np.cumsum(sizes)[:-1])
    anchor_fdr = chunks[0].reshape(anchor_p.shape)
    cond_fdr = {
        c: chunks[1 + i].reshape(cond_p[c].shape) for i, c in enumerate(COND_ORDER)
    }
    return anchor_fdr, cond_fdr


def _draw_panel(
    ax,
    data,
    present,
    fdr,
    fdr_thr,
    cmap,
    vmin,
    vmax,
    row_labels,
    col_labels,
    show_col_labels,
    star_color,
    star_size,
    title,
    divider_at=None,
    show_group_labels=False,
):
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad("white")
    plot = np.where(present, data, np.nan)
    im = ax.imshow(plot, cmap=cmap_obj, vmin=vmin, vmax=vmax, aspect="auto")

    n_rows, n_cols = data.shape
    for i in range(n_rows):
        for j in range(n_cols):
            if not present[i, j]:
                ax.add_patch(
                    mpatches.Rectangle(
                        (j - 0.5, i - 0.5),
                        1,
                        1,
                        facecolor="lightgrey",
                        edgecolor="white",
                        hatch="///",
                        linewidth=0,
                    )
                )
                continue
            if not np.isnan(fdr[i, j]) and fdr[i, j] < fdr_thr:
                ax.text(
                    j,
                    i,
                    "★",
                    ha="center",
                    va="center",
                    fontsize=star_size,
                    color="white",
                    path_effects=[pe.withStroke(linewidth=2, foreground="black")],
                )

    if divider_at is not None and 0 < divider_at < n_cols:
        ax.axvline(divider_at - 0.5, color="black", linewidth=2.0)
        if show_group_labels:
            ax.text(
                divider_at / 2 - 0.5,
                -0.9,
                "Immune-mediated diseases",
                ha="center",
                va="bottom",
                fontsize=12,
                fontweight="bold",
            )
            ax.text(
                (divider_at + n_cols) / 2 - 0.5,
                -0.9,
                "Other diseases",
                ha="center",
                va="bottom",
                fontsize=12,
                fontweight="bold",
            )

    ax.set_xticks(range(n_cols))
    if show_col_labels:
        ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=14)
    else:
        ax.set_xticklabels([])
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_labels, fontsize=15)
    ax.tick_params(axis="both", length=0)
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.5)

    ax.set_title(title, fontsize=14, loc="left", pad=4, fontweight="bold")
    return im


def plot(
    df,
    out_path,
    label_col,
    anchor,
    regulators,
    targets,
    anchor_label,
    fdr_thr,
    metric,
    vmin,
    vmax,
    cmap,
    star_color,
    star_size,
    group_imd,
):
    cfg = METRICS[metric]
    agg = aggregate(df, label_col, cfg["value_col"], cfg["p_col"])
    diseases = sorted(agg[label_col].unique().to_list())
    if group_imd:
        diseases, n_imd = order_diseases_imd(df, label_col, diseases)
    else:
        n_imd = 0
    divider = n_imd if 0 < n_imd < len(diseases) else None
    display_labels = [pretty_label(d) for d in diseases]

    anchor_v, anchor_p, cond_v, cond_p = build_panel_matrices(
        agg,
        label_col,
        diseases,
        anchor,
        regulators,
        targets,
        has_p=cfg["p_col"] is not None,
    )
    anchor_fdr, cond_fdr = joint_fdr(anchor_p, cond_p)

    if cfg["log2"]:
        floor = 2**vmin

        def transform(arr):
            return np.where(np.isnan(arr), np.nan, np.log2(np.clip(arr, floor, None)))

    else:

        def transform(arr):
            return arr  # z-scores are plotted on a linear, signed scale

    anchor_disp = transform(anchor_v)
    cond_disp = {c: transform(cond_v[c]) for c in COND_ORDER}
    anchor_present = ~np.isnan(anchor_v)
    cond_present = {c: ~np.isnan(cond_v[c]) for c in COND_ORDER}

    N = len(diseases)
    fig_w = max(11, 0.55 * N + 3)
    fig = plt.figure(figsize=(fig_w, 10))
    gs = fig.add_gridspec(
        4,
        2,
        width_ratios=[1, 0.025],
        height_ratios=[1.0, 2.0, 2.0, 2.0],
        hspace=0.25,
        wspace=0.02,
    )
    ax_anchor = fig.add_subplot(gs[0, 0])
    ax_conds = [fig.add_subplot(gs[i, 0]) for i in range(1, 4)]
    cax = fig.add_subplot(gs[:, 1])

    im = _draw_panel(
        ax_anchor,
        anchor_disp,
        anchor_present,
        anchor_fdr,
        fdr_thr,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        row_labels=[anchor_label],
        col_labels=display_labels,
        show_col_labels=False,
        star_color=star_color,
        star_size=star_size,
        title=f"{anchor_label}",
        divider_at=divider,
        show_group_labels=divider is not None,
    )

    for ax, cond in zip(ax_conds, COND_ORDER):
        is_last = cond == COND_ORDER[-1]
        _draw_panel(
            ax,
            cond_disp[cond],
            cond_present[cond],
            cond_fdr[cond],
            fdr_thr,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            row_labels=["Regulators", "Targets"],
            col_labels=display_labels,
            show_col_labels=is_last,
            star_color=star_color,
            star_size=star_size,
            title=cond,
            divider_at=divider,
            show_group_labels=False,
        )

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(cfg["cbar"], fontsize=15)
    cbar.ax.tick_params(labelsize=13)

    fig.suptitle(
        f"S-LDSC {cfg['cbar']} — {anchor_label} anchor (|logFC| > 0.5)   ★ = FDR < {fdr_thr}",
        fontsize=16,
        fontweight="bold",
        y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}  (metric={metric})")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--summary", required=True, help="CSV from parse_ldsc_results.py"
    )
    parser.add_argument(
        "--metric",
        choices=list(METRICS),
        default="enrichment",
        help="Which S-LDSC statistic to plot (default enrichment). "
        "coefficient_z is the conditional metric.",
    )
    parser.add_argument(
        "--anchor-annot",
        default="H4_genes",
        help="Anchor annotation name (default H4_genes).",
    )
    parser.add_argument("--regulators-annot", default="regulators_H4")
    parser.add_argument("--targets-annot", default="targets_H4")
    parser.add_argument(
        "--anchor-label",
        default="coloc genes",
        help="Display name for the anchor panel " "(default 'coloc genes').",
    )
    parser.add_argument(
        "--label-col",
        default=None,
        help="Column to use as disease label "
        "(default: traitFromSource → trait_short → study).",
    )
    parser.add_argument("--fdr", type=float, default=0.05)
    # Colour-scale controls. Defaults are metric-aware when left unset; the
    # --log2-* names are kept as backward-compatible aliases.
    parser.add_argument(
        "--vmin",
        type=float,
        default=None,
        help="Colour-scale min (metric-aware default).",
    )
    parser.add_argument(
        "--vmax",
        type=float,
        default=None,
        help="Colour-scale max (metric-aware default).",
    )
    parser.add_argument(
        "--cmap", default=None, help="Matplotlib colormap (metric-aware default)."
    )
    parser.add_argument(
        "--log2-vmin", type=float, default=None, help="Deprecated alias for --vmin."
    )
    parser.add_argument(
        "--log2-vmax", type=float, default=None, help="Deprecated alias for --vmax."
    )
    parser.add_argument(
        "--log2-cmap", default=None, help="Deprecated alias for --cmap."
    )
    parser.add_argument(
        "--star-color",
        default="#7CFC00",
        help="Colour for the FDR-significant stars " "(default light green).",
    )
    parser.add_argument("--star-size", type=float, default=12.0)
    parser.add_argument(
        "--no-group-imd",
        dest="group_imd",
        action="store_false",
        help="Disable grouping IMDs left of non-IMDs (default: group, with a "
        "divider and group labels, like the drug-target heatmaps).",
    )
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    cfg = METRICS[args.metric]
    # Resolve colour-scale: explicit --vmin > deprecated --log2-vmin > default.
    vmin = next(v for v in (args.vmin, args.log2_vmin, cfg["vmin"]) if v is not None)
    vmax = next(v for v in (args.vmax, args.log2_vmax, cfg["vmax"]) if v is not None)
    cmap = next(c for c in (args.cmap, args.log2_cmap, cfg["cmap"]) if c is not None)

    df = pl.read_csv(args.summary).filter(
        pl.col("annot").is_in(
            [args.anchor_annot, args.regulators_annot, args.targets_annot]
        )
    )
    if df.height == 0:
        raise SystemExit("No rows matched the anchor/regulators/targets annots")
    if cfg["value_col"] not in df.columns:
        raise SystemExit(
            f"--metric {args.metric} needs column '{cfg['value_col']}' "
            f"in {args.summary}; present columns: {df.columns}"
        )

    label_col = pick_label_col(df, args.label_col)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    plot(
        df,
        out_path,
        label_col=label_col,
        anchor=args.anchor_annot,
        regulators=args.regulators_annot,
        targets=args.targets_annot,
        anchor_label=args.anchor_label,
        fdr_thr=args.fdr,
        metric=args.metric,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        star_color=args.star_color,
        star_size=args.star_size,
        group_imd=args.group_imd,
    )


if __name__ == "__main__":
    main()
