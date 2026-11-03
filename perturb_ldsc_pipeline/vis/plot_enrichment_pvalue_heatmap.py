#!/usr/bin/env python
"""Side-by-side log2(Enrichment) / -log10(p) heatmap with FDR stars.

Reads the long-format CSV from parse_ldsc_results.py and renders one
panel showing log2(Enrichment) (diverging) and one showing -log10(p)
(sequential), with:

  - rows  = {trait label} | {condition}
  - cols  = the annotations passed via --annots
  - stars = BH-FDR < --fdr (default 0.05), computed across all visible cells
  - grey hatched cells = annotation absent for that (trait, condition)
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import polars as pl

COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]


def bh_fdr(pvals):
    """Benjamini-Hochberg adjusted p-values. Ignores NaN entries."""
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
    if "traitFromSource" in df.columns and df["traitFromSource"].null_count() < df.height:
        return "traitFromSource"
    if "trait_short" in df.columns:
        return "trait_short"
    return "study"


def build_matrices(df, label_col, annots, log2_floor):
    """Return enrichment, -log10(p), fdr, and present-mask 2-D arrays.

    Rows are sorted alphabetically by label, then by COND_ORDER.
    Columns follow the order in `annots`.
    Aggregates duplicate (label, condition, annot) cells by mean.
    """
    agg = df.group_by([label_col, "condition", "annot"]).agg(
        [pl.col("enrichment").mean(), pl.col("enrichment_p").mean()]
    )

    labels = sorted(agg[label_col].unique().to_list())
    row_ids = [(lbl, cond) for lbl in labels for cond in COND_ORDER]

    n_rows, n_cols = len(row_ids), len(annots)
    enr = np.full((n_rows, n_cols), np.nan)
    pval = np.full((n_rows, n_cols), np.nan)

    lookup = {
        (r[label_col], r["condition"], r["annot"]): (r["enrichment"], r["enrichment_p"])
        for r in agg.iter_rows(named=True)
    }
    for i, (lbl, cond) in enumerate(row_ids):
        for j, annot in enumerate(annots):
            hit = lookup.get((lbl, cond, annot))
            if hit is None:
                continue
            enr[i, j], pval[i, j] = hit

    log2_enr = np.where(np.isnan(enr), np.nan, np.log2(np.clip(enr, log2_floor, None)))
    neg_log_p = -np.log10(np.clip(pval, 1e-300, None))
    fdr = bh_fdr(pval.ravel()).reshape(pval.shape)
    present = ~np.isnan(enr)

    return row_ids, log2_enr, neg_log_p, fdr, present


def _draw_heatmap(ax, data, present, fdr, fdr_thr, cmap, vmin, vmax, cbar_label,
                  col_labels, row_labels, sep_after, show_row_labels):
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad("white")

    plot = np.where(present, data, np.nan)
    im = ax.imshow(plot, cmap=cmap_obj, vmin=vmin, vmax=vmax, aspect="auto")

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if not present[i, j]:
                ax.add_patch(
                    mpatches.Rectangle(
                        (j - 0.5, i - 0.5), 1, 1,
                        facecolor="lightgrey", edgecolor="white",
                        hatch="///", linewidth=0,
                    )
                )

    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if not present[i, j]:
                continue
            if not np.isnan(fdr[i, j]) and fdr[i, j] < fdr_thr:
                ax.text(j, i, "★", ha="center", va="center",
                        fontsize=8, color="black")

    if sep_after is not None and 0 < sep_after < data.shape[1]:
        ax.axvline(sep_after - 0.5, color="black", linewidth=1.5)

    ax.set_xticks(range(data.shape[1]))
    ax.set_xticklabels(col_labels, rotation=45, ha="right")
    ax.set_yticks(range(data.shape[0]))
    if show_row_labels:
        ax.set_yticklabels(row_labels, fontsize=7)
    else:
        ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(-0.5, data.shape[1] - 0.5)
    ax.set_ylim(data.shape[0] - 0.5, -0.5)

    cbar = plt.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
    cbar.set_label(cbar_label)
    return im


def plot(df, out_path, label_col, annots, display_labels, sep_after, fdr_thr,
         log2_vmin, log2_vmax, p_vmax, log2_cmap):
    # Anything below log2_vmin saturates at the floor of the colour scale.
    # With log2_vmin = 0 (default) this hides depletion (Enrichment < 1)
    # at the lightest shade; the right panel still flags it via p-value.
    log2_floor = 2 ** log2_vmin
    row_ids, log2_enr, neg_log_p, fdr, present = build_matrices(
        df, label_col, annots, log2_floor
    )
    row_labels = [f"{lbl} | {cond}" for lbl, cond in row_ids]
    col_labels = display_labels or annots

    if p_vmax is None:
        p_vmax = float(np.nanmax(neg_log_p[present])) if present.any() else 1.0

    n_rows = len(row_labels)
    fig_h = max(6, 0.28 * n_rows + 2)
    fig, axes = plt.subplots(
        1, 2, figsize=(11, fig_h),
        gridspec_kw={"wspace": 0.28},
    )

    _draw_heatmap(
        axes[0], log2_enr, present, fdr, fdr_thr,
        cmap=log2_cmap, vmin=log2_vmin, vmax=log2_vmax,
        cbar_label="log2(Enrichment)",
        col_labels=col_labels, row_labels=row_labels,
        sep_after=sep_after, show_row_labels=True,
    )
    axes[0].set_title(f"log2(Enrichment)   ★ = FDR < {fdr_thr}", fontsize=10)

    _draw_heatmap(
        axes[1], neg_log_p, present, fdr, fdr_thr,
        cmap="YlOrRd", vmin=0, vmax=p_vmax,
        cbar_label="-log10(p)",
        col_labels=col_labels, row_labels=row_labels,
        sep_after=sep_after, show_row_labels=False,
    )
    axes[1].set_title(f"-log10(p)   ★ = FDR < {fdr_thr}", fontsize=10)

    fig.suptitle(
        "LDSC Enrichment -- Custom Gene-Set Annotations\n"
        "(grey hatching = annotation absent from this study)",
        fontsize=12, fontweight="bold",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", required=True,
                        help="CSV from parse_ldsc_results.py")
    parser.add_argument("--annots", nargs="+", required=True,
                        help="Annotation column order (left → right). "
                             "E.g. H4_genes targets_H4 regulators_H4")
    parser.add_argument("--display-labels", nargs="+", default=None,
                        help="Pretty names for --annots (same length).")
    parser.add_argument("--label-col", default=None,
                        help="Row label column (default: traitFromSource → "
                             "trait_short → study).")
    parser.add_argument("--sep-after", type=int, default=None,
                        help="Draw a vertical separator after this many "
                             "columns from the left.")
    parser.add_argument("--fdr", type=float, default=0.05,
                        help="FDR threshold for stars (default 0.05).")
    parser.add_argument("--log2-vmin", type=float, default=0.0,
                        help="Floor for log2(Enrichment) colour scale. "
                             "Default 0 hides depletion at the lightest shade.")
    parser.add_argument("--log2-vmax", type=float, default=6.0)
    parser.add_argument("--log2-cmap", default="Reds",
                        help="Matplotlib colormap for log2(Enrichment). "
                             "Use a sequential map ('Reds', 'YlOrRd') with "
                             "log2-vmin=0, or 'RdBu_r' with a negative "
                             "log2-vmin to show depletion in blue.")
    parser.add_argument("--p-vmax", type=float, default=None,
                        help="Cap for the -log10(p) colour scale "
                             "(default: data max).")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    if args.display_labels and len(args.display_labels) != len(args.annots):
        parser.error("--display-labels must match --annots in length")

    df = pl.read_csv(args.summary)
    df = df.filter(pl.col("annot").is_in(args.annots))
    if df.height == 0:
        raise SystemExit(f"No rows match --annots {args.annots}")

    label_col = pick_label_col(df, args.label_col)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    plot(
        df, out_path,
        label_col=label_col,
        annots=args.annots,
        display_labels=args.display_labels,
        sep_after=args.sep_after,
        fdr_thr=args.fdr,
        log2_vmin=args.log2_vmin,
        log2_vmax=args.log2_vmax,
        p_vmax=args.p_vmax,
        log2_cmap=args.log2_cmap,
    )


if __name__ == "__main__":
    main()
