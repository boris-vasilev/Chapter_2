#!/usr/bin/env python
"""Figure 2.1 - study design and anchors, as one file at print size (6.3 in).

  A. Schematic: chapter2/assets/ch2_schematic.png by default; --schematic PATH
     inserts another raster image (PNG/JPG; convert an SVG to PNG first) and
     --schematic none draws a light placeholder box. The panel is full width;
     its height follows the image's aspect ratio.
  B. Anchor coverage per disease and condition (06_anchor_coverage). The
     asterisk on disease names is explained in the caption, not the figure.
  C. Jaccard overlap of anchor sets between diseases (04_jaccard_heatmaps),
     rows named and numbered, columns numbered.
"""

import argparse
import importlib
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from common import (HERE, PRINT_W, PT, add_variant_args, axes_in,
                    load_diseases, load_edges, load_gene_sets, panel_label,
                    print_style, resolve_variant, save)

coverage = importlib.import_module("06_anchor_coverage")
jaccard = importlib.import_module("04_jaccard_heatmaps")

PLACEHOLDER = "Schematic: anchor → regulators / targets → S-LDSC, drug targets"
DEFAULT_SCHEMATIC = HERE.parent / "assets" / "ch2_schematic.png"
ROW_IN = 0.125          # height of one disease row in B and C (inches)


def load_schematic(path):
    if path is None or str(path).lower() == "none":
        return None
    path = Path(path)
    if path.suffix.lower() == ".svg":
        raise SystemExit("--schematic: convert the SVG to PNG first")
    if not path.exists():
        raise SystemExit(f"--schematic: {path} not found")
    return plt.imread(path)


def draw_schematic(ax, img=None):
    ax.set_axis_off()
    if img is not None:
        ax.imshow(img, interpolation="lanczos")
        return
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.add_patch(FancyBboxPatch((0.01, 0.06), 0.98, 0.88,
                                boxstyle="round,pad=0,rounding_size=0.02",
                                facecolor="#f6f6f4", edgecolor="#c9c8c2",
                                linestyle="--", linewidth=0.8))
    ax.text(0.5, 0.5, PLACEHOLDER, ha="center", va="center", fontsize=PT["title"],
            color="#6f6e69")


def draw_figure(img, sets, edges, diseases):
    W = PRINT_W
    n = len(diseases)
    side = n * ROW_IN
    a_h = W * img.shape[0] / img.shape[1] if img is not None else 1.6

    # Vertical layout, bottom-up (inches).
    c_y = 0.22                       # C matrix bottom (column numbers below)
    c_top = c_y + side + 0.08
    b_y = c_top + 0.75               # B axes bottom (x label + C title above C)
    b_top = b_y + side
    a_y = b_top + 0.68               # A image bottom (facet titles + legend)
    H = a_y + a_h + 0.22
    fig = plt.figure(figsize=(W, H))

    ax_a = axes_in(fig, 0, a_y, W, a_h)
    draw_schematic(ax_a, img)
    panel_label(fig, 0.02, H - 0.02, "A", "Study design")

    # B - anchor coverage facets.
    facet_w, gap, x0 = 1.42, 0.1, 1.68
    axes_b = [axes_in(fig, x0 + k * (facet_w + gap), b_y, facet_w, side)
              for k in range(3)]
    *_, handles = coverage.draw_anchor_coverage(axes_b, sets, edges, diseases,
                                                label_groups=False)
    for ax in axes_b:
        ax.set_xlabel("Anchors")
    fig.legend(handles=handles, loc="upper left", ncol=2,
               bbox_to_anchor=(x0 / W, (b_top + 0.5) / H), borderaxespad=0,
               handlelength=1.2, columnspacing=1.0)
    panel_label(fig, 0.02, b_top + 0.66, "B", "Anchor coverage in the perturbation graph")

    # C - anchor Jaccard heatmap.
    cx = 1.8
    ax_c = axes_in(fig, cx, c_y, side, side)
    cax = axes_in(fig, cx + side + 0.12, c_y, 0.08, side)
    jaccard.draw_anchor_jaccard_print(ax_c, cax, sets, diseases)
    fig.legend(handles=jaccard.group_handles(short=True), loc="center left",
               bbox_to_anchor=((cx + side + 0.72) / W, (c_y + side / 2) / H),
               borderaxespad=0, handlelength=1.2)
    panel_label(fig, 0.02, c_top + 0.3, "C", "Overlap of anchor sets between diseases")
    return fig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--schematic", default=str(DEFAULT_SCHEMATIC),
                    help="Raster image for panel A (PNG/JPG), or 'none' for the "
                         "placeholder. Default: assets/ch2_schematic.png.")
    add_variant_args(ap)
    args = resolve_variant(ap.parse_args())

    diseases = load_diseases()
    sets = load_gene_sets(args.gene_sets_dir)
    edges = load_edges(min_abs_logfc=args.min_abs_logfc)
    img = load_schematic(args.schematic)

    with print_style():
        fig = draw_figure(img, sets, edges, diseases)
    save(fig, "fig_2_1_overview", args.out_dir, exact=True)


if __name__ == "__main__":
    main()
