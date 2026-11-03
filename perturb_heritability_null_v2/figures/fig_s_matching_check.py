#!/usr/bin/env python
"""Supplementary: does the matched null match the real anchors? (Rest)

Rows = AITD, ASTCO. Columns:
  expression ECDF and in-degree ECDF of anchors (real, pilot, k-NN null);
  the knocked-down anchor's out-degree against the knocked-down pool genes in
  its (in-degree bin, expression bin) cell, with its k-NN neighbours marked;
  regulator and target set sizes (pilot, quartile-fixed, k-NN null; real).

Pilot and k-NN draws: source_data/ from scripts/matching_check_knn.py.
Quartile-fixed set sizes: the pre-k-NN matching check
(data/permutations/matching_check/matching_check_draws.csv, method "fixed").

  uv run python figures/fig_s_matching_check.py
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from matplotlib.lines import Line2D

from fig_common import (
    DISEASE_NAMES, FIG_WIDTH, MATCH, SET_COLOURS, SOURCE, plain_log_axis, save,
    setup_style,
)

METHOD_COLOURS = {"real": "#000000", "old": "#E69F00", "quartile": "#56B4E9",
                  "knn": "#0072B2"}
METHOD_LABELS = {"real": "Real", "old": "Pilot matching",
                 "quartile": "Corrected, quartile bins", "knn": "Corrected, k-NN"}


def ecdf(ax, x, colour, lw=1.0, ls="-"):
    x = np.sort(np.asarray(x, dtype=float))
    ax.step(x, np.arange(1, x.size + 1) / x.size, where="post", color=colour,
            lw=lw, ls=ls)


def load():
    vals = pl.read_csv(SOURCE / "matching_check_knn_anchor_values.csv")
    knn = pl.read_csv(SOURCE / "matching_check_knn_draws.csv")
    quart = (pl.read_csv(MATCH / "matching_check_draws.csv")
             .filter(pl.col("method") == "fixed")
             .with_columns(pl.lit("quartile").alias("method"))
             .select(knn.columns))
    neigh = pl.read_csv(SOURCE / "knn_neighbours.csv")
    return vals, pl.concat([knn, quart]), neigh


def anchor_ecdf(ax, vals, tr, var):
    for m in ("old", "knn", "real"):
        x = (vals.filter((pl.col("study") == tr) & (pl.col("method") == m))[var]
             .drop_nulls().to_numpy().astype(float))
        ecdf(ax, x + 1, METHOD_COLOURS[m], lw=1.4 if m == "real" else 1.0,
             ls="--" if m == "real" else "-")
    plain_log_axis(ax, "x", [1, 10, 100] if var == "in_degree"
                   else [1, 100, 10000])
    ax.set_ylim(0, 1.02)


def neighbour_panel(ax, neigh, tr):
    sub = neigh.filter(pl.col("study") == tr)
    a_out = sub["anchor_out_degree"].item(0)
    rng = np.random.default_rng(0)
    cell = sub.filter(pl.col("in_cell"))
    nb = sub.filter(pl.col("is_neighbour"))
    y = rng.uniform(-0.3, 0.3, cell.height)
    ax.scatter(cell["out_degree"], y, s=7, color="#c8c8c8", lw=0, zorder=1)
    yn = rng.uniform(-0.3, 0.3, nb.height)
    ax.scatter(nb["out_degree"], yn, s=9, facecolor="none",
               edgecolor=METHOD_COLOURS["knn"], lw=0.8, zorder=2)
    ax.axvline(a_out, ymax=0.62, color=SET_COLOURS["anchors"], lw=1.3, zorder=3)
    below = int((nb["out_degree"] < a_out).sum())
    ax.text(0.02, 0.98, f"anchor: {a_out}\ncell: {cell.height} genes\n"
            f"k-NN: {nb.height}, {below} below", transform=ax.transAxes, va="top",
            fontsize=6, linespacing=1.2)
    plain_log_axis(ax, "x", [1, 30, 1000])
    ax.set_xlim(1, max(cell["out_degree"].max(), a_out) * 1.6)
    ax.set_ylim(-0.6, 1.4)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)


def size_panel(ax, draws, tr, var):
    sub = draws.filter(pl.col("study") == tr)
    real = sub.filter(pl.col("method") == "real")[var].item(0)
    null = sub.filter(pl.col("method") != "real")
    x_all = null[var].to_numpy().astype(float)
    pos = x_all[x_all > 0]
    bins = np.geomspace(max(1, pos.min()), max(pos.max(), real) * 1.05, 31)
    notes = []
    for m in ("old", "quartile", "knn"):
        x = null.filter(pl.col("method") == m)[var].to_numpy().astype(float)
        ax.hist(x[x > 0], bins=bins, histtype="step", color=METHOD_COLOURS[m], lw=1.0)
        if (x == 0).any():
            notes.append(f"{int((x == 0).sum())} empty\n({m if m != 'old' else 'pilot'})")
    ax.axvline(real, color="k", lw=1.2, ls="--")
    ticks = [10, 1000] if var == "n_targets" else [100, 200, 300, 500]
    plain_log_axis(ax, "x", ticks)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.3)
    if notes:
        ax.text(0.02, 0.98, "\n".join(notes), transform=ax.transAxes, ha="left",
                va="top", fontsize=6)


def main():
    setup_style()
    vals, draws, neigh = load()
    traits = ["AITD", "ASTCO"]
    titles = ["Expression", "In-degree", "Knocked-down anchor", "Regulator set size",
              "Target set size"]
    xlabels = ["Screen baseMean + 1", "In-degree + 1", "Out-degree",
               "Regulators", "Targets"]
    fig, axes = plt.subplots(2, 5, figsize=(FIG_WIDTH, 3.9), squeeze=False)
    for i, tr in enumerate(traits):
        anchor_ecdf(axes[i][0], vals, tr, "expression")
        anchor_ecdf(axes[i][1], vals, tr, "in_degree")
        neighbour_panel(axes[i][2], neigh, tr)
        size_panel(axes[i][3], draws, tr, "n_regulators")
        size_panel(axes[i][4], draws, tr, "n_targets")
        axes[i][0].set_ylabel(f"{DISEASE_NAMES[tr]}\n\nECDF of anchors".replace(
            "thyroid disease", "thyroid\ndisease").replace(" (childhood", "\n(childhood"))
        axes[i][3].set_ylabel("Draws")
        for j in range(5):
            if i == 0:
                axes[i][j].set_title(titles[j], fontsize=8)
            axes[i][j].set_xlabel(xlabels[j], fontsize=7)
    handles = [Line2D([], [], color=METHOD_COLOURS[m], lw=1.2,
                      ls="--" if m == "real" else "-", label=METHOD_LABELS[m])
               for m in ("real", "old", "quartile", "knn")]
    handles += [Line2D([], [], color=SET_COLOURS["anchors"], lw=1.3,
                       label="Knocked-down anchor"),
                Line2D([], [], marker="o", ls="", ms=3, color="#c8c8c8",
                       label="Pool genes in cell"),
                Line2D([], [], marker="o", ls="", ms=3, mfc="none",
                       mec=METHOD_COLOURS["knn"], label="k-NN neighbours (k = 20)")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=6.5, bbox_to_anchor=(0.5, -0.07))
    fig.tight_layout(w_pad=0.9, h_pad=1.0)
    save(fig, "fig_s_matching_check")


if __name__ == "__main__":
    main()
