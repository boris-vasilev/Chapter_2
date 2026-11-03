#!/usr/bin/env python
"""Bipartite network: H4 anchors ↔ regulators for one (study, condition).

Node universe is taken from the gene-set files under data/gene_sets/<module>/
(so the network reflects whatever K-filter was applied at make_gene_sets
time — e.g. coloc K=1, coloc_strict K=3). The DE table is used only to
recover which regulator → anchor edges to draw.

  - Anchors (coloc H4 genes) are pinned to a circle, large red nodes.
  - Regulators (cis perturbations whose KD significantly changes an
    anchor's expression) are blue nodes around them; shade darkens
    with the number of distinct anchors they hit (a regulator hitting
    many anchors is a convergent / hub regulator).
  - Edges connect a regulator to every anchor it modulates.
  - Spring layout, with anchors held fixed — regulators settle near
    the anchor(s) they connect to.

Use it to ask: which anchors share regulators? Which regulators are
hubs touching multiple anchors at once?
"""

import argparse
from collections import defaultdict
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import polars as pl

CONDS = ["Rest", "Stim8hr", "Stim48hr"]


def read_gene_set(path):
    """Read a one-ENSG-per-line gene set file into a set."""
    if not path.exists():
        raise SystemExit(f"Gene set file not found: {path}")
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def load_node_sets(gene_sets_dir, study, condition, anchor_annot, regulator_annot):
    """Read anchors + regulators from the gene-set files."""
    anchors_path = Path(gene_sets_dir) / f"{study}_{condition}_{anchor_annot}.txt"
    regs_path    = Path(gene_sets_dir) / f"{study}_{condition}_{regulator_annot}.txt"
    anchors = read_gene_set(anchors_path)
    regs    = read_gene_set(regs_path)
    if not anchors:
        raise SystemExit(f"Anchor set empty: {anchors_path}")
    return anchors, regs, anchors_path, regs_path


def load_edges(de_path, condition, anchors, regulators):
    """Select DE rows whose (cis, trans) pair is a regulator → anchor edge.

    Returns (edges, ensg_to_symbol). Edges are filtered to only those
    where the regulator is in our regulator universe AND the anchor is
    in our anchor universe — i.e. the same set of relationships that
    drove the gene-set selection upstream.
    """
    df = pl.read_csv(de_path).filter(
        pl.col("culture_condition") == condition,
        pl.col("cis_ensg").is_in(list(regulators)),
        pl.col("trans_ensg").is_in(list(anchors)),
        pl.col("cis_ensg") != pl.col("trans_ensg"),
    )
    edges = []
    symbols = {}
    z_col = "zscore" if "zscore" in df.columns else None
    for r in df.iter_rows(named=True):
        z = float(r[z_col]) if z_col else 0.0
        edges.append((r["cis_ensg"], r["trans_ensg"], z))
        if "cis_gene" in df.columns and r["cis_gene"]:
            symbols[r["cis_ensg"]] = r["cis_gene"]
        if "trans_gene" in df.columns and r["trans_gene"]:
            symbols[r["trans_ensg"]] = r["trans_gene"]
    return edges, symbols


def build_graph(anchors, regulators, edges, min_anchor_hits):
    hits = defaultdict(set)
    for cis, trans, _ in edges:
        hits[cis].add(trans)
    keep_regs = {
        cis for cis in regulators
        if len(hits[cis]) >= min_anchor_hits
    }

    G = nx.Graph()
    connected_anchors = set()
    for cis, trans, z in edges:
        if cis in keep_regs and trans in anchors:
            G.add_edge(cis, trans, zscore=abs(z))
            connected_anchors.add(trans)

    for a in anchors:
        G.add_node(a, kind="anchor")
    for cis in keep_regs:
        G.add_node(cis, kind="regulator", n_anchors=len(hits[cis]))

    return G, hits, connected_anchors


def compute_layout(G, anchors, seed):
    anchor_list = sorted(a for a in anchors if a in G.nodes)
    pos = {}
    n = len(anchor_list)
    if n > 0:
        radius = 1.0
        for i, a in enumerate(anchor_list):
            angle = 2 * np.pi * i / n
            pos[a] = (radius * np.cos(angle), radius * np.sin(angle))

    rng = np.random.default_rng(seed)
    for node in G.nodes:
        if node not in pos:
            pos[node] = (rng.normal(0, 0.1), rng.normal(0, 0.1))

    if len(G.nodes) > 0:
        pos = nx.spring_layout(
            G, pos=pos, fixed=anchor_list,
            k=1.6 / np.sqrt(max(len(G.nodes), 1)),
            iterations=200, seed=seed,
        )
    return pos


def plot(G, pos, symbols, study, condition, out_path,
         label_min_hits, figsize, anchor_label_size, reg_label_size):
    fig, ax = plt.subplots(figsize=figsize)

    anchors = [n for n in G.nodes if G.nodes[n]["kind"] == "anchor"]
    regs = [n for n in G.nodes if G.nodes[n]["kind"] == "regulator"]
    isolated_anchors = [a for a in anchors if G.degree(a) == 0]

    if G.number_of_edges() > 0:
        edge_widths = []
        for _, _, d in G.edges(data=True):
            z = d.get("zscore", 1.0)
            edge_widths.append(0.3 + min(z * 0.15, 2.0))
        nx.draw_networkx_edges(
            G, pos, ax=ax, width=edge_widths,
            edge_color="grey", alpha=0.35,
        )

    if regs:
        reg_hits = [G.nodes[r]["n_anchors"] for r in regs]
        reg_sizes = [80 + 50 * h for h in reg_hits]
        sc = nx.draw_networkx_nodes(
            G, pos, nodelist=regs, ax=ax,
            node_color=reg_hits, cmap="Blues",
            vmin=1, vmax=max(max(reg_hits), 2),
            node_size=reg_sizes,
            edgecolors="white", linewidths=0.6,
        )
        cbar = fig.colorbar(sc, ax=ax, fraction=0.025, pad=0.02, shrink=0.4)
        cbar.set_label("Anchors hit by regulator", fontsize=10)
        cbar.ax.tick_params(labelsize=9)

    nx.draw_networkx_nodes(
        G, pos, nodelist=anchors, ax=ax,
        node_color="#e74c3c", node_size=550,
        edgecolors="black", linewidths=1.2,
    )

    anchor_labels = {a: symbols.get(a, a) for a in anchors}
    nx.draw_networkx_labels(
        G, pos, labels=anchor_labels, ax=ax,
        font_size=anchor_label_size, font_weight="bold",
    )

    reg_labels = {
        r: symbols.get(r, r) for r in regs
        if G.nodes[r]["n_anchors"] >= label_min_hits
    }
    if reg_labels:
        nx.draw_networkx_labels(
            G, pos, labels=reg_labels, ax=ax,
            font_size=reg_label_size, font_color="#222222",
        )

    trait = study.split("_")[0]
    title = (
        f"{trait}  |  {condition}  —  H4 anchors and their regulators\n"
        f"{len(anchors)} anchors ({len(isolated_anchors)} with no regulator), "
        f"{len(regs)} regulators shown, {G.number_of_edges()} edges"
    )
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)

    legend_handles = [
        mpatches.Patch(color="#e74c3c", label="H4 anchor (coloc gene)"),
        mpatches.Patch(color="#6baed6", label="Regulator (perturbation)"),
    ]
    ax.legend(
        handles=legend_handles, loc="upper left",
        bbox_to_anchor=(0.0, 1.0), frameon=False, fontsize=11,
    )

    ax.set_axis_off()
    pad = 0.18
    xs = [p[0] for p in pos.values()]
    ys = [p[1] for p in pos.values()]
    if xs and ys:
        ax.set_xlim(min(xs) - pad, max(xs) + pad)
        ax.set_ylim(min(ys) - pad, max(ys) + pad)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--de-results", required=True,
                        help="DE results CSV (cis_ensg/trans_ensg columns). "
                             "Used only to determine which regulator → anchor "
                             "edges to draw.")
    parser.add_argument("--gene-sets-dir", required=True,
                        help="Directory with {study}_{cond}_{annot}.txt files "
                             "(e.g. data/gene_sets/coloc or "
                             "data/gene_sets/coloc_strict). The anchor and "
                             "regulator node sets are read from here, so the "
                             "network reflects whatever K-filter the upstream "
                             "module applied.")
    parser.add_argument("--study", required=True,
                        help="Study identifier matching the gene-set "
                             "filenames (e.g. ASTAO_Ferreira_30929738_1).")
    parser.add_argument("--condition", required=True, choices=CONDS)
    parser.add_argument("--anchor-annot", default="H4_genes",
                        help="Annotation name for anchors in the gene-set "
                             "filenames (default H4_genes; use l2g_genes for "
                             "the L2G module).")
    parser.add_argument("--regulator-annot", default="regulators_H4",
                        help="Annotation name for regulators (default "
                             "regulators_H4; use regulators_l2g for L2G).")
    parser.add_argument("--min-anchor-hits", type=int, default=1,
                        help="Only draw regulators hitting at least this many "
                             "anchors (default 1 = all regulators in the "
                             "gene-set file).")
    parser.add_argument("--label-min-hits", type=int, default=2,
                        help="Only label regulators hitting at least this many "
                             "anchors (default 2 = hub regulators only).")
    parser.add_argument("--figsize", nargs=2, type=float, default=(13.0, 13.0),
                        metavar=("W", "H"))
    parser.add_argument("--anchor-label-size", type=float, default=11.0)
    parser.add_argument("--reg-label-size", type=float, default=8.5)
    parser.add_argument("--seed", type=int, default=42,
                        help="Layout seed for reproducibility.")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    anchors, regulators, anchors_path, regs_path = load_node_sets(
        args.gene_sets_dir, args.study, args.condition,
        args.anchor_annot, args.regulator_annot,
    )
    print(f"  Anchors:    {len(anchors)}   from {anchors_path}")
    print(f"  Regulators: {len(regulators)}   from {regs_path}")

    edges, symbols = load_edges(
        args.de_results, args.condition, anchors, regulators
    )

    G, hits, connected_anchors = build_graph(
        anchors, regulators, edges, args.min_anchor_hits
    )

    n_regs = sum(1 for n in G.nodes if G.nodes[n].get("kind") == "regulator")
    print(f"  In-graph anchors with ≥1 regulator: {len(connected_anchors)}")
    print(f"  In-graph regulators (≥ {args.min_anchor_hits} anchor hits): {n_regs}")
    print(f"  Edges: {G.number_of_edges()}")

    pos = compute_layout(G, anchors, args.seed)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plot(
        G, pos, symbols, args.study, args.condition, out_path,
        args.label_min_hits, tuple(args.figsize),
        args.anchor_label_size, args.reg_label_size,
    )


if __name__ == "__main__":
    main()
