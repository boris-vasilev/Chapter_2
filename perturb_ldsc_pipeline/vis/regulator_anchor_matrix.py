#!/usr/bin/env python
"""Adjacency-matrix heatmap: H4 anchors × regulators for one (study, condition).

Rectangular heatmap of regulator → anchor relationships:

  - Rows = anchors (coloc H4 genes), from the gene-set file
  - Cols = regulators (perturbed cis genes), from the gene-set file
  - Cell = max |z-score| of the perturbation effect; blank if no edge
  - Both axes biclustered (average linkage on Euclidean distance) so
    anchors with shared regulators sit next to each other and reveal
    "modules" as visual blocks. --no-cluster falls back to alphabetical.

Node universes are read from data/gene_sets/<module>/{study}_{cond}_{annot}.txt
so the matrix reflects whatever K-filter the upstream module applied.
DE results CSV is the edge source only.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

CONDS = ["Rest", "Stim8hr", "Stim48hr"]


def read_gene_set(path):
    if not path.exists():
        raise SystemExit(f"Gene set file not found: {path}")
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def load_node_sets(gene_sets_dir, study, condition, anchor_annot, regulator_annot):
    a_path = Path(gene_sets_dir) / f"{study}_{condition}_{anchor_annot}.txt"
    r_path = Path(gene_sets_dir) / f"{study}_{condition}_{regulator_annot}.txt"
    anchors = read_gene_set(a_path)
    regs = read_gene_set(r_path)
    if not anchors:
        raise SystemExit(f"Anchor set empty: {a_path}")
    return anchors, regs, a_path, r_path


def load_edge_table(de_path, condition, anchors, regulators):
    return pl.read_csv(de_path).filter(
        pl.col("culture_condition") == condition,
        pl.col("cis_ensg").is_in(list(regulators)),
        pl.col("trans_ensg").is_in(list(anchors)),
        pl.col("cis_ensg") != pl.col("trans_ensg"),
    )


def build_matrix(df, anchors, regulators):
    """Returns (M, anchor_list, reg_list, symbols).
    M[i, j] = max |zscore| from regulator reg_list[j] to anchor anchor_list[i],
    0 if no edge.
    """
    anchor_list = sorted(anchors)
    reg_list = sorted(regulators)
    anchor_idx = {a: i for i, a in enumerate(anchor_list)}
    reg_idx = {r: j for j, r in enumerate(reg_list)}

    M = np.zeros((len(anchor_list), len(reg_list)))
    symbols = {}
    z_col = "zscore" if "zscore" in df.columns else None

    for r in df.iter_rows(named=True):
        cis, trans = r["cis_ensg"], r["trans_ensg"]
        i, j = anchor_idx[trans], reg_idx[cis]
        z = abs(float(r[z_col])) if z_col else 1.0
        if z > M[i, j]:
            M[i, j] = z
        if "cis_gene" in df.columns and r["cis_gene"]:
            symbols[cis] = r["cis_gene"]
        if "trans_gene" in df.columns and r["trans_gene"]:
            symbols[trans] = r["trans_gene"]

    return M, anchor_list, reg_list, symbols


def filter_cols(M, reg_list, min_hits):
    hits = (M > 0).sum(axis=0)
    keep = hits >= min_hits
    if not keep.any():
        raise SystemExit(f"No regulators hit ≥ {min_hits} anchors.")
    return M[:, keep], [r for r, k in zip(reg_list, keep) if k]


def filter_rows(M, anchor_list, min_regs):
    hits = (M > 0).sum(axis=1)
    keep = hits >= min_regs
    if not keep.any():
        raise SystemExit(f"No anchors have ≥ {min_regs} regulators.")
    return M[keep, :], [a for a, k in zip(anchor_list, keep) if k]


def cluster_order(matrix, axis=0, metric="jaccard", method="complete"):
    """Bicluster ordering for sparse connectivity matrices.

    Default = Jaccard distance + complete linkage on the binarised matrix.
    Jaccard ignores joint absences (which dominate sparse data and make
    Euclidean meaningless); complete linkage avoids chaining and produces
    tight, interpretable clusters.

    Pass metric='euclidean' + method='average' to recover the previous
    (worse) behaviour.
    """
    try:
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import pdist
    except ImportError:
        return np.arange(matrix.shape[axis])

    data = matrix if axis == 0 else matrix.T
    if data.shape[0] < 2:
        return np.arange(data.shape[0])

    if metric in ("jaccard", "hamming", "dice"):
        data = (data > 0).astype(int)

    # All-zero rows have undefined Jaccard distance — exclude them from
    # the clustering and append at the end.
    rowsum = data.sum(axis=1)
    nonzero = np.where(rowsum > 0)[0]
    zero    = np.where(rowsum == 0)[0]
    if len(nonzero) < 2:
        return np.arange(data.shape[0])

    sub = data[nonzero]
    d = pdist(sub, metric=metric)
    if not np.any(np.isfinite(d)) or not np.any(d):
        return np.arange(data.shape[0])

    Z = linkage(d, method=method)
    return np.concatenate([nonzero[leaves_list(Z)], zero])


def degree_order(matrix, axis=0):
    """Sort rows (or cols) by descending edge count — produces an L-shape."""
    data = matrix if axis == 0 else matrix.T
    return np.argsort(-(data > 0).sum(axis=1))


def leiden_membership(matrix, axis=0, resolution=1.0, seed=42):
    """Compute Leiden community labels on the co-occurrence projection.

    Returns a list of community IDs, one per row (or col) of the input
    matrix, or None if Leiden isn't usable (no edges, no library, etc).
    """
    data = matrix if axis == 0 else matrix.T
    binary = (data > 0).astype(int)
    n = binary.shape[0]
    if n < 3:
        return list(range(n))

    co = binary @ binary.T
    np.fill_diagonal(co, 0)
    edges, weights = [], []
    for i in range(n):
        for j in range(i + 1, n):
            w = int(co[i, j])
            if w > 0:
                edges.append((i, j))
                weights.append(w)

    if not edges:
        return list(range(n))
    return _run_leiden(n, edges, weights, resolution, seed)


def leiden_order(matrix, axis=0, resolution=1.0, seed=42):
    """Group rows (or cols) into Leiden communities and return ordering.

    Sort: same-community items grouped together, communities sorted by
    size descending, items within each community sorted by degree
    descending.
    """
    data = matrix if axis == 0 else matrix.T
    binary = (data > 0).astype(int)
    n = binary.shape[0]
    if n < 3:
        return np.arange(n)

    membership = leiden_membership(matrix, axis=axis,
                                   resolution=resolution, seed=seed)
    if membership is None:
        return cluster_order(matrix, axis=axis,
                             metric="jaccard", method="complete")

    deg = binary.sum(axis=1)
    by_comm = {}
    for i, c in enumerate(membership):
        by_comm.setdefault(c, []).append(i)
    ordered_comms = sorted(by_comm.values(), key=lambda xs: -len(xs))
    out = []
    for nodes in ordered_comms:
        out.extend(sorted(nodes, key=lambda i: -deg[i]))
    return np.array(out)


def write_edges_csv(M, anchor_list, reg_list, symbols,
                    anchor_membership, reg_membership,
                    study, condition, out_path):
    """One row per (regulator, anchor) edge, with Leiden communities.

    Columns: study, condition, anchor_ensg, anchor_symbol, regulator_ensg,
    regulator_symbol, abs_zscore, anchor_n_regulators, regulator_n_anchors,
    anchor_community, regulator_community.
    """
    anchor_deg = (M > 0).sum(axis=1)
    reg_deg    = (M > 0).sum(axis=0)
    rows = []
    for i, a in enumerate(anchor_list):
        for j, r in enumerate(reg_list):
            if M[i, j] <= 0:
                continue
            rows.append({
                "study": study,
                "condition": condition,
                "anchor_ensg": a,
                "anchor_symbol": symbols.get(a, a),
                "regulator_ensg": r,
                "regulator_symbol": symbols.get(r, r),
                "abs_zscore": float(M[i, j]),
                "anchor_n_regulators": int(anchor_deg[i]),
                "regulator_n_anchors": int(reg_deg[j]),
                "anchor_community": (
                    int(anchor_membership[i]) if anchor_membership else -1
                ),
                "regulator_community": (
                    int(reg_membership[j]) if reg_membership else -1
                ),
            })
    pl.DataFrame(rows).write_csv(out_path)
    print(f"Wrote {out_path}  ({len(rows)} edges)")


def _run_leiden(n_nodes, edges, weights, resolution, seed):
    """Try leidenalg first, fall back to networkx greedy modularity."""
    try:
        import igraph as ig
        import leidenalg
        g = ig.Graph(n=n_nodes, edges=edges, directed=False)
        g.es["weight"] = weights
        partition = leidenalg.find_partition(
            g, leidenalg.RBConfigurationVertexPartition,
            weights="weight",
            resolution_parameter=resolution,
            seed=seed,
        )
        return list(partition.membership)
    except ImportError:
        pass
    try:
        import networkx as nx
        from networkx.algorithms.community import greedy_modularity_communities
        G = nx.Graph()
        G.add_nodes_from(range(n_nodes))
        for (u, v), w in zip(edges, weights):
            G.add_edge(u, v, weight=w)
        comms = greedy_modularity_communities(G, weight="weight",
                                              resolution=resolution)
        membership = [0] * n_nodes
        for c_idx, comm in enumerate(comms):
            for node in comm:
                membership[node] = c_idx
        return membership
    except ImportError:
        return None


def plot(M, anchor_list, reg_list, symbols, study, condition, out_path,
         cmap, figsize, show_col_labels):
    n_rows, n_cols = M.shape
    fig, ax = plt.subplots(figsize=figsize)

    M_masked = np.ma.masked_where(M == 0, M)
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad("white")
    vmax = float(M.max()) if M.max() > 0 else 1.0
    im = ax.imshow(M_masked, cmap=cmap_obj, vmin=0, vmax=vmax, aspect="auto")

    row_labels = [symbols.get(a, a) for a in anchor_list]
    col_labels = [symbols.get(r, r) for r in reg_list]
    row_fs = max(6, min(11, 220 / max(n_rows, 1)))
    # Allow regulator labels down to 5pt so we can fit ~200 in a wide
    # figure rather than auto-hide them.
    col_fs = max(5, min(10, 280 / max(n_cols, 1)))

    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_labels, fontsize=row_fs)
    if show_col_labels:
        ax.set_xticks(range(n_cols))
        ax.set_xticklabels(col_labels, rotation=90, fontsize=col_fs)
    else:
        ax.set_xticks([])
    ax.tick_params(axis="both", length=0)
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)

    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="lightgrey", linewidth=0.25)
    ax.tick_params(which="minor", length=0)
    for s in ax.spines.values():
        s.set_linewidth(0.5)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("|z-score| of perturbation effect", fontsize=12)
    cbar.ax.tick_params(labelsize=10)

    trait = study.split("_")[0]
    n_edges = int((M > 0).sum())
    ax.set_title(
        f"{trait}  |  {condition}  —  H4 anchors × regulators\n"
        f"{n_rows} anchors, {n_cols} regulators, {n_edges} edges",
        fontsize=14, fontweight="bold", pad=10,
    )
    ax.set_xlabel("Regulators (perturbations)", fontsize=12)
    ax.set_ylabel("H4 anchors (coloc genes)", fontsize=12)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--de-results", required=True)
    parser.add_argument("--gene-sets-dir", required=True)
    parser.add_argument("--study", required=True)
    parser.add_argument("--condition", required=True, choices=CONDS)
    parser.add_argument("--anchor-annot", default="H4_genes")
    parser.add_argument("--regulator-annot", default="regulators_H4")
    parser.add_argument("--min-anchor-hits", type=int, default=1,
                        help="Drop regulators (columns) hitting fewer than "
                             "this many anchors (default 1).")
    parser.add_argument("--min-regulators", type=int, default=0,
                        help="Drop anchors (rows) with fewer than this many "
                             "regulators in the visible set (default 0).")
    parser.add_argument("--order",
                        choices=["leiden", "cluster", "degree", "alphabetical"],
                        default="leiden",
                        help="Row/column ordering. "
                             "'leiden' = community detection on the "
                             "co-occurrence projection (recommended; tries "
                             "leidenalg, falls back to networkx greedy "
                             "modularity); "
                             "'cluster' = Jaccard + complete-linkage "
                             "hierarchical; 'degree' = sort by edge count "
                             "descending (L-shape); 'alphabetical' = preserve "
                             "gene-symbol order.")
    parser.add_argument("--leiden-resolution", type=float, default=1.0,
                        help="Resolution parameter for Leiden — higher = "
                             "more, smaller communities (default 1.0).")
    parser.add_argument("--leiden-seed", type=int, default=42)
    parser.add_argument("--cluster-metric", default="jaccard",
                        choices=["jaccard", "hamming", "dice", "euclidean", "cosine"],
                        help="Distance metric for clustering "
                             "(default jaccard; ignored unless --order cluster).")
    parser.add_argument("--cluster-linkage", default="complete",
                        choices=["complete", "average", "single", "ward"],
                        help="Linkage method (default complete; ward requires "
                             "metric=euclidean).")
    parser.add_argument("--hide-col-labels", action="store_true",
                        help="Suppress regulator labels (useful when there "
                             "are too many to read).")
    parser.add_argument("--cmap", default="YlOrRd")
    parser.add_argument("--figsize", nargs=2, type=float, default=None,
                        metavar=("W", "H"),
                        help="Override auto-sized figure. Default scales "
                             "with matrix dimensions.")
    parser.add_argument("--out", required=True)
    parser.add_argument("--csv-out", default=None,
                        help="Optional path for a long-format edges CSV "
                             "(one row per regulator → anchor edge) with "
                             "Leiden community labels, gene symbols, and "
                             "degrees. Useful for handing off to another "
                             "analysis chat.")
    args = parser.parse_args()

    anchors, regulators, a_path, r_path = load_node_sets(
        args.gene_sets_dir, args.study, args.condition,
        args.anchor_annot, args.regulator_annot,
    )
    print(f"  Anchors:    {len(anchors)}   from {a_path}")
    print(f"  Regulators: {len(regulators)}   from {r_path}")

    df = load_edge_table(args.de_results, args.condition, anchors, regulators)
    M, anchor_list, reg_list, symbols = build_matrix(df, anchors, regulators)

    M, reg_list = filter_cols(M, reg_list, args.min_anchor_hits)
    if args.min_regulators > 0:
        M, anchor_list = filter_rows(M, anchor_list, args.min_regulators)
    print(f"  Matrix after filtering: {M.shape[0]} anchors × {M.shape[1]} regulators "
          f"({(M > 0).sum()} edges)")

    # Always compute Leiden membership in the *original* index space —
    # used by the CSV and (optionally) the plot ordering. This way the
    # CSV is meaningful even when --order is set to degree or alphabetical.
    anchor_membership = leiden_membership(
        M, axis=0,
        resolution=args.leiden_resolution, seed=args.leiden_seed,
    )
    reg_membership = leiden_membership(
        M, axis=1,
        resolution=args.leiden_resolution, seed=args.leiden_seed,
    )

    if args.csv_out:
        csv_path = Path(args.csv_out)
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        write_edges_csv(
            M, anchor_list, reg_list, symbols,
            anchor_membership, reg_membership,
            args.study, args.condition, csv_path,
        )

    if args.order == "leiden":
        row_order = leiden_order(M, axis=0,
                                 resolution=args.leiden_resolution,
                                 seed=args.leiden_seed)
        col_order = leiden_order(M, axis=1,
                                 resolution=args.leiden_resolution,
                                 seed=args.leiden_seed)
    elif args.order == "cluster":
        row_order = cluster_order(M, axis=0,
                                  metric=args.cluster_metric,
                                  method=args.cluster_linkage)
        col_order = cluster_order(M, axis=1,
                                  metric=args.cluster_metric,
                                  method=args.cluster_linkage)
    elif args.order == "degree":
        row_order = degree_order(M, axis=0)
        col_order = degree_order(M, axis=1)
    else:
        row_order = np.arange(M.shape[0])
        col_order = np.arange(M.shape[1])

    M = M[row_order, :][:, col_order]
    anchor_list = [anchor_list[i] for i in row_order]
    reg_list = [reg_list[j] for j in col_order]

    n_rows, n_cols = M.shape
    # Wider per-column scaling when many regulators, so labels stay legible.
    figsize = tuple(args.figsize) if args.figsize else (
        max(8, 0.22 * n_cols + 4),
        max(6, 0.30 * n_rows + 2),
    )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    show_cols = (not args.hide_col_labels) and n_cols <= 250
    plot(M, anchor_list, reg_list, symbols, args.study, args.condition,
         out_path, args.cmap, figsize, show_cols)


if __name__ == "__main__":
    main()
