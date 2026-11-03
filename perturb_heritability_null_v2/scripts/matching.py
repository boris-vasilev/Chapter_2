"""Matched pseudo-anchor sampling for the coloc-anchor permutation null (v2).

Shared by make_permutation_gene_sets.py (runs in the Python 3.9 LDSC
container), the Stage 1 matching check, and the drug-target permutation copy,
so the matching is defined in exactly one place. polars + numpy only.

What changed from the pilot (perturb_ldsc_permutation/), see NOTES_null_v2.md:

  1. Stratum assignment. The pilot compared each real anchor's RAW out-degree /
     mean expression with the strata's BIN INDICES (0..n_bins-1), so every
     expressed anchor landed in the top expression bin. Here pool genes and
     anchors are binned by the same `Binning` object, and an anchor is sampled
     from its exact stratum.
  2. Out-degree: bin 0 = out-degree 0, then quantile bins among the non-zero
     pool genes (most pool genes are never knocked down, so plain quantiles
     collapse into one bin).
  3. In-degree (number of knockdowns that change the gene) is a third matching
     variable, binned the same way (zero bin + quantile bins among non-zero).
  4. Degrees come from the SAME edges that build the gene sets: one culture
     condition, |logFC| >= min_abs_logfc, of DE_results_strong.csv (already
     strong-perturbation filtered), self-edges dropped.
  5. All coloc genes (any trait, CD4, H4 > threshold) are excluded from the
     pool.

  6. Anchors with out-degree > 0 (knocked down in the screen) are matched by
     nearest neighbour on out-degree instead of an out-degree quartile: sample
     from the k (default 20) pool genes with out-degree > 0 closest in
     out-degree, among genes in the anchor's (in-degree bin, expression bin);
     if that cell has < k available genes, cells within L1 bin distance 1, 2,
     ... are added and the relaxation distance is logged. A quartile bin let a
     single knocked-down anchor be replaced by genes with very different
     out-degree, which made the target-set size null lumpy.

Anchors not measured in the screen (no expression, degree 0 by construction)
cannot be placed in a measured stratum; they are matched to genes that are
also not measured in the screen, are in the study's gene-coordinate file, are
not coloc genes, and share the anchor's coding status (protein-coding or not).
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np
import polars as pl


# --- edges and node properties ---------------------------------------------

def condition_edges(de: pl.DataFrame, cond: str, min_abs_logfc: float) -> pl.DataFrame:
    """Unique (cis_ensg, trans_ensg) edges used to build regulator/target sets
    for one condition: |logFC| >= min_abs_logfc, self-edges dropped."""
    e = de.filter(pl.col("culture_condition") == cond)
    if min_abs_logfc > 0:
        e = e.filter(pl.col("logFC").abs() >= min_abs_logfc)
    return (e.filter(pl.col("cis_ensg") != pl.col("trans_ensg"))
             .select("cis_ensg", "trans_ensg").unique())


def degrees(edges: pl.DataFrame):
    """(out_degree, in_degree) dicts. Out = distinct genes a knockdown changes;
    in = distinct knockdowns that change the gene."""
    out = edges.group_by("cis_ensg").agg(pl.len().alias("n"))
    inn = edges.group_by("trans_ensg").agg(pl.len().alias("n"))
    return (dict(zip(out["cis_ensg"].to_list(), out["n"].to_list())),
            dict(zip(inn["trans_ensg"].to_list(), inn["n"].to_list())))


def load_expression(path: str, cond: str) -> dict:
    """gene -> per-condition screen baseMean (from build_screen_expression.py).
    Genes with a missing value are omitted (treated as unmeasured)."""
    ex = (pl.read_csv(path, separator="\t")
            .filter(pl.col("condition") == cond)
            .drop_nulls("base_mean")
            .filter(pl.col("base_mean").is_not_nan()))
    return dict(zip(ex["gene"].to_list(), ex["base_mean"].to_list()))


def load_coloc_genes(coloc_path: str, cell_type: str, h4_threshold: float,
                     study: str | None = None) -> list:
    """Coloc H4 genes (version suffix stripped): one study, or all studies when
    study is None (used to exclude every coloc gene from the pool)."""
    q = (pl.scan_csv(coloc_path, separator="\t", infer_schema_length=10000)
           .filter(pl.col("cell_type") == cell_type)
           .filter(pl.col("H4") > h4_threshold))
    if study is not None:
        q = q.filter(pl.col("imd") == study)
    q = q.with_columns(pl.col("Gene_id").str.replace(r"\..*", "").alias("Gene_id"))
    return sorted(q.collect()["Gene_id"].drop_nulls().unique().to_list())


def read_gene_ids(coord_path: str) -> set:
    """GENE column of a space-separated gene-coordinate file."""
    return set(pl.read_csv(coord_path, separator=" ", columns=["GENE"])["GENE"].to_list())


# --- binning ----------------------------------------------------------------

def _zero_plus_quantile_edges(vals: np.ndarray, n_nonzero_bins: int) -> np.ndarray:
    nz = vals[vals > 0]
    if nz.size == 0:
        return np.array([])
    return np.quantile(nz, np.linspace(0, 1, n_nonzero_bins + 1))


def _bin_zero_plus_quantile(v: float, edges: np.ndarray, n_nonzero_bins: int) -> int:
    if v <= 0 or edges.size == 0:
        return 0
    i = int(np.searchsorted(edges[1:-1], v, side="right"))
    return 1 + max(0, min(n_nonzero_bins - 1, i))


def _bin_quantile(v: float, edges: np.ndarray, n_bins: int) -> int:
    i = int(np.searchsorted(edges[1:-1], v, side="right"))
    return max(0, min(n_bins - 1, i))


class Binning:
    """Bin edges fitted on the measured pool, applied identically to pool genes
    and real anchors. Stratum key = (out_bin, in_bin, expr_bin)."""

    def __init__(self, pool, out_deg, in_deg, expr,
                 n_degree_bins: int = 4, n_expr_bins: int = 5):
        self.out_deg, self.in_deg, self.expr = out_deg, in_deg, expr
        self.n_degree_bins, self.n_expr_bins = n_degree_bins, n_expr_bins
        o = np.array([out_deg.get(g, 0) for g in pool], dtype=float)
        i = np.array([in_deg.get(g, 0) for g in pool], dtype=float)
        e = np.array([expr[g] for g in pool], dtype=float)
        self.out_edges = _zero_plus_quantile_edges(o, n_degree_bins)
        self.in_edges = _zero_plus_quantile_edges(i, n_degree_bins)
        self.expr_edges = np.quantile(e, np.linspace(0, 1, n_expr_bins + 1))

    def values(self, g):
        return (self.out_deg.get(g, 0), self.in_deg.get(g, 0), self.expr[g])

    def key(self, g):
        o, i, e = self.values(g)
        return (_bin_zero_plus_quantile(o, self.out_edges, self.n_degree_bins),
                _bin_zero_plus_quantile(i, self.in_edges, self.n_degree_bins),
                _bin_quantile(e, self.expr_edges, self.n_expr_bins))


# --- the matcher ------------------------------------------------------------

class Matcher:
    """Pseudo-anchor sampler for one study x condition.

    measured_pool:   screen-measured genes with expression, minus coloc genes.
    unmeasured_pool: dict coding_flag -> genes in the study's gene-coordinate
                     file that are not measured and not coloc genes.
    """

    def __init__(self, measured_pool, unmeasured_pool, protein_coding,
                 out_deg, in_deg, expr, n_degree_bins=4, n_expr_bins=5,
                 knn_k=20):
        self.measured_pool = sorted(measured_pool)
        self.unmeasured_pool = {k: sorted(v) for k, v in unmeasured_pool.items()}
        self.protein_coding = protein_coding
        self.expr = expr
        self.knn_k = knn_k
        self.binning = Binning(self.measured_pool, out_deg, in_deg, expr,
                               n_degree_bins, n_expr_bins)
        strata = defaultdict(list)
        # Knocked-down pool genes grouped by (in_bin, expr_bin) for k-NN.
        knock = defaultdict(list)
        for g in self.measured_pool:
            k = self.binning.key(g)
            strata[k].append(g)
            if k[0] > 0:
                knock[k[1:]].append(g)
        self.strata = dict(strata)
        self.knock_cells = dict(knock)

    def is_measured(self, g):
        return g in self.expr

    def _knn_candidates(self, g, drawn):
        """The k nearest available pool genes by out-degree (ties at the k-th
        distance kept), widening the (in_bin, expr_bin) cell set by L1 bin
        distance until >= k genes are available. Returns (cands, relax_dist)."""
        _, ib, eb = self.binning.key(g)
        o = self.binning.out_deg.get(g, 0)
        max_d = self.binning.n_degree_bins + self.binning.n_expr_bins
        avail, d = [], 0
        for d in range(max_d + 1):
            avail = [x for (ci, ce), genes in self.knock_cells.items()
                     if abs(ci - ib) + abs(ce - eb) <= d
                     for x in genes if x not in drawn]
            if len(avail) >= self.knn_k:
                break
        if not avail:
            return [], d
        dist = sorted((abs(self.binning.out_deg.get(x, 0) - o), x) for x in avail)
        cutoff = dist[min(self.knn_k, len(dist)) - 1][0]
        return sorted(x for dd, x in dist if dd <= cutoff), d

    def sample(self, real_anchors, rng):
        """One pseudo-anchor per real anchor (sorted iteration for
        reproducibility), without replacement. Returns (pseudo_set, log) where
        log counts exact-stratum draws, nearest-stratum fallbacks (with their
        bin-index L1 distances), k-NN draws for knocked-down anchors (with
        their cell relaxation distances), and unmeasured-anchor draws."""
        drawn = set()
        log = {"n_exact": 0, "n_fallback": 0, "fallback_dist": [],
               "n_knn": 0, "n_knn_relaxed": 0, "knn_relax_dist": [],
               "n_unmeasured": 0}
        for g in sorted(real_anchors):
            if not self.is_measured(g):
                cands = [x for x in self.unmeasured_pool[g in self.protein_coding]
                         if x not in drawn]
                log["n_unmeasured"] += 1
            elif self.binning.out_deg.get(g, 0) > 0:
                cands, d = self._knn_candidates(g, drawn)
                log["n_knn"] += 1
                log["knn_relax_dist"].append(d)
                if d > 0:
                    log["n_knn_relaxed"] += 1
            else:
                k = self.binning.key(g)
                cands = [x for x in self.strata.get(k, []) if x not in drawn]
                if cands:
                    log["n_exact"] += 1
                else:
                    # Nearest non-empty stratum by bin-index (L1) distance; all
                    # strata at the minimum distance are pooled.
                    best, best_d = [], None
                    for sk, genes in self.strata.items():
                        avail = [x for x in genes if x not in drawn]
                        if not avail:
                            continue
                        d = sum(abs(a - b) for a, b in zip(sk, k))
                        if best_d is None or d < best_d:
                            best, best_d = avail, d
                        elif d == best_d:
                            best = best + avail
                    cands = sorted(best)
                    log["n_fallback"] += 1
                    log["fallback_dist"].append(best_d)
            if cands:
                drawn.add(cands[int(rng.integers(len(cands)))])
        return drawn, log


def perm_rng(base_seed: int, perm_idx: int, cond: str):
    """Independent, run-order-free stream per (perm, condition): the pilot's
    seed (base_seed + perm) plus a stable per-condition component."""
    cond_code = sum((i + 1) * ord(c) for i, c in enumerate(cond))
    return np.random.default_rng([base_seed + perm_idx, cond_code])


def build_matcher(de, cond, min_abs_logfc, expression_path, all_coloc,
                  gene_coord_path, protein_coding_coord_path,
                  n_degree_bins=4, n_expr_bins=5, knn_k=20):
    """Assemble edges, degrees, expression and pools for one study x condition.
    Returns (matcher, edges)."""
    edges = condition_edges(de, cond, min_abs_logfc)
    out_deg, in_deg = degrees(edges)
    expr = load_expression(expression_path, cond)
    coloc = set(all_coloc)
    measured_pool = [g for g in expr if g not in coloc]
    coord = read_gene_ids(gene_coord_path)
    pc = read_gene_ids(protein_coding_coord_path)
    unmeasured = [g for g in coord if g not in expr and g not in coloc]
    unmeasured_pool = {True: [g for g in unmeasured if g in pc],
                       False: [g for g in unmeasured if g not in pc]}
    m = Matcher(measured_pool, unmeasured_pool, pc, out_deg, in_deg, expr,
                n_degree_bins, n_expr_bins, knn_k)
    return m, edges


# --- gene-set construction (same convergence rules as make_gene_sets_coloc) --

def convergent_regulators(edges: pl.DataFrame, anchors, min_k: int = 1) -> list:
    """Cis genes whose perturbation changes >= min_k distinct anchors."""
    return (edges.filter(pl.col("trans_ensg").is_in(list(anchors)))
                 .group_by("cis_ensg").agg(pl.col("trans_ensg").n_unique().alias("n"))
                 .filter(pl.col("n") >= min_k)["cis_ensg"].unique().to_list())


def convergent_targets(edges: pl.DataFrame, anchors, min_k: int = 1) -> list:
    """Trans genes hit by >= min_k distinct anchor perturbations."""
    return (edges.filter(pl.col("cis_ensg").is_in(list(anchors)))
                 .group_by("trans_ensg").agg(pl.col("cis_ensg").n_unique().alias("n"))
                 .filter(pl.col("n") >= min_k)["trans_ensg"].unique().to_list())
