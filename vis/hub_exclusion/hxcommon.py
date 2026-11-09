"""Shared setup for the hub-exclusion comparison (figures/hub_exclusion/).

Four variants of the deg0 design: no exclusion (deg0) and the top 1, 5 and
10 % of regulators by out-degree (n_total_de_genes) dropped per condition
(hub0.99 / hub0.95 / hub0.90). deg75 is not part of this comparison.

Hub list. The pipeline (nf_ldsc_pipeline/bin/extract_de_results.py) drops
obs with n_total_de_genes > np.quantile(n, q) within each culture_condition,
over the obs passing the on-target and n_cells >= 50 filters. The DE h5ad
is not available locally, so the obs are taken from
data/DE_stats.suppl_table.csv (ontarget_significant, n_cells_target >= 50),
plus the deg0 perturbations with edges that this table misses (RHOH,
Stim8hr). That table is not identical to the h5ad obs: at q = 0.90 the
recomputed quantile falls on the other side of an integer DEG count in Rest
and Stim8hr. HUB_CUTOFF holds the integer cutoffs that reproduce every
published gene set exactly (01_collect.py re-checks this and fails on any
mismatch); at 0.99 / 0.95 they equal the recomputed quantile. In Stim8hr
at q = 0.90 cutoffs 152 and 153 both reproduce the sets; 152 is used
(UQCC1, n = 153, treated as a hub).

Everything degree-dependent (edges, universes, expected Jaccard) is
computed on each variant's hub-filtered graph: deg0 edges minus all edges
out of that variant's hubs.
"""

import importlib.util
import os
import sys
from functools import lru_cache
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent), str(HERE.parent / "heritability")]
# common.py reads the design from CH2_VARIANT at import: point it at deg0,
# whose DE_results.csv is the unfiltered edge list.
os.environ["CH2_VARIANT"] = "deg0_cells50_ot1_p0.05_z2.0_lfc0.5_h40.8_k1-1_w10000_snp1000"

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import polars as pl  # noqa: E402

from common import (CONDS, DATA, DESIGNS, FIGURES, ROOT, ldsc_dir,  # noqa: E402,F401
                    load_coloc, load_diseases, load_edges, load_gene_sets,
                    universe)

BASE = DESIGNS["deg0"]
# (key, % of regulators excluded, results/ldsc folder)
VARIANTS = [("deg0", 0, BASE), ("hub0.99", 1, "hub0.99_" + BASE),
            ("hub0.95", 5, "hub0.95_" + BASE), ("hub0.90", 10, "hub0.90_" + BASE)]
VKEYS = [v for v, _, _ in VARIANTS]
PCT = {v: p for v, p, _ in VARIANTS}
FOLDER = {v: f for v, _, f in VARIANTS}
QUANTILE = {"hub0.99": 0.99, "hub0.95": 0.95, "hub0.90": 0.90}
HUB_CUTOFF = {  # drop n_total_de_genes > cutoff (see module docstring)
    "hub0.99": {"Rest": 1355, "Stim8hr": 1996, "Stim48hr": 1566},
    "hub0.95": {"Rest": 374, "Stim8hr": 479, "Stim48hr": 426},
    "hub0.90": {"Rest": 130, "Stim8hr": 152, "Stim48hr": 131},
}
# Ordinal ramp, light -> dark = more hubs excluded (validated: adjacent CVD dE >= 13).
VCOLOR = dict(zip(VKEYS, ["#93c0ea", "#3f86cc", "#1d5ba0", "#0c2f5c"]))
VLABEL = {v: f"{PCT[v]} %" for v in VKEYS}

MIN_ABS_LOGFC = 0.5
H4 = 0.8
MAIN_COND = "Stim8hr"              # condition of main Fig. 2.6A
SET_ANNOTS = ["regulators_H4", "targets_H4"]
# Related-GWAS pairs (shared samples), by trait code in diseases.csv.
RELATED = [("ASTCO", "ASTAO"), ("RA_1", "RA_2"), ("IBD", "CD"), ("IBD", "UC")]
CONSERVATIVE = "conservative: ignores positive covariance between nested annotations"

OUT = FIGURES / "hub_exclusion"
MAIN, SUPP, DOUT = OUT / "main", OUT / "supp", OUT / "data"


def out_dirs():
    for d in (MAIN, SUPP, DOUT):
        d.mkdir(parents=True, exist_ok=True)


def collect_module():
    """vis/heritability/01_collect_results.py as a module (collect, add_fdr)."""
    p = HERE.parent / "heritability" / "01_collect_results.py"
    spec = importlib.util.spec_from_file_location("collect_results", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@lru_cache
def deg0_edges():
    """deg0 edges as the gene sets use them (|logFC| >= 0.5, no self-edges)."""
    return load_edges(MIN_ABS_LOGFC)


@lru_cache
def obs_table():
    """culture_condition, ensg, symbol, n_degs for the obs passing the
    pre-hub filters (see module docstring)."""
    s = pl.read_csv(DATA / "DE_stats.suppl_table.csv", infer_schema_length=100000)
    s = s.filter(pl.col("ontarget_significant"), pl.col("n_cells_target") >= 50).select(
        "culture_condition", pl.col("target_contrast").alias("ensg"),
        pl.col("target_contrast_gene_name").alias("symbol"),
        pl.col("n_total_de_genes").alias("n_degs"))
    de = pl.read_csv(ldsc_dir(BASE) / "DE_results.csv",
                     columns=["culture_condition", "cis_ensg", "cis_gene", "n_total_de_genes"])
    de = de.unique(["culture_condition", "cis_ensg"]).select(
        "culture_condition", pl.col("cis_ensg").alias("ensg"),
        pl.col("cis_gene").alias("symbol"), pl.col("n_total_de_genes").alias("n_degs"))
    return pl.concat([s, de]).unique(["culture_condition", "ensg"], keep="first")


@lru_cache
def hubs(variant, cond):
    """Ensembl ids dropped as hubs in this variant and condition."""
    if variant == "deg0":
        return frozenset()
    o = obs_table().filter(pl.col("culture_condition") == cond,
                           pl.col("n_degs") > HUB_CUTOFF[variant][cond])
    return frozenset(o["ensg"])


@lru_cache
def edges(variant):
    """Hub-filtered edge list for this variant."""
    e = deg0_edges()
    if variant == "deg0":
        return e
    parts = [e.filter(pl.col("culture_condition") == c,
                      ~pl.col("cis_ensg").is_in(list(hubs(variant, c)))) for c in CONDS]
    return pl.concat(parts)


@lru_cache
def anchors():
    """study -> set of anchor Ensembl ids (CD4 coloc PP.H4 > 0.8)."""
    c = load_coloc().filter(pl.col("H4") > H4)
    return {s: frozenset(g["Gene_id"]) for (s,), g in c.group_by("imd")}


@lru_cache
def gene_sets(variant):
    return load_gene_sets(ldsc_dir(FOLDER[variant]) / "gene_sets")


def perturbed_anchors(variant, study, cond):
    """Anchors with >= 1 out-edge in this variant's graph (they seed the target set)."""
    e = edges(variant).filter(pl.col("culture_condition") == cond)
    return set(anchors().get(study, ())) & set(e["cis_ensg"])


def lost_anchors(variant, study, cond):
    """Anchors seeding the deg0 target set that this variant removed as hubs."""
    return perturbed_anchors("deg0", study, cond) - perturbed_anchors(variant, study, cond)


def studies():
    """Traits with an anchor set, in diseases.csv order."""
    d = load_diseases()
    return [s for s in d.index if (s, "Rest", "H4_genes") in gene_sets("deg0")]


def trait_code(diseases, s):
    return diseases.loc[s, "trait"]


def pair_class(diseases, a, b):
    """IMD-IMD / IMD-control / control-control, or 'related' for same-GWAS pairs."""
    ta, tb = trait_code(diseases, a), trait_code(diseases, b)
    if (ta, tb) in RELATED or (tb, ta) in RELATED:
        return "related"
    g = sorted(diseases.loc[[a, b], "group"].astype(str))
    return {("IMD", "IMD"): "IMD–IMD", ("IMD", "non-IMD"): "IMD–control",
            ("non-IMD", "non-IMD"): "control–control"}[tuple(g)]


PAIR_CLASSES = ["IMD–IMD", "IMD–control", "control–control"]


def read(name):
    return pd.read_csv(DOUT / name)
