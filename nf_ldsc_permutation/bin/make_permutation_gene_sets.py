#!/usr/bin/env python
"""Gene sets for one S-LDSC permutation (v2 matched null).

For one study x condition x perm, write regulators_H4 / targets_H4 / H4_genes
gene-set files (one Ensembl ID per line, no header) built exactly as
make_gene_sets_coloc.py builds them, but with the anchor set replaced by a
matched pseudo-anchor set.

  --perm real   anchors = the real coloc H4 genes. If --reference-dir is given
                (e.g. data/sensitivity/thr0.5/gene_sets, the coloc_thr0.5 sets),
                the three sets are compared with it; any difference is a hard error.
  --perm <int>  anchors = pseudo set drawn by matching.Matcher, matched to the
                real anchors on in-degree, out-degree and expression stratum in
                THIS condition's edge set; seed derived from (base_seed + perm,
                condition) so each draw is reproducible in isolation.

A one-line matching log (exact draws, fallbacks, unmeasured draws) is written
next to the gene sets as {study}_{cond}_matching.tsv.

Runs under the Python 3.9 LDSC container: polars + numpy only.
"""

# Keep annotations lazy so this runs under the Python 3.9 LDSC container.
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from matching import (  # noqa: E402
    build_matcher, convergent_regulators, convergent_targets,
    load_coloc_genes, perm_rng,
)

ANNOTS = ("regulators_H4", "targets_H4", "H4_genes")


def write_gene_set(out_dir, study, cond, annot, genes):
    fname = out_dir / f"{study}_{cond}_{annot}.txt"
    # Sorted for byte-stable, run-order-independent output.
    genes = sorted({g for g in genes if g})
    if genes:
        pl.DataFrame({"GENE": genes}).write_csv(
            fname, separator="\t", include_header=False)
    else:
        fname.write_text("")  # empty set -> backend's low-SNP filter drops it
    print(f"  {fname.name}\t{len(genes)} genes", file=sys.stderr)
    return genes


def read_set(path):
    with open(path) as f:
        return sorted({line.strip() for line in f if line.strip()})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--de-results", required=True)
    p.add_argument("--coloc", required=True)
    p.add_argument("--expression", required=True,
                   help="screen_expression.tsv from build_screen_expression.py")
    p.add_argument("--gene-coord", required=True,
                   help="The study's gene-coordinate file (unmeasured pool).")
    p.add_argument("--protein-coding-coord", required=True)
    p.add_argument("--study", required=True)
    p.add_argument("--condition", required=True)
    p.add_argument("--cell-type", default="CD4")
    p.add_argument("--h4-threshold", type=float, default=0.8)
    p.add_argument("--min-anchors-per-regulator", type=int, default=1)
    p.add_argument("--min-perturbations-per-target", type=int, default=1)
    p.add_argument("--min-abs-logfc", type=float, default=0.5)
    p.add_argument("--n-degree-bins", type=int, default=4)
    p.add_argument("--n-expr-bins", type=int, default=5)
    p.add_argument("--knn-k", type=int, default=20,
                   help="k for out-degree nearest-neighbour matching of "
                        "knocked-down anchors.")
    p.add_argument("--perm", required=True,
                   help="'real' for the observed test, or an integer perm index")
    p.add_argument("--base-seed", type=int, default=42)
    p.add_argument("--reference-dir", default=None,
                   help="For --perm real: sets to compare with (optional).")
    p.add_argument("--out-dir", required=True)
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    study, cond = args.study, args.condition

    de = pl.read_csv(args.de_results,
                     columns=["cis_ensg", "trans_ensg", "culture_condition", "logFC"])
    real = load_coloc_genes(args.coloc, args.cell_type, args.h4_threshold, study)
    if not real:
        raise SystemExit(f"No coloc H4 anchors for {study} "
                         f"({args.cell_type}, H4>{args.h4_threshold})")
    all_coloc = load_coloc_genes(args.coloc, args.cell_type, args.h4_threshold)
    matcher, edges = build_matcher(
        de, cond, args.min_abs_logfc, args.expression, all_coloc,
        args.gene_coord, args.protein_coding_coord,
        args.n_degree_bins, args.n_expr_bins, args.knn_k)

    if args.perm == "real":
        anchors = real
        log = {"n_exact": None, "n_fallback": None, "fallback_dist": [],
               "n_knn": None, "n_knn_relaxed": None, "knn_relax_dist": [],
               "n_unmeasured": None}
        print(f"[real] {len(anchors)} coloc anchors", file=sys.stderr)
    else:
        perm_idx = int(args.perm)
        anchors, log = matcher.sample(real, perm_rng(args.base_seed, perm_idx, cond))
        print(f"[perm {perm_idx}] {len(anchors)} pseudo-anchors (target "
              f"{len(real)}): exact={log['n_exact']} fallback={log['n_fallback']} "
              f"knn={log['n_knn']} (relaxed {log['n_knn_relaxed']}) "
              f"unmeasured={log['n_unmeasured']}", file=sys.stderr)

    sets = {
        "regulators_H4": convergent_regulators(edges, anchors,
                                               args.min_anchors_per_regulator),
        "targets_H4": convergent_targets(edges, anchors,
                                         args.min_perturbations_per_target),
        "H4_genes": anchors,
    }
    written = {a: write_gene_set(out_dir, study, cond, a, sets[a]) for a in ANNOTS}

    pl.DataFrame([{
        "study": study, "condition": cond, "perm": args.perm,
        "n_real_anchors": len(real), "n_anchors": len(written["H4_genes"]),
        "n_exact": log["n_exact"], "n_fallback": log["n_fallback"],
        "max_fallback_dist": max(log["fallback_dist"]) if log["fallback_dist"] else None,
        "n_knn": log["n_knn"], "n_knn_relaxed": log["n_knn_relaxed"],
        "max_knn_relax_dist": max(log["knn_relax_dist"]) if log["knn_relax_dist"] else None,
        "n_unmeasured": log["n_unmeasured"],
        "n_measured_pool": len(matcher.measured_pool),
        "n_regulators": len(written["regulators_H4"]),
        "n_targets": len(written["targets_H4"]),
    }]).write_csv(out_dir / f"{study}_{cond}_matching.tsv", separator="\t")

    if args.perm == "real" and args.reference_dir:
        bad = []
        for a in ANNOTS:
            ref = read_set(Path(args.reference_dir) / f"{study}_{cond}_{a}.txt")
            if ref != written[a]:
                bad.append(f"{a}: new {len(written[a])} vs reference {len(ref)}")
        if bad:
            for f in ANNOTS:
                (out_dir / f"{study}_{cond}_{f}.txt").unlink()
            raise SystemExit(f"[real] {study} {cond} differs from "
                             f"{args.reference_dir}: " + "; ".join(bad))
        print(f"[real] identical to {args.reference_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
