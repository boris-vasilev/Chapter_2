#!/usr/bin/env python
"""Scatter plot: gene-level vs SNP-level overlap of H4 gene sets.

For every (study, condition) the three H4 annotation sets

  G = H4_genes          R = regulators_H4         T = targets_H4

are compared pairwise (G-T, G-R, T-R) at two scales (gene IDs vs bp
footprint) with either of two coefficients:

  overlap coefficient (Szymkiewicz–Simpson):  |A ∩ B| / min(|A|, |B|)
  Jaccard index:                              |A ∩ B| / |A ∪ B|

Both are computed and written to --csv-out; --metric selects which one
populates the scatter axes.

The SNP-level metric is a bp-overlap proxy: each gene set is expanded to
its genomic footprint (CDS ± window-bp, merged across chromosomes) and
treated as a continuous interval set. Under approximately uniform SNP
density this is monotonic with actual SNP-count overlap, so it stands in
for "how similar do the LDSC annotation columns look" without needing a
1000G bim file.

Points where either axis > --label-threshold (default 0.5) are labelled
as "{trait}-{cond_code}:{pair}", e.g. "AITD-S8:T-R".
"""

import argparse
import gzip
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

COND_CODE = {"Rest": "R", "Stim8hr": "S8", "Stim48hr": "S48"}
DEFAULT_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y"}

ANNOTS = ["H4_genes", "regulators_H4", "targets_H4"]
CODES = {"H4_genes": "G", "regulators_H4": "R", "targets_H4": "T"}
PAIRS = [("H4_genes", "targets_H4"),
         ("H4_genes", "regulators_H4"),
         ("targets_H4", "regulators_H4")]


def load_gene_coords(gtf_path, chroms, protein_coding_only):
    """gene_id -> (chrom, start, end) from an Ensembl GTF."""
    coords = {}
    gid_re = re.compile(r'gene_id "([^"]+)"')
    bt_re = re.compile(r'gene_biotype "([^"]+)"')
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            chrom = cols[0]
            if chrom not in chroms:
                continue
            attr = cols[8]
            m = gid_re.search(attr)
            if not m:
                continue
            if protein_coding_only:
                mb = bt_re.search(attr)
                if mb is None or mb.group(1) != "protein_coding":
                    continue
            coords[m.group(1)] = (chrom, int(cols[3]), int(cols[4]))
    return coords


def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    out = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def intersect_sorted(a, b):
    """Intersection of two sorted lists of disjoint (start, end)."""
    i = j = 0
    out = []
    while i < len(a) and j < len(b):
        lo = max(a[i][0], b[j][0])
        hi = min(a[i][1], b[j][1])
        if lo < hi:
            out.append((lo, hi))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def total_bp(per_chrom):
    return sum(e - s for ivals in per_chrom.values() for s, e in ivals)


def build_intervals(gene_ids, coords, window):
    """gene set -> {chrom: sorted disjoint [(start, end), ...]}."""
    per_chrom = defaultdict(list)
    for g in gene_ids:
        hit = coords.get(g)
        if hit is None:
            continue
        chrom, start, end = hit
        per_chrom[chrom].append((max(0, start - window), end + window))
    return {c: merge_intervals(ivals) for c, ivals in per_chrom.items()}


def bp_overlap(a, b):
    """Returns (intersection_bp, |A|_bp, |B|_bp)."""
    inter_bp = 0
    for chrom in set(a) & set(b):
        inter_bp += sum(e - s for s, e in intersect_sorted(a[chrom], b[chrom]))
    return inter_bp, total_bp(a), total_bp(b)


def coeffs(inter, size_a, size_b):
    """Returns (overlap_coefficient, jaccard)."""
    if size_a == 0 or size_b == 0:
        return 0.0, 0.0
    overlap = inter / min(size_a, size_b)
    union = size_a + size_b - inter
    jaccard = inter / union if union else 0.0
    return overlap, jaccard


def discover_panels(gene_sets_dir):
    """Yield ((study, cond), {annot: path}) for every full H4 triple."""
    annot_re = "|".join(re.escape(a) for a in sorted(ANNOTS, key=len, reverse=True))
    pat = re.compile(rf"^(?P<study>.+)_(?P<cond>[^_]+)_(?P<annot>{annot_re})\.txt$")
    bag = defaultdict(dict)
    for p in sorted(Path(gene_sets_dir).glob("*.txt")):
        m = pat.match(p.name)
        if not m:
            continue
        bag[(m.group("study"), m.group("cond"))][m.group("annot")] = p
    for key, paths in bag.items():
        if set(paths) == set(ANNOTS):
            yield key, paths


def read_genes(path):
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def plot(rows, out_path, label_threshold, window_kb, metric):
    fig, ax = plt.subplots(figsize=(20, 17))
    x_key, y_key = f"snp_{metric}", f"gene_{metric}"
    xs = [r[x_key] for r in rows]
    ys = [r[y_key] for r in rows]
    hot = [(x > label_threshold or y > label_threshold) for x, y in zip(xs, ys)]

    cold_x = [x for x, h in zip(xs, hot) if not h]
    cold_y = [y for y, h in zip(ys, hot) if not h]
    hot_x = [x for x, h in zip(xs, hot) if h]
    hot_y = [y for y, h in zip(ys, hot) if h]

    ax.scatter(cold_x, cold_y, color="#1f77b4", alpha=0.55, s=260, edgecolor="none")
    ax.scatter(hot_x, hot_y, color="#ff7f0e", alpha=0.95, s=360, edgecolor="none")

    lim = max(max(xs + [1.0]), max(ys + [1.0]))
    ax.plot([0, lim], [0, lim], linestyle="--", color="#4a90d9", alpha=0.6,
            linewidth=4.0)

    for r, is_hot in zip(rows, hot):
        if not is_hot:
            continue
        ax.annotate(
            r["label"], (r[x_key], r[y_key]),
            xytext=(12, 0), textcoords="offset points",
            fontsize=34, va="center",
        )

    metric_label = "overlap coefficient" if metric == "overlap" else "Jaccard index"
    ax.set_xlabel(f"SNP-level {metric_label}", fontsize=44)
    ax.set_ylabel(f"Gene-level {metric_label}", fontsize=44)
    ax.set_title("Gene vs SNP Overlap", fontsize=48, fontweight="bold")
    ax.tick_params(axis="both", labelsize=36)
    ax.text(0.02, 0.97, f"±{window_kb:g}KB", transform=ax.transAxes,
            color="red", fontsize=48, fontweight="bold", va="top")
    ax.set_xlim(-0.02, lim + 0.02)
    ax.set_ylim(-0.02, lim + 0.02)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-sets-dir", required=True,
                        help="e.g. data/gene_sets/coloc")
    parser.add_argument("--gtf", required=True,
                        help="Ensembl GTF (gzipped OK).")
    parser.add_argument("--window-kb", type=float, default=10.0,
                        help="bp window each side of every gene's CDS for "
                             "the SNP-level metric (default 10kb).")
    parser.add_argument("--protein-coding-only", action="store_true")
    parser.add_argument("--label-threshold", type=float, default=0.5)
    parser.add_argument("--metric", choices=["overlap", "jaccard"],
                        default="overlap",
                        help="Coefficient shown on the plot axes (both are "
                             "always written to --csv-out).")
    parser.add_argument("--out", required=True)
    parser.add_argument("--csv-out", default=None,
                        help="Optional path for the per-pair table.")
    args = parser.parse_args()

    print(f"Loading GTF {args.gtf}…")
    coords = load_gene_coords(args.gtf, DEFAULT_CHROMS, args.protein_coding_only)
    print(f"  {len(coords)} gene coordinates loaded")

    window_bp = int(args.window_kb * 1000)
    rows = []

    for (study, cond), paths in discover_panels(args.gene_sets_dir):
        if cond not in COND_CODE:
            print(f"  skip {study}/{cond}: unknown condition")
            continue
        genes = {a: read_genes(paths[a]) for a in ANNOTS}
        intervals = {a: build_intervals(genes[a], coords, window_bp) for a in ANNOTS}
        trait_short = study.split("_")[0]
        cond_code = COND_CODE[cond]

        for a, b in PAIRS:
            pair = f"{CODES[a]}-{CODES[b]}"
            gene_inter = len(genes[a] & genes[b])
            gene_overlap, gene_jaccard = coeffs(
                gene_inter, len(genes[a]), len(genes[b])
            )
            snp_inter, snp_a, snp_b = bp_overlap(intervals[a], intervals[b])
            snp_overlap, snp_jaccard = coeffs(snp_inter, snp_a, snp_b)
            rows.append({
                "study": study, "trait_short": trait_short,
                "condition": cond, "pair": pair,
                "n_genes_a": len(genes[a]), "n_genes_b": len(genes[b]),
                "gene_overlap": gene_overlap, "gene_jaccard": gene_jaccard,
                "snp_overlap": snp_overlap, "snp_jaccard": snp_jaccard,
                "label": f"{trait_short}-{cond_code}:{pair}",
            })

    if not rows:
        raise SystemExit("No (study, condition) triples found.")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plot(rows, out_path, args.label_threshold, args.window_kb, args.metric)

    if args.csv_out:
        import csv
        keys = ["study", "trait_short", "condition", "pair",
                "n_genes_a", "n_genes_b",
                "gene_overlap", "gene_jaccard",
                "snp_overlap", "snp_jaccard", "label"]
        with open(args.csv_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in rows:
                w.writerow(r)
        print(f"Wrote {args.csv_out}")


if __name__ == "__main__":
    main()
