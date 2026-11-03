#!/usr/bin/env python
"""Programme over-representation of coloc / regulator / target gene sets.

One figure per culture condition. Within each figure, three stacked panels:

    Panel 1  Coloc (H4) genes
    Panel 2  Regulators
    Panel 3  Targets

Rows = diseases (studies), columns = Perturb-seq programmes. Cell colour =
log2(odds ratio) of a one-sided Fisher exact test asking whether that
disease's gene set is over-represented in that programme. ``*``/``**``/``***``
mark nominal p < 0.05 / 0.01 / 0.001 and a gold border marks BH-FDR < 0.05 —
the same aesthetic as approved_drug_enrichment_all.py. ``·`` = not tested
(no overlap).

------------------------------------------------------------------------
STATISTICS  (everything is done within one condition c)

  universe N : perturbed genes (``perturb_ensg``) assigned to *any*
               programme in condition c. Programme membership is only
               defined for perturbed genes, so this is the only valid
               background.
  programme P: a cluster in condition c;  K = |universe ∩ P|.
  gene set S : the disease's genes for (annot, condition c);
               foreground = S ∩ universe,  n = |foreground|.
  overlap k  : |foreground ∩ P|.
  test       : Fisher exact on [[k, n-k], [K-k, N-K-(n-k)]] with
               alternative="greater". This equals the hypergeometric
               P(X >= k) your old scripts used, but also returns the odds
               ratio used for the colour.

Only (cell, programme) pairs with k >= 1 are tested; everything else is
"no overlap" and drawn as ``·``. BH-FDR is computed across *every* test in
the whole run, so it is invariant to how the output is split into files.

Caveat worth knowing: regulators are themselves perturbed genes, so they
lie fully inside the universe. Coloc (H4) and target genes are trans/effect
genes — only the subset that is *also* perturbed-and-assigned is testable,
so those panels test fewer genes per disease.
------------------------------------------------------------------------
"""

import argparse
import re
from collections import defaultdict
from math import log10
from pathlib import Path

import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import pdist
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]


# ---- programme universe ------------------------------------------------

def build_program_tables(programs_csv):
    """Per-condition lookups derived from program_assignments.csv.

    Returns:
      gene_to_cluster[cond][ensg] -> cluster id
      cluster_size[cond][cluster] -> # assigned perturbed genes
      cluster_meta[cluster]       -> (program_number, manual_annotation)
    """
    df = pl.read_csv(programs_csv)
    needed = {"perturb_ensg", "condition", "cluster", "manual_annotation"}
    missing = needed - set(df.columns)
    if missing:
        raise SystemExit(f"program_assignments missing columns: {missing}")
    has_program = "program" in df.columns

    gene_to_cluster = defaultdict(dict)
    cluster_size = defaultdict(lambda: defaultdict(int))
    cluster_meta = {}
    for r in df.iter_rows(named=True):
        c = r["condition"]
        if c not in COND_ORDER:
            continue
        cl = int(r["cluster"])
        ensg = r["perturb_ensg"]
        if ensg in gene_to_cluster[c]:
            continue  # first assignment per gene per condition wins
        gene_to_cluster[c][ensg] = cl
        cluster_size[c][cl] += 1
        prog = int(r["program"]) if has_program and r["program"] is not None else cl
        cluster_meta.setdefault(cl, (prog, r["manual_annotation"] or f"cluster {cl}"))
    return gene_to_cluster, cluster_size, cluster_meta


# ---- gene sets ---------------------------------------------------------

def read_gene_set(path):
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def discover_studies(gene_sets_dir, annots, conds):
    """Study names that have at least one (cond, annot) file present."""
    studies = set()
    for annot in annots:
        for p in Path(gene_sets_dir).glob(f"*_{annot}.txt"):
            stem = p.name[: -len(f"_{annot}.txt")]  # "{study}_{cond}"
            for c in conds:
                if stem.endswith(f"_{c}"):
                    studies.add(stem[: -len(f"_{c}")])
                    break
    return sorted(studies)


# ---- testing -----------------------------------------------------------

def run_tests(gene_sets_dir, studies, annot_map, conds,
              gene_to_cluster, cluster_size, cluster_meta,
              min_overlap, min_program_size):
    """One Fisher test per (study, condition, annot, programme) with k >= 1."""
    rows = []
    gsdir = Path(gene_sets_dir)
    for cond in conds:
        g2c = gene_to_cluster[cond]
        universe = set(g2c.keys())
        N = len(universe)
        if N == 0:
            continue
        for annot_key, annot in annot_map.items():
            for study in studies:
                path = gsdir / f"{study}_{cond}_{annot}.txt"
                if not path.exists():
                    continue
                fg = read_gene_set(path) & universe
                n = len(fg)
                if n < min_overlap:
                    continue
                k_by_cluster = defaultdict(int)
                for g in fg:
                    k_by_cluster[g2c[g]] += 1
                for cl, k in k_by_cluster.items():
                    K = cluster_size[cond][cl]
                    if K < min_program_size:
                        continue
                    # 2x2: in-set/out-of-set  x  in-programme/out-of-programme
                    a = k
                    b = n - k
                    c_ = K - k
                    d = N - K - (n - k)
                    odds_ratio, p = fisher_exact([[a, b], [c_, d]],
                                                 alternative="greater")
                    expected = n * K / N
                    log2_fe = (log10(k / expected) / log10(2)
                               if expected > 0 and k > 0 else float("nan"))
                    prog, annotation = cluster_meta.get(cl, (cl, f"cluster {cl}"))
                    rows.append({
                        "condition": cond,
                        "annot_key": annot_key,
                        "annot": annot,
                        "study": study,
                        "cluster": cl,
                        "program": prog,
                        "programme": annotation,
                        "k": k, "n": n, "K": K, "N": N,
                        "expected": expected,
                        "odds_ratio": float(odds_ratio),
                        "log2_or": float(np.log2(odds_ratio))
                        if odds_ratio not in (0.0, np.inf) else
                        (np.inf if odds_ratio == np.inf else -np.inf),
                        "log2_fold_enrichment": log2_fe,
                        "p_value": float(p),
                    })
    return rows


# ---- labelling / ordering ----------------------------------------------

# Readable disease names, keyed by study code — kept in sync with
# approved_drug_enrichment_all.py / disease_specificity_matrix.py so the
# heritability and drug-target figures all use the same labels.
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
    "MS": "Multiple sclerosis",
}


def make_short_labels(studies):
    """Readable disease labels from the study code (e.g. ``RA_Ishigaki_..`` ->
    ``Rheum. arthritis``), matching the heritability and drug-target figures.
    Studies that share a code (the two RA GWAS) are disambiguated by appending
    the trailing token, e.g. ``Rheum. arthritis (1)`` / ``(2)``."""
    code = {s: s.split("_")[0] for s in studies}
    counts = defaultdict(int)
    for c in code.values():
        counts[c] += 1
    label = {}
    for s in studies:
        c = code[s]
        name = SHORT_NAMES.get(c, c)
        label[s] = f"{name} ({s.split('_')[-1]})" if counts[c] > 1 else name
    return label


def order_studies(studies, groups_path):
    """Return (ordered_studies, label_of, imd_set). IMD first if groups given."""
    if not groups_path:
        return sorted(studies), make_short_labels(studies), set()
    g = pl.read_csv(groups_path)
    label = {r["study"]: r["short_label"] for r in g.iter_rows(named=True)}
    group = {r["study"]: r["group"] for r in g.iter_rows(named=True)}
    imd = {s for s in studies if group.get(s) == "IMD"}
    ordered = sorted(studies, key=lambda s: (group.get(s, "zzz") != "IMD", s))
    label_of = {s: label.get(s, s.split("_")[0]) for s in studies}
    return ordered, label_of, imd


def select_programs(rows, cond, max_programs, sig_p):
    """Clusters with >=1 cell at nominal p < sig_p in this condition,
    ranked by max -log10(p); fall back to best-p clusters if none."""
    best = defaultdict(lambda: 1.0)
    for r in rows:
        if r["condition"] == cond:
            best[r["cluster"]] = min(best[r["cluster"]], r["p_value"])
    if not best:
        return []
    sig = {cl: p for cl, p in best.items() if p < sig_p}
    pool = sig or best  # fall back to all tested clusters if none significant
    ordered = sorted(pool, key=lambda cl: pool[cl])
    return ordered[:max_programs]


# ---- plotting ----------------------------------------------------------

COLOR_SPECS = {
    "log2or": ("log2_or", "RdBu_r", "log2(odds ratio)"),
    "log2fe": ("log2_fold_enrichment", "RdBu_r", "log2(fold enrichment)"),
    "neglog10p": (None, "YlOrRd", "-log10(p)"),
}


def _value(r, color):
    if color == "neglog10p":
        return -log10(max(r["p_value"], 1e-300))
    field = COLOR_SPECS[color][0]
    return r[field]


def draw_panel(ax, mat, pmat, fmat, prog_labels, row_labels, title,
               cmap, vmin, vmax, fdr_thr, show_xticks, imd_boundary):
    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad("white")
    im = ax.imshow(mat, cmap=cmap_obj, vmin=vmin, vmax=vmax, aspect="auto")

    n_rows, n_cols = mat.shape
    for i in range(n_rows):
        for j in range(n_cols):
            p = pmat[i, j]
            if np.isnan(p):
                ax.text(j, i, "·", ha="center", va="center",
                        color="grey", fontsize=12)
                continue
            stars = ("***" if p < 0.001 else "**" if p < 0.01
                     else "*" if p < 0.05 else "")
            if stars:
                t = ax.text(j, i, stars, ha="center", va="center",
                            color="white", fontsize=13, fontweight="bold")
                t.set_path_effects([path_effects.Stroke(linewidth=2,
                                    foreground="black"), path_effects.Normal()])
            if not np.isnan(fmat[i, j]) and fmat[i, j] < fdr_thr:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1.0, 1.0,
                             fill=False, edgecolor="#FFD700", linewidth=3.0))

    if imd_boundary and 0 < imd_boundary < n_rows:
        ax.axhline(imd_boundary - 0.5, color="black", linewidth=1.5)

    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(row_labels, fontsize=11)
    ax.set_xticks(range(n_cols))
    if show_xticks:
        ax.set_xticklabels(prog_labels, rotation=45, ha="right", fontsize=10)
    else:
        ax.set_xticklabels([])
    ax.tick_params(axis="both", length=0)
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)
    for s in ax.spines.values():
        s.set_linewidth(0.5)
    ax.set_title(title, fontsize=13, loc="left", pad=4, fontweight="bold")
    return im


def cluster_rows_by_regulators(rows, cond, studies, vmax):
    """Cluster diseases on their *regulator* enrichment profile — the dense,
    well-populated panel (the only one where every disease has signal).

    Feature vector per disease = log2(odds ratio) over *every* programme
    tested in this condition (no-overlap -> 0, i.e. no enrichment), so the
    full programme space is used, not just the displayed columns. Distance =
    correlation (groups diseases enriched in the *same* programmes), average
    linkage. Returns (perm, Z) for reordering + the inline dendrogram, or
    (None, None) when there is too little to cluster."""
    reg = [r for r in rows if r["condition"] == cond
           and r["annot_key"] == "regulators" and r["study"] in studies]
    if not reg or len(studies) < 3:
        return None, None
    clusters = sorted({r["cluster"] for r in reg})
    cidx = {c: j for j, c in enumerate(clusters)}
    sidx = {s: i for i, s in enumerate(studies)}
    X = np.zeros((len(studies), len(clusters)))
    for r in reg:
        v = r["log2_or"]
        if not np.isfinite(v):
            v = np.sign(v) * vmax
        X[sidx[r["study"]], cidx[r["cluster"]]] = float(np.clip(v, -vmax, vmax))
    if np.allclose(X, 0.0):
        return None, None
    D = np.nan_to_num(pdist(X, metric="correlation"), nan=1.0)  # flat rows -> max dist
    Z = linkage(D, method="average")
    return list(leaves_list(Z)), Z


def plot_condition(rows, cond, annot_keys, annot_titles, ordered_studies,
                   label_of, imd_set, programs, cluster_meta,
                   color, vmax, fdr_thr, out_path,
                   cluster_rows=False):
    field, cmap, cbar_label = COLOR_SPECS[color]
    studies = [s for s in ordered_studies
               if any(r["study"] == s and r["condition"] == cond for r in rows)]
    if not studies:
        print(f"  [{cond}] no studies with tests; skipping figure")
        return
    n_imd = sum(1 for s in studies if s in imd_set)

    prog_idx = {cl: j for j, cl in enumerate(programs)}
    study_idx = {s: i for i, s in enumerate(studies)}

    # disambiguate duplicate annotation labels by appending the programme #
    seen = defaultdict(int)
    for cl in programs:
        seen[cluster_meta[cl][1]] += 1
    prog_labels = []
    for cl in programs:
        prog_no, ann = cluster_meta[cl]
        ann = ann if len(ann) <= 28 else ann[:27] + "…"
        prog_labels.append(f"{ann} (P{prog_no})" if seen[cluster_meta[cl][1]] > 1
                           else ann)

    shape = (len(studies), len(programs))
    mats = {a: np.full(shape, np.nan) for a in annot_keys}
    pmats = {a: np.full(shape, np.nan) for a in annot_keys}
    fmats = {a: np.full(shape, np.nan) for a in annot_keys}
    for r in rows:
        if (r["condition"] != cond or r["annot_key"] not in mats
                or r["cluster"] not in prog_idx or r["study"] not in study_idx):
            continue
        i, j = study_idx[r["study"]], prog_idx[r["cluster"]]
        val = _value(r, color)
        if color != "neglog10p":
            val = float(np.clip(val, -vmax, vmax))
        mats[r["annot_key"]][i, j] = val
        pmats[r["annot_key"]][i, j] = r["p_value"]
        fmats[r["annot_key"]][i, j] = r["fdr"]

    clustered = False
    dendro_Z = None
    if cluster_rows:
        perm, dendro_Z = cluster_rows_by_regulators(rows, cond, studies, vmax)
        if perm is not None:
            studies = [studies[i] for i in perm]
            for a in annot_keys:
                mats[a] = mats[a][perm]
                pmats[a] = pmats[a][perm]
                fmats[a] = fmats[a][perm]
            clustered = True

    if color == "neglog10p":
        finite = np.concatenate([m[np.isfinite(m)] for m in mats.values()])
        vlo, vhi = 0.0, (float(np.max(finite)) if finite.size else 1.0)
    else:
        vlo, vhi = -vmax, vmax

    n_panels = len(annot_keys)
    fig_w = max(9, 0.45 * len(programs) + 4)
    fig_h = max(6, 0.32 * len(studies) * n_panels + 2)
    fig, axes = plt.subplots(n_panels, 1, figsize=(fig_w, fig_h),
                             squeeze=False)
    axes = axes[:, 0]
    cax = fig.add_axes([0.92, 0.25, 0.012, 0.5])

    im = None
    for idx, a in enumerate(annot_keys):
        im = draw_panel(
            axes[idx], mats[a], pmats[a], fmats[a], prog_labels,
            [label_of[s] for s in studies], annot_titles[a],
            cmap, vlo, vhi, fdr_thr,
            show_xticks=(idx == n_panels - 1),
            imd_boundary=0 if clustered else (n_imd if imd_set else 0),
        )

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(cbar_label, fontsize=12)
    cbar.ax.tick_params(labelsize=10)

    cluster_note = ("\nrows clustered by regulator enrichment profile "
                    "(dendrogram, left)") if clustered else ""
    fig.suptitle(
        f"Programme over-representation — {cond}\n"
        f"Fisher exact: * p<0.05  ** p<0.01  *** p<0.001;  "
        f"gold border: BH-FDR < {fdr_thr};  · = no overlap" + cluster_note,
        fontsize=14, fontweight="bold",
    )
    left = 0.30 if clustered else 0.16
    fig.subplots_adjust(left=left, right=0.90, top=0.90, bottom=0.16, hspace=0.22)

    # Inline dendrogram strip, aligned to the (regulators) panel it is built
    # from. All panels share this row order.
    if clustered and dendro_Z is not None:
        reg_idx = (annot_keys.index("regulators")
                   if "regulators" in annot_keys else 0)
        pos = axes[reg_idx].get_position()
        dax = fig.add_axes([0.045, pos.y0, 0.085, pos.height])
        dendrogram(dendro_Z, orientation="left", ax=dax, no_labels=True,
                   color_threshold=0, above_threshold_color="#555555")
        # leaves_list order runs bottom->top in a 'left' dendrogram; the
        # heatmap has row 0 at the top, so flip to match.
        dax.invert_yaxis()
        dax.axis("off")
        dax.set_title("regulator\nclustering", fontsize=8, pad=2)

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


# ---- CLI ---------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--programs", required=True,
                   help="data/programs_zhu/program_assignments.csv")
    p.add_argument("--gene-sets-dir", required=True)
    p.add_argument("--studies", nargs="+", default=None,
                   help="Subset of studies; default = all in gene-sets-dir.")
    p.add_argument("--conditions", nargs="+", default=COND_ORDER)
    p.add_argument("--coloc-annot", default="H4_genes")
    p.add_argument("--regulators-annot", default="regulators_H4")
    p.add_argument("--targets-annot", default="targets_H4")
    p.add_argument("--coloc-label", default="Coloc (H4) genes")
    p.add_argument("--regulators-label", default="Regulators")
    p.add_argument("--targets-label", default="Targets")
    p.add_argument("--disease-groups", default=None,
                   help="Optional disease_groups.csv for short labels + "
                        "IMD-first row ordering (from make_gene_sets_tenk10k.py).")
    p.add_argument("--color", choices=list(COLOR_SPECS), default="log2or",
                   help="Cell colour metric (default log2or, like the drug plot).")
    p.add_argument("--vmax", type=float, default=4.0,
                   help="Symmetric colour cap for log2 metrics (ignored for "
                        "neglog10p, which auto-scales).")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--min-overlap", type=int, default=2,
                   help="Skip a (disease, annot) cell whose gene set has fewer "
                        "than this many programme-assigned genes.")
    p.add_argument("--min-program-size", type=int, default=3)
    p.add_argument("--max-programs", type=int, default=40,
                   help="Max programme columns per figure (top by significance).")
    p.add_argument("--sig-p", type=float, default=0.05,
                   help="Nominal p a programme must reach in >=1 cell to be "
                        "shown as a column.")
    p.add_argument("--cluster-rows", action=argparse.BooleanOptionalAction,
                   default=True,
                   help="Cluster the disease rows by their regulator enrichment "
                        "profile (correlation distance over all programmes) and "
                        "draw an inline dendrogram (default on). Replaces the "
                        "fixed/IMD-first ordering and drops the IMD divider; the "
                        "coloc and target panels inherit the same row order. Use "
                        "--no-cluster-rows to keep the grouped order.")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--format", default="png", choices=["png", "pdf"])
    p.add_argument("--csv-out", default=None,
                   help="Optional path for the full per-test table.")
    args = p.parse_args()

    g2c, csize, cmeta = build_program_tables(args.programs)

    annot_map = {
        "coloc": args.coloc_annot,
        "regulators": args.regulators_annot,
        "targets": args.targets_annot,
    }
    annot_titles = {
        "coloc": args.coloc_label,
        "regulators": args.regulators_label,
        "targets": args.targets_label,
    }
    annot_keys = list(annot_map)

    studies = args.studies or discover_studies(
        args.gene_sets_dir, list(annot_map.values()), args.conditions)
    if not studies:
        raise SystemExit("No studies found in gene-sets-dir.")
    print(f"studies: {len(studies)}  conditions: {args.conditions}")

    rows = run_tests(
        args.gene_sets_dir, studies, annot_map, args.conditions,
        g2c, csize, cmeta, args.min_overlap, args.min_program_size,
    )
    if not rows:
        raise SystemExit("No tests run — check gene-set / programme overlap.")

    # One global BH-FDR across every test in the run.
    fdr = multipletests([r["p_value"] for r in rows], method="fdr_bh")[1]
    for r, q in zip(rows, fdr):
        r["fdr"] = float(q)
    print(f"{len(rows)} tests; {sum(q < args.fdr for q in fdr)} pass BH-FDR<{args.fdr}")

    ordered_studies, label_of, imd_set = order_studies(studies, args.disease_groups)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for cond in args.conditions:
        programs = select_programs(rows, cond, args.max_programs, args.sig_p)
        if not programs:
            print(f"  [{cond}] no programmes to show; skipping")
            continue
        out_path = out_dir / f"program_enrichment_{cond}.{args.format}"
        plot_condition(
            rows, cond, annot_keys, annot_titles, ordered_studies, label_of,
            imd_set, programs, cmeta, args.color, args.vmax, args.fdr, out_path,
            cluster_rows=args.cluster_rows,
        )

    csv_out = args.csv_out or (out_dir / "program_enrichment.csv")
    pl.DataFrame(rows).sort(
        ["condition", "annot_key", "study", "p_value"]
    ).write_csv(csv_out)
    print(f"wrote {csv_out}")


if __name__ == "__main__":
    main()
