#!/usr/bin/env python
"""S-LDSC heritability enrichment of Perturb-seq programmes, per disease.

Companion to ``program_enrichment.py`` (which tests gene-set over-representation
in programmes). Here the signal is the **S-LDSC partitioned-heritability
enrichment** of each Zhu programme annotation for each disease GWAS, parsed into
``all_results_parsed.csv`` (one row per programme × disease × condition).

One figure per culture condition. Rows = diseases (studies), columns =
programmes. Cell colour = heritability ``enrichment`` (prop_h2 / prop_snps),
diverging around 1 (= no enrichment). ``*``/``**``/``***`` mark the one-sided
enrichment test p < 0.05 / 0.01 / 0.001 from z = (enrichment - 1) / SE, and a
gold border marks BH-FDR < 0.05 (computed across every test in the run).

Disease rows are hierarchically clustered by their enrichment profile across
*all* programmes (correlation distance, average linkage) and an inline
dendrogram is drawn on the left — the same treatment as program_enrichment.py.
Row labels use the shared SHORT_NAMES, with the two RA GWAS disambiguated.

The matrix is dense (every study is tested against every programme), so unlike
the gene-set version there are no "no-overlap" cells.

Example:
    python perturb_ldsc_pipeline/vis/ldsc_program_enrichment.py \
        --results plots/ldsc_zhu_diagnostic/all_results_parsed.csv \
        --out-dir plots/ldsc_program_enrichment
"""

import argparse
from collections import defaultdict
from math import log10
from pathlib import Path

import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import pdist
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

COND_ORDER = ["Rest", "Stim8hr", "Stim48hr"]

# Readable disease names, keyed by study code — kept in sync with
# approved_drug_enrichment_all.py / program_enrichment.py.
SHORT_NAMES = {
    "AST": "Asthma", "ASTAO": "Asthma (AO)", "ASTCO": "Asthma (CO)",
    "ATD": "Atopic derm.", "AITD": "Autoimm. thyroid", "CD": "Crohn's",
    "IBD": "IBD", "UC": "Ulcer. colitis", "RA": "Rheum. arthritis",
    "SLE": "SLE", "T1D": "T1D", "IS": "Ischemic stroke",
    "MYOI": "Myocard. infarct.", "OA": "Osteoarthritis", "T2D": "T2D",
    "ADFH": "Alzheimer's", "SCZ": "Schizophrenia", "MS": "Multiple sclerosis",
}

# colour metric -> (csv column, colourbar label)
COLOR_SPECS = {
    "enrichment": ("enrichment", "heritability enrichment"),
    "tau_z": ("tau_z", "tau z-score"),
}


def make_short_labels(studies):
    """Readable disease labels from the study code; studies that share a code
    (the two RA GWAS) get a trailing ``(1)`` / ``(2)``."""
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


def load_results(path, sig_source):
    """Read the parsed S-LDSC table and attach a one-sided (enriched) p-value
    and a global BH-FDR, computed from ``sig_source``:

      * ``tau_z``      — z-score of the S-LDSC coefficient (the standard,
                         conservative test that the programme contributes to
                         heritability *conditional on the baseline model*).
                         Recommended: it isolates disease-specific signal.
      * ``enrichment`` — z = (enrichment - 1) / SE. Flags almost every cell,
                         because programme annotations are broadly enriched and
                         the ratio is confounded by annotation size.
    """
    df = pl.read_csv(path)
    if sig_source == "tau_z":
        z = pl.col("tau_z")
    else:
        z = (pl.col("enrichment") - 1.0) / pl.col("enrich_se")
    df = df.with_columns(z.alias("_z"))
    p = norm.sf(df["_z"].to_numpy())  # one-sided upper tail
    df = df.with_columns(pl.Series("p_value", p))
    fdr = multipletests(df["p_value"].to_numpy(), method="fdr_bh")[1]
    df = df.with_columns(pl.Series("fdr", fdr))
    return df


def program_labels(df):
    """program number -> display annotation (most common annotation wins)."""
    meta = {}
    for prog, sub in df.group_by("program"):
        prog = prog[0] if isinstance(prog, tuple) else prog
        ann = sub["annotation"].drop_nulls().mode()
        meta[int(prog)] = (ann[0] if len(ann) else f"programme {prog}")
    return meta


def select_programs(df_cond, max_programs, sig_p):
    """Programmes with >=1 cell at p < sig_p in this condition, ranked by max
    -log10(p); fall back to the best programmes if none reach sig_p."""
    best = (df_cond.group_by("program")
            .agg(pl.col("p_value").min().alias("p"))
            .sort("p"))
    sig = best.filter(pl.col("p") < sig_p)
    pool = sig if sig.height else best
    return [int(x) for x in pool["program"][:max_programs].to_list()]


def cluster_rows(M):
    """Leaf order + linkage for the disease rows from the value matrix M
    (studies x all programmes), correlation distance + average linkage."""
    if M.shape[0] < 3 or np.allclose(np.nan_to_num(M), 0.0):
        return None, None
    X = np.nan_to_num(M, nan=1.0)
    D = np.nan_to_num(pdist(X, metric="correlation"), nan=1.0)
    Z = linkage(D, method="average")
    return list(leaves_list(Z)), Z


def plot_condition(df, cond, studies, label_of, prog_meta, color, sig_source,
                   vmin, vmax, fdr_thr, max_programs, sig_p, do_cluster,
                   mask_nonsig, out_path):
    field, cbar_label = COLOR_SPECS[color]
    d = df.filter(pl.col("condition") == cond)
    if d.height == 0:
        print(f"  [{cond}] no rows; skipping")
        return
    cond_studies = [s for s in studies if s in set(d["study"].to_list())]

    # full (studies x all programmes) value matrix — used for clustering
    all_progs = sorted({int(p) for p in d["program"].to_list()})
    sidx = {s: i for i, s in enumerate(cond_studies)}
    pidx_all = {p: j for j, p in enumerate(all_progs)}
    full = np.full((len(cond_studies), len(all_progs)), np.nan)
    val_lookup, p_lookup, f_lookup = {}, {}, {}
    for r in d.iter_rows(named=True):
        if r["study"] not in sidx:
            continue
        i, j = sidx[r["study"]], pidx_all[int(r["program"])]
        full[i, j] = r[field]
        val_lookup[(r["study"], int(r["program"]))] = r[field]
        p_lookup[(r["study"], int(r["program"]))] = r["p_value"]
        f_lookup[(r["study"], int(r["program"]))] = r["fdr"]

    dendro_Z = None
    if do_cluster:
        perm, dendro_Z = cluster_rows(full)
        if perm is not None:
            cond_studies = [cond_studies[i] for i in perm]

    # displayed programmes (top by significance)
    programs = select_programs(d, max_programs, sig_p)

    # disambiguate duplicate annotation labels by appending the programme #
    seen = defaultdict(int)
    for p in programs:
        seen[prog_meta.get(p, "")] += 1
    prog_labels = []
    for p in programs:
        ann = prog_meta.get(p, f"programme {p}")
        ann = ann if len(ann) <= 28 else ann[:27] + "…"
        prog_labels.append(f"{ann} (P{p})" if seen[prog_meta.get(p, "")] > 1
                           else ann)

    shape = (len(cond_studies), len(programs))
    mat = np.full(shape, np.nan)
    pmat = np.full(shape, np.nan)
    fmat = np.full(shape, np.nan)
    for i, s in enumerate(cond_studies):
        for j, p in enumerate(programs):
            if (s, p) in val_lookup:
                mat[i, j] = val_lookup[(s, p)]
                pmat[i, j] = p_lookup[(s, p)]
                fmat[i, j] = f_lookup[(s, p)]

    if color == "enrichment":
        norm_obj = TwoSlopeNorm(vmin=vmin, vcenter=1.0, vmax=vmax)
        plot_mat = np.clip(mat, vmin, vmax)
    else:
        norm_obj = Normalize(vmin=-vmax, vmax=vmax)
        plot_mat = np.clip(mat, -vmax, vmax)

    if mask_nonsig:  # grey out cells that don't pass the FDR border
        plot_mat = np.where(np.isnan(fmat) | (fmat >= fdr_thr), np.nan, plot_mat)

    n_rows, n_cols = shape
    fig_w = max(9, 0.45 * n_cols + 4)
    fig_h = max(5, 0.34 * n_rows + 2.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#dddddd" if mask_nonsig else "white")
    im = ax.imshow(plot_mat, cmap=cmap, norm=norm_obj, aspect="auto")

    for i in range(n_rows):
        for j in range(n_cols):
            p = pmat[i, j]
            if np.isnan(p):
                continue
            stars = ("***" if p < 0.001 else "**" if p < 0.01
                     else "*" if p < 0.05 else "")
            if stars:
                t = ax.text(j, i, stars, ha="center", va="center",
                            color="white", fontsize=12, fontweight="bold")
                t.set_path_effects([path_effects.Stroke(linewidth=2,
                                    foreground="black"), path_effects.Normal()])
            if not np.isnan(fmat[i, j]) and fmat[i, j] < fdr_thr:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1.0, 1.0,
                             fill=False, edgecolor="#FFD700", linewidth=3.0))

    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([label_of[s] for s in cond_studies], fontsize=11)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(prog_labels, rotation=45, ha="right", fontsize=10)
    ax.tick_params(axis="both", length=0)
    for s in ax.spines.values():
        s.set_linewidth(0.5)

    cluster_note = ("\nrows clustered by enrichment profile (dendrogram, left)"
                    if dendro_Z is not None else "")
    sig_lab = "coefficient z" if sig_source == "tau_z" else "enrichment"
    ax.set_title(
        f"S-LDSC heritability enrichment — {cond}\n"
        f"{sig_lab} test: * p<0.05  ** p<0.01  *** p<0.001;  "
        f"gold border: BH-FDR < {fdr_thr}" + cluster_note,
        fontsize=13, fontweight="bold", loc="left", pad=8,
    )

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label(cbar_label, fontsize=12)
    if color == "enrichment":
        cbar.ax.axhline(cbar.norm(1.0), color="black", linewidth=0.8)

    left = 0.30 if dendro_Z is not None else 0.16
    fig.subplots_adjust(left=left, right=0.99, top=0.86, bottom=0.22)

    if dendro_Z is not None:
        pos = ax.get_position()
        dax = fig.add_axes([0.045, pos.y0, 0.085, pos.height])
        dendrogram(dendro_Z, orientation="left", ax=dax, no_labels=True,
                   color_threshold=0, above_threshold_color="#555555")
        dax.invert_yaxis()   # match imshow row 0 at top
        dax.axis("off")

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results", required=True,
                   help="all_results_parsed.csv from the S-LDSC zhu run.")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--conditions", nargs="+", default=COND_ORDER)
    p.add_argument("--studies", nargs="+", default=None,
                   help="Subset of studies; default = all in the table.")
    p.add_argument("--color", choices=list(COLOR_SPECS), default="tau_z",
                   help="Cell-colour metric. Default tau_z (the conditional "
                        "S-LDSC coefficient z) — it is centred at 0 and "
                        "discriminates disease-specific programmes. "
                        "'enrichment' shows the raw enrichment magnitude, which "
                        "is broadly positive and far less informative.")
    p.add_argument("--sig-source", choices=["tau_z", "enrichment"],
                   default="tau_z",
                   help="What drives the significance stars / FDR border and "
                        "programme selection. Default tau_z (the conditional "
                        "S-LDSC coefficient test) — far more disease-specific "
                        "than the broadly-positive enrichment ratio.")
    p.add_argument("--vmax", type=float, default=None,
                   help="Upper colour cap (enrichment) or symmetric cap (tau_z). "
                        "Default: 4 for tau_z, 15 for enrichment.")
    p.add_argument("--vmin", type=float, default=-5.0,
                   help="Lower colour cap for enrichment (ignored for tau_z).")
    p.add_argument("--mask-nonsig", action=argparse.BooleanOptionalAction,
                   default=False,
                   help="Grey out cells that do not pass the BH-FDR border, so "
                        "only the statistically supported signal is coloured.")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--sig-p", type=float, default=0.05,
                   help="Nominal p a programme must reach in >=1 cell to be "
                        "shown as a column.")
    p.add_argument("--max-programs", type=int, default=40)
    p.add_argument("--cluster-rows", action=argparse.BooleanOptionalAction,
                   default=True,
                   help="Cluster disease rows by enrichment profile and draw an "
                        "inline dendrogram (default on).")
    p.add_argument("--format", default="png", choices=["png", "pdf"])
    args = p.parse_args()

    vmax = args.vmax if args.vmax is not None else (
        4.0 if args.color == "tau_z" else 15.0)

    df = load_results(args.results, args.sig_source)
    if args.studies:
        df = df.filter(pl.col("study").is_in(args.studies))
    studies = sorted(df["study"].unique().to_list())
    label_of = make_short_labels(studies)
    prog_meta = program_labels(df)
    n_sig = int((df["fdr"] < args.fdr).sum())
    print(f"studies: {len(studies)}  programmes: {df['program'].n_unique()}  "
          f"tests: {df.height}  passing BH-FDR<{args.fdr}: {n_sig}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for cond in args.conditions:
        out_path = out_dir / f"ldsc_program_enrichment_{cond}.{args.format}"
        plot_condition(
            df, cond, studies, label_of, prog_meta, args.color, args.sig_source,
            args.vmin, vmax, args.fdr, args.max_programs, args.sig_p,
            args.cluster_rows, args.mask_nonsig, out_path,
        )

    csv_out = out_dir / "ldsc_program_enrichment.csv"
    df.sort(["condition", "study", "p_value"]).write_csv(csv_out)
    print(f"wrote {csv_out}")


if __name__ == "__main__":
    main()
