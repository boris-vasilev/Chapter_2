#!/usr/bin/env python
"""Collect every table behind the hub-exclusion figures (data/).

  hub_regulators.csv      hubs per threshold x condition (gene, n_degs, anchor?)
  gene_set_rebuild.csv    Check 1: every regulator / target / anchor set rebuilt
                          from the hub-filtered graph vs gene_sets/ (must be 0
                          mismatches; the script stops otherwise)
  set_sizes.csv           set sizes, universes, anchors lost as hubs
  jaccard_pairs.csv       pairwise observed / expected / excess Jaccard
  heritability_long.csv   S-LDSC estimates for all variants (01_collect_results
                          parsers; BH-FDR within each variant as in Fig. 4)
  anchor_shifts.csv       anchor estimates vs deg0 (replaces the identity check:
                          anchors are fitted jointly with regulators + targets)
  top_shared_regulators.csv  top 20 regulators by number of IMD sets, per
                          variant x condition
  hub_anchor_fraction.csv fraction of perturbed anchors that are hubs vs all
                          perturbations

Target sets are flagged when the variant removed >= 1 of their seeding anchors
(lost_anchor); ``unchanged_10`` marks sets whose seeding anchors are all
kept even at 10 % exclusion (a fixed subset, nested across thresholds since
the hub lists are nested).
"""

import itertools

import numpy as np
import pandas as pd
import polars as pl

from hxcommon import (CONDS, FOLDER, PCT, QUANTILE, SET_ANNOTS, VKEYS, anchors,
                      collect_module, edges, gene_sets, hubs, ldsc_dir,
                      load_diseases, lost_anchors, obs_table, out_dirs,
                      pair_class, perturbed_anchors, studies, universe, DOUT)
import importlib.util
from pathlib import Path

# expected_jaccard / jaccard_matrix from 04_jaccard_heatmaps.py
_p = Path(__file__).resolve().parents[1] / "04_jaccard_heatmaps.py"
_s = importlib.util.spec_from_file_location("jacc", _p)
jacc = importlib.util.module_from_spec(_s)
_s.loader.exec_module(jacc)


def hub_table(diseases):
    sym = dict(obs_table().select("ensg", "symbol").iter_rows())
    nd = {(c, g): n for c, g, n in obs_table().select(
        "culture_condition", "ensg", "n_degs").iter_rows()}
    anc_of = {}
    for s, genes in anchors().items():
        if s in diseases.index:
            for g in genes:
                anc_of.setdefault(g, []).append(diseases.loc[s, "display_name"])
    rows = [{"gene": sym.get(g, g), "ensg": g, "condition": c, "n_degs": nd[(c, g)],
             "threshold": QUANTILE[v], "pct_excluded": PCT[v],
             "is_anchor": g in anc_of, "anchor_traits": "; ".join(anc_of.get(g, []))}
            for v in VKEYS[1:] for c in CONDS for g in hubs(v, c)]
    return pd.DataFrame(rows).sort_values(["pct_excluded", "condition", "n_degs"],
                                          ascending=[True, True, False])


def rebuild(variant, study, cond, annot):
    """make_gene_sets_coloc.py logic on the hub-filtered graph."""
    e = edges(variant).filter(pl.col("culture_condition") == cond)
    anc = list(anchors().get(study, ()))
    if annot == "regulators_H4":
        return set(e.filter(pl.col("trans_ensg").is_in(anc))["cis_ensg"])
    if annot == "targets_H4":
        return set(e.filter(pl.col("cis_ensg").is_in(anc))["trans_ensg"])
    return set(anc)


def check_rebuild():
    rows = []
    for v in VKEYS:
        for (s, c, a), g in gene_sets(v).items():
            r = rebuild(v, s, c, a)
            rows.append({"variant": v, "study": s, "condition": c, "annot": a,
                         "n_published": len(g), "n_rebuilt": len(r), "match": r == g})
    df = pd.DataFrame(rows)
    bad = df[~df["match"]]
    if len(bad):
        print(bad.to_string())
        raise SystemExit(f"Check 1 failed: {len(bad)} gene sets differ from gene_sets/")
    print(f"Check 1: {len(df)} sets rebuilt, 0 mismatches")
    return df


def set_sizes(diseases):
    rows = []
    for v in VKEYS:
        for c in CONDS:
            regs, tgts = universe(edges(v), c)
            for s in studies():
                lost = lost_anchors(v, s, c)
                for a in SET_ANNOTS:
                    n_u = len(regs) if a == "regulators_H4" else len(tgts)
                    n = len(gene_sets(v)[(s, c, a)])
                    rows.append({"variant": v, "pct_excluded": PCT[v], "study": s,
                                 "disease": diseases.loc[s, "display_name"],
                                 "group": str(diseases.loc[s, "group"]),
                                 "condition": c, "annot": a, "n_genes": n,
                                 "universe": n_u, "pct_of_universe": 100 * n / n_u,
                                 "n_perturbed_anchors": len(perturbed_anchors(v, s, c)),
                                 "n_lost_anchors": len(lost),
                                 "lost_anchor": a == "targets_H4" and bool(lost)})
    df = pd.DataFrame(rows)
    # Fixed subset: target sets with no seeding anchor lost even at 10 %.
    lost10 = df[(df.variant == VKEYS[-1]) & (df.annot == "targets_H4")].set_index(
        ["study", "condition"])["lost_anchor"]
    df["unchanged_10"] = [a != "targets_H4" or not lost10[(s, c)]
                          for s, c, a in zip(df.study, df.condition, df.annot)]
    return df


def jaccard_pairs(diseases, sizes):
    st = studies()
    unch = sizes.set_index(["variant", "study", "condition", "annot"])["unchanged_10"]
    rows = []
    for v in VKEYS:
        for c in CONDS:
            regs, tgts = universe(edges(v), c)
            for a in SET_ANNOTS:
                sets = [gene_sets(v)[(s, c, a)] for s in st]
                obs = jacc.jaccard_matrix(sets)
                n_u = len(regs) if a == "regulators_H4" else len(tgts)
                exp = jacc.expected_jaccard([len(x) for x in sets], n_u)
                for i, j in itertools.combinations(range(len(st)), 2):
                    rows.append({
                        "variant": v, "pct_excluded": PCT[v], "condition": c, "annot": a,
                        "study_a": st[i], "study_b": st[j],
                        "disease_a": diseases.loc[st[i], "display_name"],
                        "disease_b": diseases.loc[st[j], "display_name"],
                        "pair_class": pair_class(diseases, st[i], st[j]),
                        "observed": obs[i, j], "expected": exp[i, j],
                        "excess": obs[i, j] - exp[i, j],
                        "both_unchanged_10": bool(unch[(v, st[i], c, a)]
                                                  and unch[(v, st[j], c, a)])})
    return pd.DataFrame(rows)


def heritability(diseases, sizes):
    coll = collect_module()
    out = []
    for v in VKEYS:
        long, _ = coll.collect(ldsc_dir(FOLDER[v]) / "ldsc_results", gene_sets(v), diseases)
        long = coll.add_fdr(long)   # BH within this variant, families as Fig. 4
        long.insert(0, "variant", v)
        long.insert(1, "pct_excluded", PCT[v])
        out.append(long)
    df = pd.concat(out, ignore_index=True)
    df.insert(3, "disease", df["study"].map(diseases["display_name"]))
    df.insert(4, "group", df["study"].map(diseases["group"].astype(str)))
    df["log2_enrichment"] = np.log2(df["enrichment"].where(df["enrichment"] > 0))
    # SE of log2 enrichment by the delta method: se(E) / (E ln 2).
    df["log2_enrichment_se"] = df["enrichment_se"] / (df["enrichment"].where(
        df["enrichment"] > 0) * np.log(2))
    flags = sizes.set_index(["variant", "study", "condition", "annot"])[
        ["lost_anchor", "n_lost_anchors", "unchanged_10"]]
    df = df.join(flags, on=["variant", "study", "condition", "annot"])
    anchor = df["annot"] == "H4_genes"
    df.loc[anchor, ["lost_anchor", "unchanged_10"]] = [False, True]
    df.loc[anchor, "n_lost_anchors"] = 0
    return df


def anchor_shifts(h):
    a = h[h.annot == "H4_genes"].set_index(["study", "condition"])
    base = a[a.variant == "deg0"]
    rows = []
    for v in VKEYS[1:]:
        x = a[a.variant == v].reindex(base.index)
        d = pd.DataFrame({
            "variant": v, "pct_excluded": PCT[v], "disease": base["disease"],
            "status": x["status"], "prop_snps_diff": x["prop_snps"] - base["prop_snps"],
            "enrichment_deg0": base["enrichment"], "enrichment": x["enrichment"],
            "enrichment_shift_in_se": (x["enrichment"] - base["enrichment"]) / base["enrichment_se"],
            "tau_star_deg0": base["tau_star"], "tau_star": x["tau_star"],
            "tau_star_shift_in_se": (x["tau_star"] - base["tau_star"]) / base["tau_star_se"],
            "enrichment_sig_deg0": base["enrichment_sig"], "enrichment_sig": x["enrichment_sig"],
            "tau_star_sig_deg0": base["tau_star_sig"], "tau_star_sig": x["tau_star_sig"]})
        rows.append(d.reset_index())
    return pd.concat(rows, ignore_index=True)


def top_shared(diseases):
    sym = dict(obs_table().select("ensg", "symbol").iter_rows())
    imd = [s for s in studies() if diseases.loc[s, "group"] == "IMD"]
    ctl = [s for s in studies() if diseases.loc[s, "group"] != "IMD"]
    rows = []
    for v in VKEYS:
        for c in CONDS:
            deg = edges(v).filter(pl.col("culture_condition") == c).group_by(
                "cis_ensg").len()
            deg = dict(deg.iter_rows())
            cnt = {}
            for grp, ss in [("n_imd", imd), ("n_control", ctl)]:
                for s in ss:
                    for g in gene_sets(v)[(s, c, "regulators_H4")]:
                        cnt.setdefault(g, {"n_imd": 0, "n_control": 0})[grp] += 1
            df = pd.DataFrame([{"gene": sym.get(g, g), "ensg": g, **k,
                                "out_degree": deg.get(g, 0)} for g, k in cnt.items()])
            df = df.sort_values(["n_imd", "n_control", "gene"],
                                ascending=[False, True, True]).head(20)
            df.insert(0, "rank", range(1, len(df) + 1))
            df.insert(0, "condition", c)
            df.insert(0, "pct_excluded", PCT[v])
            df.insert(0, "variant", v)
            rows.append(df)
    return pd.concat(rows, ignore_index=True)


def hub_anchor_fraction(diseases):
    """Per threshold x condition: hubs among perturbed IMD anchors, control
    anchors and all perturbations (all with >= 1 edge in deg0)."""
    groups = {"IMD": [s for s in studies() if diseases.loc[s, "group"] == "IMD"],
              "control": [s for s in studies() if diseases.loc[s, "group"] != "IMD"]}
    rows = []
    for v in VKEYS[1:]:
        for c in CONDS + ["all"]:
            cs = CONDS if c == "all" else [c]
            r = {"variant": v, "pct_excluded": PCT[v], "condition": c}
            for name, pool in [("IMD_anchors", groups["IMD"]),
                               ("control_anchors", groups["control"]), ("all", None)]:
                n = k = 0
                for cc in cs:
                    pert = set(edges("deg0").filter(
                        pl.col("culture_condition") == cc)["cis_ensg"])
                    genes = pert if pool is None else pert & set().union(
                        *[anchors().get(s, set()) for s in pool])
                    n += len(genes)
                    k += len(genes & hubs(v, cc))
                r[f"n_{name}"], r[f"hubs_{name}"] = n, k
                r[f"frac_{name}"] = k / n if n else np.nan
            rows.append(r)
    return pd.DataFrame(rows)


def main():
    out_dirs()
    diseases = load_diseases()
    hub_table(diseases).to_csv(DOUT / "hub_regulators.csv", index=False)
    check_rebuild().to_csv(DOUT / "gene_set_rebuild.csv", index=False)
    sizes = set_sizes(diseases)
    sizes.to_csv(DOUT / "set_sizes.csv", index=False)
    jaccard_pairs(diseases, sizes).to_csv(DOUT / "jaccard_pairs.csv", index=False)
    h = heritability(diseases, sizes)
    h.to_csv(DOUT / "heritability_long.csv", index=False)
    anchor_shifts(h).to_csv(DOUT / "anchor_shifts.csv", index=False)
    top_shared(diseases).to_csv(DOUT / "top_shared_regulators.csv", index=False)
    hub_anchor_fraction(diseases).to_csv(DOUT / "hub_anchor_fraction.csv", index=False)
    print(f"Wrote tables to {DOUT}")


if __name__ == "__main__":
    main()
