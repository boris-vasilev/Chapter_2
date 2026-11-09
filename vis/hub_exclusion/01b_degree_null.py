#!/usr/bin/env python
"""Degree-aware expected Jaccard for regulator sets (fig_h1 E).

The repo has no degree-preserving overlap null, so one is built here. For
each variant x condition, the regulator universe (perturbations with >= 1
edge in that variant's hub-filtered graph) is split into out-degree bins
(out-degree = number of distinct targets in the hub-filtered graph): 20
quantile bins of log out-degree, adjacent bins merged until every trait's
count in a bin fits in it. Each null draw replaces every trait's regulator
set by a random set with the same count per bin, drawn independently per
trait, so set sizes and out-degree profiles are kept and the trait-to-trait
dependence is removed. B = 1,000 draws, seed fixed per variant x condition.

Per pair: expected_deg = mean null Jaccard; excess_deg = observed - expected.
Per pair class (IMD-IMD, IMD-control, control-control; related-GWAS pairs
excluded; pairs involving an empty set excluded; the related pairs get the same
summary but are not in the BH family): statistic = median over
pairs of excess_deg; null = the same median computed on each draw;
empirical p (one-sided, H1: observed > null) = (1 + #{null >= obs}) / (B + 1).
BH within each variant across the 9 class x condition tests (as the other
panels: within each variant, all tests in a panel).

Writes data/jaccard_degree_null_pairs.csv and data/h1e_degree_null_summary.csv.
"""

import numpy as np
import pandas as pd
import polars as pl
from statsmodels.stats.multitest import multipletests

from hxcommon import (CONDS, DOUT, PAIR_CLASSES, PCT, VKEYS, edges, gene_sets,
                      load_diseases, pair_class, studies)

B = 1000
N_BINS = 20
ANNOT = "regulators_H4"


def degree_bins(deg, counts):
    """Bin index per universe gene: quantile bins of log degree, merged until
    each trait's count per bin <= bin size. ``counts``: (T, N) membership."""
    q = np.unique(np.quantile(np.log(deg), np.linspace(0, 1, N_BINS + 1)[1:-1]))
    b = np.searchsorted(q, np.log(deg), side="right")
    while True:
        labels = np.unique(b)
        size = np.array([(b == l).sum() for l in labels])
        need = np.array([[counts[t, b == l].sum() for l in labels] for t in range(len(counts))])
        bad = np.where((need > size).any(axis=0))[0]
        if not len(bad):
            return np.searchsorted(labels, b)
        i = bad[-1]                         # merge with its lower neighbour
        lo = labels[i - 1] if i > 0 else labels[i + 1]
        b[b == labels[i]] = lo


def null_draws(member, bins, rng):
    """Generator of (T, N) boolean membership matrices, B draws."""
    T, N = member.shape
    nb = bins.max() + 1
    order0 = np.argsort(bins, kind="stable")
    sorted_bins = bins[order0]
    start = np.searchsorted(sorted_bins, np.arange(nb))
    rank = np.arange(N) - start[sorted_bins]               # rank within bin
    k = np.stack([np.bincount(bins[m], minlength=nb) for m in member])  # (T, nb)
    take = rank[None, :] < k[:, sorted_bins]               # (T, N) in bin-sorted space
    for _ in range(B):
        # random order within each bin: sort by bin + uniform key
        keys = sorted_bins[None, :] + rng.random((T, N))
        perm = np.argsort(keys, axis=1)                    # positions in bin-sorted space
        out = np.zeros((T, N), bool)
        np.put_along_axis(out, order0[perm], take, axis=1)
        yield out


def jaccard(m):
    m = m.astype(np.float32)
    inter = m @ m.T
    size = m.sum(1)
    union = size[:, None] + size[None, :] - inter
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(union > 0, inter / union, 0.0)


def main():
    diseases = load_diseases()
    st = studies()
    iu = np.triu_indices(len(st), 1)
    cls = np.array([pair_class(diseases, st[i], st[j]) for i, j in zip(*iu)])
    pair_rows, summ = [], []
    for v in VKEYS:
        for ci, c in enumerate(CONDS):
            e = edges(v).filter(pl.col("culture_condition") == c)
            deg = dict(e.group_by("cis_ensg").agg(pl.col("trans_ensg").n_unique()).iter_rows())
            genes = sorted(deg)
            idx = {g: i for i, g in enumerate(genes)}
            d = np.array([deg[g] for g in genes], float)
            member = np.zeros((len(st), len(genes)), bool)
            for t, s in enumerate(st):
                member[t, [idx[g] for g in gene_sets(v)[(s, c, ANNOT)]]] = True
            bins = degree_bins(d, member)
            rng = np.random.default_rng(42 + 10 * VKEYS.index(v) + ci)
            obs = jaccard(member)[iu]
            draws = np.array([jaccard(m)[iu] for m in null_draws(member, bins, rng)])
            exp = draws.mean(0)
            empty = member.sum(1) == 0
            ok = ~(empty[iu[0]] | empty[iu[1]])
            exc = np.where(ok, obs - exp, np.nan)
            for p, (i, j) in enumerate(zip(*iu)):
                pair_rows.append({"variant": v, "pct_excluded": PCT[v], "condition": c,
                                  "study_a": st[i], "study_b": st[j], "pair_class": cls[p],
                                  "observed": obs[p], "expected_degree": exp[p] if ok[p] else np.nan,
                                  "excess_degree": exc[p]})
            for k in PAIR_CLASSES + ["related"]:
                sel = (cls == k) & ok
                stat = np.median(obs[sel] - exp[sel])
                null = np.median(draws[:, sel] - exp[sel], axis=1)
                summ.append({"variant": v, "pct_excluded": PCT[v], "condition": c,
                             "pair_class": k, "n_pairs": int(sel.sum()),
                             "median_excess_degree": stat,
                             "q25": np.quantile(obs[sel] - exp[sel], 0.25),
                             "q75": np.quantile(obs[sel] - exp[sel], 0.75),
                             "null_median_q95": np.quantile(null, 0.95),
                             "p_empirical": (1 + (null >= stat).sum()) / (B + 1),
                             "n_draws": B, "n_degree_bins": int(bins.max() + 1)})
            print(f"{v} {c}: {bins.max() + 1} degree bins", flush=True)
    s = pd.DataFrame(summ)
    s["fdr"] = np.nan
    for v, g in s[s.pair_class != "related"].groupby("variant"):
        s.loc[g.index, "fdr"] = multipletests(g.p_empirical, method="fdr_bh")[1]
    s["sig"] = s.fdr < 0.05
    pd.DataFrame(pair_rows).to_csv(DOUT / "jaccard_degree_null_pairs.csv", index=False)
    s.to_csv(DOUT / "h1e_degree_null_summary.csv", index=False)
    print(f"Wrote {DOUT / 'h1e_degree_null_summary.csv'}")


if __name__ == "__main__":
    main()
