# Heritability enrichment (Fig. 4, S5–S7) — notes

Variant **coloc_thr0.5** (anchors: CD4 coloc PP.H4 > 0.8; edges |logFC| ≥ 0.5).
Scripts: `chapter2/vis/heritability/` (run from `chapter2/vis`, e.g.
`uv run heritability/01_collect_results.py`; `uv run 00_build_all.py` runs all).
Nothing here re-runs S-LDSC or touches gene sets, annotations or results.

---

## Step 0 — what is in `data/ldsc_results/coloc_thr0.5/`

### 0.1 File layout

One folder per study; per condition `{Rest,Stim8hr,Stim48hr}_combined.{results,log,cov}`.
52 `.log`, 49 `.results`, 49 `.cov` (no `.part_delete` or other files).
Per-run inventory with files present: `data/heritability/ldsc_runs.csv`.

| Study | Rest | 8 h | 48 h |
|---|---|---|---|
| 15 studies (all but the three below) | results, log, cov | results, log, cov | results, log, cov |
| ADFH_Jansen_30617256_1 (Alzheimer's proxy) | — (no files) | — (no files) | results, log, cov |
| SCZ_Pardiñas_29483656_1 (schizophrenia) | log only (failed) | log only (failed) | log only (failed) |

### 0.2 Model

- ldsc v3.0.1, `--h2 … --overlap-annot --print-coefficients --print-cov`.
- `--ref-ld-chr baselineLD_v2.2/baselineLD.,<study>/<cond>_combined.`: **baselineLD
  v2.2 (97 annotations incl. base) plus the custom annotations, all in one model**.
  Anchors, regulators and targets are fitted **jointly** in the same regression per
  trait × condition, so each set's τ is conditional on the baseline and on the other
  two sets.
- Category names in `.results`: `H4_genesL2_1` (anchors), `regulators_H4L2_1`,
  `targets_H4L2_1`; baseline categories end in `L2_0`.
- Annotations: binary, gene body ± 10 kb (`window_size: 10000` in
  `perturb_ldsc_pipeline/config.yaml`).
- **Single unnamed column `L2_1`.** `combine_annots.py` drops any set covering
  < 1,000 SNPs (summed over chr 1–22) before LD scores are computed. When only one
  set survives, ldsc names its LD-score column plain `L2`, so the category loses its
  set name. This happens in 4 runs; mapped explicitly in `01_collect_results.py`
  (`SINGLE_COLUMN`, asserted):
  - Osteoarthritis, all conditions → **regulators**. Anchors (5 genes) dropped,
    targets empty. Check: Prop._SNPs follows the regulator count (108 / 111 / 49
    genes → 0.22 / 0.21 / 0.09 %).
  - Alzheimer's proxy, 48 h → **targets**.
    - Spans (gene ± 10 kb, Ensembl 115): anchors ≈ 141 kb, regulators ≈ 84 kb,
      targets ≈ 1.0 Mb. 1000G has about 1 SNP per 300 bp, so only the targets
      clear 1,000 SNPs.
    - The same anchors and the 8 h regulators (≈ 221 kb, 5 genes) were dropped in
      the other conditions, which left nothing to fit at Rest and 8 h. That is why
      those runs have no files.
  - This mapping is an inference from the gene sets: the `.annot` files are on the
    cluster. Confirm with `head -1 …/ADFH*/Stim48hr_combined.1.annot` if needed.

### 0.3 Missing runs and empty annotations

Status of the 162 trait × condition × annotation cells (`status` column):

| Status | n | Cells |
|---|---|---|
| fitted | 129 | |
| empty set | 15 | targets: Asthma (EAS), Ulcerative colitis, Ischaemic stroke, Osteoarthritis (all conditions); Alzheimer's proxy (Rest, 8 h); Type 2 diabetes (Rest) |
| dropped (< 1,000 SNPs) | 5 | anchors: Osteoarthritis (all conditions), Alzheimer's proxy (48 h); regulators: Alzheimer's proxy (48 h) |
| not run | 4 | Alzheimer's proxy Rest / 8 h anchors and regulators (all sets < 1,000 SNPs → no run) |
| run failed | 9 | Schizophrenia, all conditions and annotations |

Schizophrenia failure (all three logs):
`ValueError: After merging with reference panel LD, 0 SNPs remain.`
None of the sumstats SNP IDs matched the GRCh37 reference. This is likely a build
or SNP-ID format problem in `SCZ_Pardiñas_29483656_1.sumstats.gz`. It needs a re-run
on the cluster; it was not attempted here.

### 0.4 τ*

**Not computed anywhere in the repo before this.** Earlier scripts
(`perturb_programme_heritability/6_visualise_ldsc_zhu.py`, `7_…`) label the
coefficient z-score "tau*", but they never standardise τ. They are not reused.

**Formula** (Gazal et al. 2017):

    τ*  = τ · sd(a) · M / h²g           sd(a) = sqrt(p (1 − p)),  p = Prop._SNPs
    SE(τ*) = SE(τ) · sd(a) · M / h²g    (h²g, M treated as fixed)
    z = τ / SE(τ) (= Coefficient_z-score),  p = Φ(−z)  (one-sided, H1: τ* > 0)

The two-sided p = 2 Φ(−|z|) and its FDR/significance are kept in
`heritability_long.csv` as `tau_star_p_two_sided`, `tau_star_fdr_two_sided`,
`tau_star_sig_two_sided`.

**Inputs.**
- τ, SE(τ): `Coefficient`, `Coefficient_std_error` from `.results`.
- p: `Prop._SNPs`. For these binary annotations this is the fraction of the
  M_5_50 reference SNPs in the set.
- M / h²g: the `.l2.M_5_50` files are not stored locally (they are in
  `data/sensitivity/thr0.5/ld_scores/` on the cluster), so M / h²g is taken from
  the `.results` file itself. The derivation uses ldsc's code
  (`ldscore/regressions.py`, `Hsq._cat`, `_tot`, `_overlap_output`):
  - ldsc computes h²g = Σₖ τₖ Mₖ and writes Prop._SNPsₖ = Mₖ / M, summing over
    all categories (baseline + custom; M = base-annotation M_5_50).
  - Hence **h²g / M = Σₖ τₖ · Prop._SNPsₖ**, and τ* = τ · sd(a) / Σₖ τₖ pₖ.
  - This is exact and avoids the 4-significant-figure rounding of h²g in the log.

**Checks.**
- *Implied M.* M = h²g(log) / Σₖ τₖ pₖ comes out at 5.95–5.97 million for every
  GRCh37 EUR run. The published M_5_50 for baselineLD v2.2 (1000G EUR) is
  5,961,159.
  - The spread is the rounding of h²g in the log. Ischaemic stroke is the largest
    deviation (6.00 M, h²g = 0.0064).
  - EAS runs give ≈ 5.47 M; the GRCh38 EUR runs give ≈ 5.95 M.
- *Hand check, childhood-onset asthma, Rest, regulators.*
  - Inputs: τ = 4.0488e-8, SE = 1.8236e-8, p = 0.0046623, so sd = 0.068122.
    M = 5,961,159 (published); h²g = 0.0734 (log).
  - By hand: τ* = 4.0488e-8 × 0.068122 × 5,961,159 / 0.0734 = **0.2240**,
    SE = 0.1009, z = 2.220, p = 0.0264.
  - Script: τ* = 0.2241, SE = 0.1009, z = 2.220, two-sided p = 0.0264 (one-sided
    0.0132). They agree to 0.05 % (the rounding of h²g).

**Enrichment p and z.**
- `Enrichment_p` is ldsc's own test: per-SNP h² inside vs outside the annotation,
  two-sided t with 200 jackknife blocks.
- The signed "enrichment z" on Fig. S7 is that p converted back with the t(200)
  distribution and given the sign of (enrichment − 1). It matches the p exactly.

**FDR.** Benjamini–Hochberg within two families (`fdr_family`):
- **expanded sets**: regulators + targets, all traits × conditions (84 tests);
- **anchors**: 45 tests.

Each family is corrected separately for enrichment p (two-sided) and τ* p
(one-sided). Significance is in the `enrichment_sig` and `tau_star_sig` columns
(q < 0.05).

**Negative Prop._h2.** Flagged in `neg_prop_h2` and never counted as significant
in either direction. Such sets stay in the BH correction (as one test each).
There is only one: **Systemic lupus erythematosus, Rest, regulators**
(Prop._h2 = −0.0047, enrichment = −11.8, enrichment p = 8.5e-5, τ* z = −4.26).
Under the old scheme it was "significant" on both enrichment and τ*.

### 0.5 LD reference and ancestry

The reference is chosen per study (`--ref-ld-chr`, `--w-ld-chr`, `--frqfile-chr` all
from the same panel):

| Panel | Traits |
|---|---|
| 1000G EUR, GRCh37 | 13 traits, incl. schizophrenia (failed) |
| 1000G EUR, GRCh38 | Autoimmune thyroid disease, Type 1 diabetes |
| 1000G **EAS**, GRCh37 | Asthma (EAS), Rheumatoid arthritis (EAS) |

- **The brief expected RA (EAS) and asthma (EAS) to be mismatched. They are not:**
  both were run with an East Asian baselineLD, weights and frequencies.
- The real mismatches are the **five multi-ancestry GWAS run with the EUR
  reference**: Type 1 diabetes, Ischaemic stroke, Myocardial infarction,
  Osteoarthritis, Type 2 diabetes (ancestry from `diseases.csv`).
- The **asterisk (*)** on heritability figures marks these five, not the EAS traits.
  The rule is "GWAS ancestry ≠ LD reference ancestry", computed in
  `01_collect_results.py`. To flag the EAS traits too, change `flagged_name` in
  `hcommon.py`.
- Note that Fig. 2.1B also uses * for a different meaning ("no targets in any
  condition"). Each caption defines its own *.

---

## Measured results

All counts are from `data/numbers_for_text.md`. Significant = q < 0.05 in its BH
family (expanded / anchors), with one-sided τ* and negative-Prop._h2 sets excluded.

- **Anchors.** Significant in 7/12 IMDs and 2/3 controls with anchor estimates,
  on both enrichment and τ* (same 27/45 sets for both).
  - Median enrichment: 78 (IMD), 68 (control).
  - Not significant: Asthma (EAS), Crohn's disease, Ulcerative colitis, RA (EAS),
    SLE and Ischaemic stroke. These have wide SEs.
- **Regulators, IMDs:**
  - Enrichment significant in 7/12 (Rest), 8/12 (8 h) and 10/12 (48 h).
  - τ* significant in 0/12 in every condition.
  - Median enrichment 5.2 (log2 2.39).
- **Regulators, controls:** 1/12 enriched (Ischaemic stroke at 8 h), 0/12 τ*.
  Median enrichment 1.4 (log2 0.45).
- **Targets, IMDs:** enrichment 6/10, 8/10, 8/10; τ* 0/10 throughout.
  Median enrichment 3.3.
- **Targets, controls:** 1/6 enriched (Type 2 diabetes at 8 h), 0/6 τ*.
- **Expanded sets overall:**
  - IMD: 47/66 enriched, 0/66 τ*-significant, so all 47 are enriched but not
    τ*-significant.
  - Controls: 2/18 enriched, 0/18 τ*.
- **Smallest τ* q among expanded sets:** 0.057, shared by Autoimmune thyroid
  disease 48 h regulators, RA (EAS) 8 h regulators and SLE 8 h targets. These are
  close, but none passes.
- **Enrichment vs size (Fig. S6).** Enrichment falls with annotation size:
  - anchors cover 0.01–0.26 % of reference SNPs, with log2 enrichment 4.4–8;
  - regulators cover 0.03–1.1 %, with log2 enrichment −1.5 to 3.7;
  - targets cover 0.03–6 %, with log2 enrichment 0.2–3.1.
- **FDR 5% hurdles for expanded sets:**
  - enrichment |z| ≥ 2.19;
  - τ* z ≥ 3.24. No set passes, so this is the BH rank-1 hurdle, p ≤ 0.05/84.
  - The largest expanded τ* z is 2.97.

### How the follow-up changed the significant counts

| Scheme | Enrichment, expanded IMD / control | τ*, expanded IMD / control | Anchors (enr. / τ*) |
|---|---|---|---|
| Before: two-sided τ*, one BH over all 129 sets | 48/66 / 2/18 | **7/66** / 0/18 | 27 / 27 |
| BH by family, two-sided τ* | 48/66 / 2/18 | 0/66 / 0/18 | 27 / 27 |
| **Now:** BH by family, one-sided τ*, negative Prop._h2 excluded | **47/66** / 2/18 | **0/66** / 0/18 | 27 / 27 |

- **τ* 7 → 0 is caused by the family split, not by the one-sided test.**
  - Previously, the 27 very small anchor p-values were in the same BH correction
    and let the expanded sets through at q < 0.05.
  - Within the 84 expanded tests alone, the 6 positive sets have q = 0.057–0.074
    (two-sided gives none either).
  - One-sided halves those six p-values, but that still is not enough.
- **The 7th set (SLE regulators at Rest, negative τ*)** drops out three ways:
  one-sided p ≈ 1, it is flagged for negative Prop._h2, and it is no longer in
  the anchors' family.
- **Enrichment 48 → 47** is the same SLE set, now excluded for negative Prop._h2.
  The family split changes no enrichment call.
- **Anchors are unchanged** (27/45 on both).

## Interpretation

- **The brief's message holds, and more strongly than before.** Regulator and
  target sets carry more IMD heritability per SNP than the genome average, and
  more than in controls. Once baselineLD v2.2 and the anchors are in the model,
  none of the 47 enriched IMD expanded sets has a significant τ*.
  - Six positive τ* sets sit just above the threshold (q ≈ 0.06–0.07). They
    should be described as suggestive at most.
  - The marginal enrichment is therefore largely explained by what the expanded
    genes have in common with other functional annotations (e.g. promoter/enhancer
    density, conservation and eQTL annotations in baselineLD), and by LD with the
    anchors.
  - It is not a separate heritability component unique to the perturbation
    expansion.
- **Anchors are not an independent positive control.** Anchors are selected on the
  same GWAS signal (coloc PP.H4 > 0.8), so their enrichment is expected by
  construction. They are useful as a scale reference, not as validation.
- **Regulators look slightly more specific than targets.** IMD regulator sets are
  more enriched than control regulator sets (median 5.2 vs 1.4). For targets the
  IMD/control gap is small, though control targets are few.
- **Enrichment depends on annotation size (Fig. S6).** Comparing enrichment across
  set types mostly compares sizes; τ* is the size-adjusted comparison.
- The IMD/control comparison is unbalanced: 12 IMDs against 4–5 usable controls.
  Two control traits contribute no regulators or targets (schizophrenia failed;
  Alzheimer's proxy is mostly not run), so control counts are small.

## Not done / open

- **Schizophrenia S-LDSC failed** (0 SNPs after merging with the reference). The
  sumstats need checking (build / SNP IDs) and a re-run on the cluster.
- **Alzheimer's proxy at Rest and 8 h, and the osteoarthritis anchors, are absent**
  by design: the 1,000-SNP filter (`min_snps_per_annot`).
- **`L2_1` mapping** for osteoarthritis and Alzheimer's proxy 48 h is inferred,
  not read from the `.annot` header (not available locally); see 0.2.
- **M_5_50 not read from files.** τ* uses the exact h²g/M identity from `.results`
  instead (0.4). If you want it read from the files, copy
  `data/sensitivity/thr0.5/ld_scores/*/*.l2.M_5_50` from the cluster.
- **Fig. 4C** (enrichment with an all-measured-genes annotation added) is reserved
  in INDEX.md. S-LDSC for it has not been run.
- **Fig. 2.2C** now reads coverage from `heritability_long.csv`, so it includes
  the four `L2_1` sets (osteoarthritis regulators, Alzheimer's 48 h targets).
  `heritability/01_collect_results.py` therefore runs first in `00_build_all.py`.
- No new statistical tests beyond BH FDR.
