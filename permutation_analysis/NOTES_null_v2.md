# Corrected matched-null S-LDSC permutation (Chapter 2) — notes

Session 2026-10-28. Code: `perturb_heritability_null_v2/` (copy of
`perturb_ldsc_permutation/`; the pilot folders and `data/ldsc_permutation*`
were not modified). Outputs: `data/permutations/`.

---------------------------------------------------------------------------
## Measured facts

### Inputs and reference
- The pilot pipeline is `perturb_ldsc_permutation/` (tracked). Pilot runs used
  untracked copies differing only in config:
  ASTCO×Stim48 = `perturb_ldsc_permutation0.5`, ASTCO×Rest = `…0.5_rest`,
  AITD×{Rest,Stim8,Stim48} = `…0.5_{rest,stim8,stim48}_AITD`; all with
  `min_abs_logfc: 0.5`, seed 42, 500 perms. All 5 have 501 `.results`.
- `data/gene_sets/coloc_thr0.5` does not exist. `data/gene_sets/coloc0.5`
  is the **H4 > 0.5** variant (predates the logFC filter; AITD 84 anchors vs
  63). The coloc_thr0.5 sets (|logFC| ≥ 0.5, H4 > 0.8) are
  `data/sensitivity/thr0.5/gene_sets`; the pilot `real` sets are byte-identical
  to them (15/15 files). This is the reference used for fix 6. On the
  cluster, the thr0.5 reference gene sets are in
  `data/sensitivity/thr0.5/gene_sets`.
- `DE_results_strong.csv` is already strong-perturbation filtered (on-target
  KD, min DEGs, min cells, significant on-target; padj < 0.05, |z| ≥ 2,
  self-edges removed). Edges per condition: Rest 286,654; Stim8hr 447,619;
  Stim48hr 346,700 (before |logFC| ≥ 0.5).
- Every cis gene (1,564) also appears as trans, so the pilot "union" pool =
  the 10,267 trans genes. The screen measures 10,282 genes (h5ad var).
- 277 CD4 coloc genes (H4 > 0.8, any trait); 187 are measured in the screen.
  AITD: 63 anchors, 44 measured; ASTCO: 27, 21 measured. Unmeasured coloc
  genes: 35 protein-coding, 40 lncRNA, 15 other.
- `base_DESeq2` is ~constant within gene × condition (median relative SD
  0.7%). Screen baseMean per gene × condition: median over all contrasts of
  the h5ad `baseMean` layer (Rest 11,287, Stim8hr 11,415, Stim48hr 11,281
  contrasts; no NaN) → `data/permutations/inputs/screen_expression.tsv`.

### Pilot matching bug (confirmed)
- `sample_pseudo_anchors` computed L1 distance between raw (out-degree,
  mean expression) and stratum *bin indices*. Reproduced exactly: the OLD
  sampler as imported from the pilot regenerates pilot perms 1–5 for both
  traits.

### Stage 1 (AITD, ASTCO; Rest; 500 draws each; v2 variables, |logFC| ≥ 0.5)
`data/permutations/matching_check/` (summary, per-draw table, figure).

| | real | OLD (median of draw medians) | FIXED |
|---|---|---|---|
| AITD expression (measured anchors) | 29.5 | 63.7 (46.1–79.5) | 25.7 (21.2–30.1) |
| AITD in-degree median / frac zero | 7 / 0.317 | 6 / 0.101 | 7 / 0.317 |
| AITD out-degree frac zero | 0.984 | 0.972 | 0.984 |
| ASTCO expression | 27.1 | 78.9 (58.9–102.8) | 26.1 (20.2–31.9) |
| ASTCO in-degree median / frac zero | 6 / 0.222 | 5 / 0.115 | 6 / 0.222 |
| ASTCO out-degree frac zero | 0.963 | 0.891 | 0.963 |

(ranges = 2.5–97.5% of draw medians). Pooled KS vs real, FIXED: expression
0.09 / 0.07, in-degree 0.05 / 0.09, out-degree 0.02 / 0.02 (AITD / ASTCO);
OLD: 0.18 / 0.37, 0.22 / 0.13, 0.01 / 0.07.
- FIXED: 0 nearest-stratum fallbacks in 1,000 draws; 0 coloc genes in any
  pseudo set (OLD: up to 4). Measured pseudo-anchors = 44 / 21 exactly.
- Set sizes (real; FIXED median, 95% range): AITD regulators 371; 406
  (356–465). AITD targets 162; 443 (189–1,728). ASTCO regulators 217; 232
  (189–285). ASTCO targets 417; 301 (149–1,353).
- Overlap with real sets (FIXED, fraction of real set recovered): regulators
  0.66 (AITD) / 0.49 (ASTCO); targets 0.08 / 0.08.

### Container / env
- LDSC container: Python 3.9.23, polars 1.36.1, numpy 2.0.2.
  `make_permutation_gene_sets.py` + `matching.py` run inside it (verified by
  `verify_real_sets`).
- uv project: Python 3.12 (system), polars 1.31.0, numpy 2.2.6, pandas 2.3.1,
  matplotlib 3.10.3, scipy 1.16.0, h5py 3.14.0, pyyaml 6.0.2 (`uv.lock`).
  `/home` is at quota (52.4/52.4 GB, grace to 2026-11-04): `env.sh` puts uv's
  cache on RDS.
- Snakemake 9.20.0 from the existing `env_nf` conda env (with the slurm
  executor plugin); not added to uv.

### Stage 2 checks before launch
- Real gene sets for all 5 runs (built in the container) identical to
  `data/sensitivity/thr0.5/gene_sets`.
- Test chains via slurm (8 CPUs): wall 7.0–8.75 min (mean 7.8), CPU
  0.51–0.79 CPU-h (mean 0.63), max RSS ~10 GB. AITD Rest `real` `.results`
  identical to the pilot's real `.results`. Annotations / LD scores removed
  after `.results` (kept per chain ≈ 0.33 MB; intermediates ≈ 70 MB).
- Pilot sacct: 501 perms per trait×condition used 20.8–22.1 CPU-days (≈1.0–1.06
  CPU-h per chain) on one 76-core interactive node, ~7.3 h wall each.
- Disk: project quota 215,301 / 220,000 GB used (≈4.7 TB free; filesystem
  98% full, 4.3 TB available). Peak need ≈ 1.7 GB. Not tight.
- First launch (21:09, Snakemake in tmux on login-q-1) built the DAG but
  submitted no jobs (slurm-plugin hang from the login node; not debugged).
  That driver and a second login-node attempt were stopped; their empty
  slurm log dirs were removed. No chain was submitted by either.
- After the k-NN change (below) the 3 null test chains (AITD perm1/perm2 Rest,
  perm1 Stim8hr) produced different gene sets and were deleted; the AITD Rest
  `real` chain was unchanged and kept.
- **Launched 2026-10-28 ~21:45 as slurm job 39513711** (`run_driver.sbatch`:
  1 CPU, 4 GB, 20 h; icelake max is 7 days), existing profile, 20
  concurrent chains at 8 CPUs / 16 GB / 90 min each, `--rerun-triggers
  mtime`. 1,004 chains to run. Expected ≈ 630 CPU-h used (test basis;
  ≈1,045 on pilot basis / allocated) and ≈10 h wall. Driver log:
  `data/permutations/heritability/driver_39513711.log`.

### Stage 1 re-check after the out-degree k-NN change (targets only; Rest; 500 draws)
`data/permutations/matching_check/target_size_check_knn.csv`

| | knocked-down anchor (out-degree) | real targets | null median | null 95% range | null ≥ real |
|---|---|---|---|---|---|
| AITD | ENSG00000172575 (162) | 162 | 58 | 10–189 | 4.8% |
| ASTCO | ENSG00000131323 (417) | 417 | 237 | 139–658 | 19% |

Both real sizes are inside the null's 95% range. 0 k-NN relaxations and 0
fallbacks in 1,000 draws. For comparison, the quartile-bin matching gave
null medians of 443 (AITD; 189–1,728) and 301 (ASTCO; 149–1,353).
- The AITD anchor's (in-degree bin, expression bin) cell has 30 knocked-down
  pool genes; its 20 nearest by out-degree are 10–189, with 18 of the 20 below
  162. The ASTCO cell has 70; its 20 nearest are 139–658.

### Stage 3 costs (`data/permutations/plan_costs.csv`; 20 concurrent jobs)
| plan | N | jobs | CPU-h used (test / pilot basis) | wall h | peak disk GB |
|---|---|---|---|---|---|
| A: 18 traits × 3 cond | 200 | 10,854 | 6,838 / 11,288 | 107 | 4.9 |
| A | 500 | 27,054 | 17,044 / 28,136 | 266 | 10.1 |
| B: 6 IMDs + 6 controls, Stim8hr | 200 | 2,412 | 1,520 / 2,509 | 24 | 2.1 |
| B | 500 | 6,012 | 3,788 / 6,253 | 59 | 3.3 |

B = IBD, RA (EUR = RA_Ishigaki_36333501_2), ASTCO, T1D (Robertson), AITD,
SLE + controls IS, MYOI, OA, T2D (Xue), ADFH, SCZ (the repo's `NON_IMD`).
Wall assumes 7.8 min per job + 4 min scheduling/poll overhead, 20 at once;
wall scales ~inversely with concurrency (`--jobs`).

---------------------------------------------------------------------------
## Changes (fixes 1–8 and other decisions)

1. Stratum assignment: `matching.Binning` fits edges on the pool; the same
   object bins pool genes and anchors; anchors sample from their exact stratum
   without replacement. Fallback: nearest non-empty stratum by L1 bin-index
   distance (ties pooled), counted per draw (`n_fallback`, `max_fallback_dist`
   in `{study}_{cond}_matching.tsv`).
2. Out-degree: bin 0 = 0, then 4 quantile bins among non-zero pool genes.
3. In-degree added, same zero + 4-quantile scheme (22–32% of real anchors
   have in-degree 0 at Rest, so the zero bin keeps them exact).
4. Degrees per condition from the gene-set edges (|logFC| ≥ 0.5, strong DE
   table, self-edges dropped). Consequently pseudo-anchors are drawn per
   condition (seed = (42 + perm, condition code)); the pilot drew one set for
   all conditions.
5. Pool = screen-measured genes minus all 277 coloc genes (10,095 genes).
6. `min_abs_logfc: 0.5`; `verify_real_sets` must pass before any chain runs,
   and every `real` chain re-checks against the reference.
7. Collected: prop_snps, prop_h2 (+SE), enrichment (+SE, p), coefficient
   (+SE, z), set size, matching log.
9. (2026-10-28, approved) Anchors with out-degree > 0 are matched by k-NN on
   out-degree: the k = 20 knocked-down pool genes nearest in out-degree
   within the anchor's (in-degree bin, expression bin) cell; if the cell has
   < k available genes, cells within L1 bin distance 1, 2, … are added
   (`n_knn`, `n_knn_relaxed`, `max_knn_relax_dist` in the matching log).
   Ties at the k-th distance are kept.
8. `perm_chain` deletes `annotations/{cond}` and `ld_scores/{cond}` after
   `.results` exists.

Other decisions (not in the fix list; first two approved 2026-10-28):
- **Expression source**: screen baseMean from the h5ad (per condition) instead
  of the mean `base_DESeq2` over DE rows, so that anchors with no DE edges
  have a value.
- **Unmeasured anchors** (AITD 19, ASTCO 6) are matched to genes in the study's
  gene-coordinate file that are unmeasured in the screen and not coloc genes,
  split by protein-coding status (`gene_coord_protein_coding*.txt`). The
  pilot gave them expression 0.
- **Job structure**: one slurm job per perm × condition chain (1,005 jobs
  instead of ~7,000) because of the profile's 180 s status polling.
- Analysis code reads `config.yaml` for runs; collection writes
  `{study}_{cond}/` summaries, `summary_all.csv`, `pooled_mean_z.csv` and
  `matching_log_summary.csv`. (The planned `pooled_tests.csv` with a
  permutation p and BH was dropped on 2026-10-29: with two traits the pooled
  mean z is descriptive only.)

### Drug-target permutation copy (`scripts/approved_drug_enrichment_permutation_v2.py`, not run)
- Vendored `build_match_strata` / `sample_pseudo_coloc` removed; uses
  `matching.build_matcher` / `Matcher.sample` (fixes 1–3, 5).
- Per condition, expansion and degree matching use the same
  |logFC| ≥ `--min-abs-logfc` (default 0.5) edges without self-edges (fix 4);
  the original expanded through all edges, including self-edges.
- New required args: `--coloc`, `--expression`, `--gene-coord`,
  `--protein-coding-coord`; removed `--pool`, `--match-by`, `--n-bins`.
- Seeds per (perm, condition) via `perm_rng`, instead of one sequential RNG.
- Output gains `n_fallback_total`, `n_unmeasured`. Fisher universes unchanged.
- `--gene-sets-dir` must be the variant matching `--min-abs-logfc`
  (coloc_thr0.5 = `data/sensitivity/thr0.5/gene_sets`).

---------------------------------------------------------------------------
## Interpretation

- The pilot null was mis-matched mainly on **expression**: OLD pseudo-anchors
  were ~2–3× more highly expressed than the real anchors and under-sampled
  zero-in-degree genes. Out-degree was roughly right by accident (most genes
  have none). The pilot's empirical p-values therefore compared real anchors
  against a null of highly expressed, densely regulated genes.
- With the fix, the anchor-level variables match closely. At Rest, target-set
  size is driven by the single knocked-down anchor (AITD and ASTCO each have
  one), so its matching dominates the targets null. Quartile bins made that
  null lumpy and placed the real AITD size below the null's 2.5th percentile;
  k-NN on out-degree fixes that.
- **Caveat (AITD Rest targets):** within the anchor's cell the k-NN
  neighbours are mostly below the real out-degree (162 sits near the top of a
  30-gene cell), so the null target size is skewed low (median 58) and only
  4.8% of draws reach the real size. This favours "real > null" for AITD Rest
  targets independently of disease biology. Read that test with this in mind;
  if it matters, a k-NN over all knocked-down genes (dropping the cell
  restriction) or symmetric neighbours would remove the skew. ASTCO's
  neighbours straddle the real value.
- Regulator sets overlap the real ones heavily under any matched null
  (0.5–0.66 of the real set), as expected when many knockdowns hit many genes.
  The regulator test is therefore weak by construction; anchors and targets
  are the more discriminating sets.

---------------------------------------------------------------------------
## Stage 2 results (collected 2026-10-29)

Outputs in `data/permutations/heritability/`: `{study}_{cond}/`
(`permutation_ldsc_perperm.csv`, `permutation_ldsc_summary.csv`,
`permutation_null_{cond}.png`), `summary_all.csv`, `old_vs_new_comparison.csv`,
`pooled_mean_z.csv`, `matching_log_summary.csv`. Scripts:
`collect_permutation_ldsc.py`, `compare_with_pilot.py`, `estimate_costs.py`
(uv). Empirical p = (1 + #null ≥ real)/(N + 1); with N = 200 the floor is
1/201 = 0.005 (pilot, N = 500: 0.002).

### Measured facts

**Run.** Driver 39513711 ran 2026-10-28 21:28 → 2026-10-29 08:46 (11.3 h),
COMPLETED, "1005 of 1005 steps (100%) done". 1,004 chain jobs, all
COMPLETED (the AITD Rest `real` chain is the kept test chain). No
annotation or LD-score files remain.

**Completeness.** All 5 trait × condition runs: 200/200 null perms, real run
present, 0 null sets dropped by LDSC and 0 empty, for anchors, regulators
and targets (summary columns `n_null_dropped`, `n_null_empty`). Pilot, for
comparison (500 perms, from `.results`): 0 dropped for anchors and
regulators; targets dropped 38 (AITD Rest), 19 (AITD Stim8), 27 (AITD
Stim48), 7 (ASTCO Rest), 4 (ASTCO Stim48), all empty target sets.

**Pilot location.** The pilot outputs are `data/ldsc_permutation0.5*/`
(no `data/permutation_results/` exists). Each has 501 `.results`, but
`data/ldsc_permutation0.5_rest/permutation_ldsc_perperm.csv` (ASTCO Rest)
was written at 00:45, before the sweep's last `.results` at 01:03, and holds
only 465 nulls. `compare_with_pilot.py` therefore re-reads all 5 pilots from
their `.results` with the v2 collector (pilot files not modified).

**Real values.** Real enrichment / coefficient z are identical in pilot and
v2 for all 15 trait × condition × set (same real sets).

**Per set** (`summary_all.csv`; z = (real − null mean)/null SD):

| trait | cond | set | real enr | null mean (SD) | z | p | real coef z | null mean (SD) | z | p |
|---|---|---|---|---|---|---|---|---|---|---|
| AITD | Rest | anchors | 27.47 | 1.91 (1.77) | 14.5 | 0.005 | 3.70 | −0.53 (1.15) | 3.68 | 0.005 |
| AITD | Rest | regulators | 3.56 | 3.40 (0.43) | 0.37 | 0.328 | 1.81 | 1.42 (0.33) | 1.20 | 0.129 |
| AITD | Rest | targets | 3.71 | 3.33 (2.62) | 0.15 | 0.343 | 1.12 | 0.11 (1.00) | 1.01 | 0.159 |
| AITD | Stim8 | anchors | 27.45 | 2.86 (3.19) | 7.72 | 0.005 | 3.69 | −0.26 (1.05) | 3.76 | 0.005 |
| AITD | Stim8 | regulators | 4.65 | 4.15 (0.56) | 0.89 | 0.189 | 2.74 | 2.33 (0.40) | 1.03 | 0.139 |
| AITD | Stim8 | targets | 3.69 | 3.26 (0.65) | 0.66 | 0.244 | 2.24 | 1.25 (0.66) | 1.51 | 0.070 |
| AITD | Stim48 | anchors | 27.20 | 2.06 (2.81) | 8.96 | 0.005 | 3.74 | −0.64 (1.19) | 3.68 | 0.005 |
| AITD | Stim48 | regulators | 5.25 | 4.35 (0.57) | 1.56 | 0.050 | 2.88 | 2.40 (0.36) | 1.33 | 0.085 |
| AITD | Stim48 | targets | 5.28 | 3.34 (0.95) | 2.06 | 0.035 | 2.56 | 0.56 (0.60) | 3.30 | 0.005 |
| ASTCO | Rest | anchors | 75.71 | 1.71 (2.97) | 24.9 | 0.005 | 4.05 | −0.79 (1.33) | 3.64 | 0.005 |
| ASTCO | Rest | regulators | 5.16 | 3.89 (0.92) | 1.38 | 0.085 | 2.22 | 1.56 (0.45) | 1.47 | 0.075 |
| ASTCO | Rest | targets | 2.86 | 2.97 (1.57) | −0.07 | 0.433 | 0.97 | 0.16 (0.96) | 0.84 | 0.234 |
| ASTCO | Stim48 | anchors | 75.76 | 1.87 (2.62) | 28.2 | 0.005 | 4.04 | −0.54 (1.07) | 4.28 | 0.005 |
| ASTCO | Stim48 | regulators | 6.47 | 4.88 (1.03) | 1.56 | 0.060 | 1.64 | 1.95 (0.42) | −0.74 | 0.826 |
| ASTCO | Stim48 | targets | 3.85 | 3.63 (1.66) | 0.13 | 0.408 | 0.92 | 0.58 (0.89) | 0.39 | 0.468 |

AITD Stim48 regulators enrichment p = 10/201 = 0.0498.

**Set sizes** (real; null median, 95% range). Anchors: 63 / 27 in every null.
| trait | cond | regulators | targets |
|---|---|---|---|
| AITD | Rest | 371; 406 (358–467) | 162; 55 (10–189) |
| AITD | Stim8 | 491; 469 (412–527) | 2,096; 658 (349–3,586) |
| AITD | Stim48 | 396; 396 (327–470) | 416; 195 (105–454) |
| ASTCO | Rest | 217; 233 (188–287) | 417; 237 (139–658) |
| ASTCO | Stim48 | 291; 273 (218–341) | 604; 298 (142–657) |

Pilot null target sizes had 2.5th percentile 0 (AITD) / 21 (ASTCO) and
97.5th 1,258–1,787.
- Distinct null target sets out of 200: AITD Rest 20, ASTCO Rest 20 (one
  knocked-down anchor, k = 20 neighbours); AITD Stim8 197, Stim48 191, ASTCO
  Stim48 200.
- Ad hoc check (not saved to a file): Spearman ρ between null target size and
  the statistic is between −0.30 and 0.35 across the 10 target tests. For AITD
  Stim48 targets, among the 19 nulls with size ≥ 0.75 × real, none reaches the
  real enrichment or coefficient z (p = 1/20 = 0.05 for both).

**Old vs new** (`old_vs_new_comparison.csv`; old = pilot, N = 500 (targets
462–496 after drops)).
- Empirical p crosses 0.05 in 4 of 30 tests:
  - AITD Rest regulators coef z: 0.0499 → 0.129
  - ASTCO Rest regulators coef z: 0.018 → 0.075
  - AITD Stim48 regulators enrichment: 0.120 → 0.0498
  - AITD Stim48 targets enrichment: 0.118 → 0.035
  On the first 200 pilot perms the old p values are 0.050, 0.030, 0.114, 0.116.
- AITD Stim48 targets coef z: 0.021 → 0.005 (< 0.05 in both).
- Anchors: p at the floor in both (0.002 / 0.005), all 10 tests. Null mean
  enrichment fell in all 5 (AITD 2.93→1.91, 2.93→2.86, 2.94→2.06; ASTCO
  2.60→1.71, 2.61→1.87); null mean coef z fell in 4/5 (by 0.07–0.34) and rose
  0.10 in AITD Stim8.
- Regulators: null mean enrichment changed by −0.08 to +0.13; null mean
  coef z rose in 4/5 (+0.02 to +0.22), fell 0.05 in AITD Stim8.
- Targets: null mean enrichment −0.63 to +0.39; coef z −0.46 to +0.45.

**Pooled mean z** (`pooled_mean_z.csv`; Stim8hr has AITD only).
| set | stat | Rest (AITD; ASTCO) | Stim8hr (AITD) | Stim48hr (AITD; ASTCO) |
|---|---|---|---|---|
| anchors | enrichment | 19.7 (14.5; 24.9) | 7.7 | 18.6 (9.0; 28.2) |
| anchors | coef z | 3.66 (3.68; 3.64) | 3.76 | 3.98 (3.68; 4.28) |
| regulators | enrichment | 0.87 (0.37; 1.38) | 0.89 | 1.56 (1.56; 1.56) |
| regulators | coef z | 1.33 (1.20; 1.47) | 1.03 | 0.30 (1.33; −0.74) |
| targets | enrichment | 0.04 (0.15; −0.07) | 0.66 | 1.09 (2.06; 0.13) |
| targets | coef z | 0.93 (1.01; 0.84) | 1.51 | 1.85 (3.30; 0.39) |

**Matching log** (`matching_log_summary.csv`, 200 null perms each).
| trait | cond | fallbacks | k-NN anchors | perms with k-NN relaxation (total; max dist) | unmeasured anchors |
|---|---|---|---|---|---|
| AITD | Rest | 0 | 1 | 0 | 19 |
| AITD | Stim8 | 0 | 3 | 0 | 19 |
| AITD | Stim48 | 0 | 3 | 0 | 19 |
| ASTCO | Rest | 0 | 1 | 0 | 6 |
| ASTCO | Stim48 | 0 | 5 | 200 (200; dist 1) | 6 |

**Cost (sacct, 1,004 chain jobs, 8 CPUs each).** CPU used (TotalCPU)
581.2 CPU-h (per chain mean 0.579, median 0.564, range 0.46–0.97);
allocated 1,155.9 CPU-h. Chain wall mean 8.63 min (median 8.2, 5.9–21.3).
Queue wait median 0.3 min (max 5.0). Max RSS 10.2 GB. Chains ran 21:28 →
08:20 (10.9 h); effective cycle 13.0 min per job slot at 20 concurrent
(≈ 4.3 min overhead per job). The driver used 10.3 CPU-h over 11.3 h. The
driver log shows the last chain as finished at 08:43, 24 min after the
last sacct end time (08:20); not investigated. Kept outputs 310 MB
(0.31 MB per chain).

**Recomputed plan costs** (`data/permutations/plan_costs.csv`; basis 0.579
CPU-h and 8.63 + 4.34 min per chain, 20 concurrent; this replaces the
Stage 3 table above, and `cpu_hours_used_test_basis` is now
`cpu_hours_used`). The model gives 11.0 h for Stage 2 (actual 11.3 h).
| plan | N | jobs | CPU-h used | CPU-h allocated | wall h | peak disk GB |
|---|---|---|---|---|---|---|
| A: 18 traits × 3 cond | 200 | 10,854 | 6,285 | 12,489 | 117 | 4.7 |
| A | 500 | 27,054 | 15,664 | 31,130 | 293 | 9.6 |
| B: 6 IMDs + 6 controls, Stim8hr | 200 | 2,412 | 1,397 | 2,775 | 26 | 2.1 |
| B | 500 | 6,012 | 3,481 | 6,918 | 65 | 3.2 |

### Interpretation

- **Anchors** remain far above the corrected null in every trait ×
  condition, for both statistics (real exceeds all 200 nulls). The fix
  lowered the null mean enrichment (the old null's pseudo-anchors were
  over-expressed), so the anchor result was conservative in the pilot and
  strengthens slightly. This is the S-LDSC signal of the coloc genes
  themselves and is expected given how they were chosen.
- **Regulators**: no set clears the null robustly. The two pilot "hits"
  (AITD Rest and ASTCO Rest coefficient z) lose significance. That is partly
  the null coefficient z rising and partly p-value resolution / Monte Carlo
  noise: AITD Rest was already borderline at 0.0499, and the first 200 pilot
  perms give 0.050 and 0.030. The new AITD Stim48 enrichment p = 0.0498 sits
  on the threshold. Real regulator enrichments are consistently a little above
  the null mean (z 0.4–1.6), as expected when 50–66% of the real set is in
  every null; a modest pooled shift, not a regulator-specific signal.
- **Targets**: the only non-anchor result that looks real is AITD Stim48hr
  (coefficient z p = 0.005, enrichment p = 0.035; pilot coef z p = 0.021).
  It does not appear to be driven by set size: ρ ≈ 0 and it holds among
  size-comparable nulls, though only 19 of those exist, so p ≥ 0.05 there.
  It is one trait × condition of ten target tests with no multiplicity
  correction, and ASTCO Stim48 targets are null (z 0.13 / 0.39). Treat it as
  a lead, not a finding.
- **Rest targets nulls have only 20 distinct sets** (one knocked-down anchor,
  k = 20), so their effective N is 20 and the p-values are coarse and lumpy
  (visible in the AITD Rest figure). The caveat above about the AITD Rest
  null skewing low does not produce a signal: AITD Rest targets p = 0.34 /
  0.16.
- **Matching** behaved as designed: no stratum fallbacks; the k-NN
  relaxation fired only for one ASTCO Stim48 anchor (always, at distance 1),
  i.e. its cell has < 20 knocked-down genes. Unmeasured-anchor counts match
  Stage 1 (19 / 6).
- **Pooled means** are descriptive (two traits, one at Stim8). They repeat
  the per-trait picture: anchors z 3.7–4.0 on coefficient z, everything
  else ≤ 1.9, and the highest non-anchor value (targets Stim48, 1.85) comes
  from AITD alone (ASTCO 0.39).
- **Costs**: the measured per-chain CPU (0.58 CPU-h) is below the test basis
  (0.63) and ~half the pilot's; wall per chain is a bit longer (8.6 vs 7.8
  min), so plan A at N = 200 is ≈ 117 h wall at 20 concurrent jobs.
  Allocated CPU-h is ~2× used (8 CPUs per chain, ~50% used); if the account
  is billed on allocation, 4 CPUs per chain would roughly halve the charge at
  some cost in wall time (not tested).
