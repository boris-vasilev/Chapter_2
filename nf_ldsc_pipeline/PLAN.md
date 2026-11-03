# Nextflow rewrite plan: `nf_ldsc_pipeline` and `nf_ldsc_permutation`

Planning only. Nothing has been built or run. `keep` = port into the new pipeline, `drop` = do not
port (the file stays where it is, nothing is deleted), `not sure` = your call (section 10).

## Findings that change the brief

1. **"One process per step" does not scale.** Nextflow has no equivalent of Snakemake's job groups,
   so each process task is one SLURM job: 7 jobs per chain, ~76,000 for 10,854 chains. It also
   conflicts with per-permutation cleanup, because deleting another task's outputs breaks `-resume`.
   Decision: annotation → combine → LD scores → h2 is one process, `LDSC_CHAIN`, used by both
   pipelines (section 4; question 1).
2. **Dropping the 75-DEG filter changes matching more than gene sets.** Upper-bound proxy from
   `data/DE_results_sig.csv` (on-target KD, padj < 0.05 only), at |logFC| ≥ 0.5: knockdowns with an
   edge go from 894–970 to 6,797–7,111 per condition, edges rise only 14–21%. Most knockdowns have
   few edges, so k-NN matching will apply to many more anchors than the current 1–5 per run.
3. **`extract_de_results.py` will probably run out of memory unfiltered.** It builds the full
   obs × gene table before filtering: ~29M rows now, ~217M without the DEG filter (from the proxy),
   and `obs_idx` is a Python list of that length (line 114). The 16 GB allocation will very likely
   fail. It needs a fix in either pipeline.
4. **"Exactly" must mean set equality for main-pipeline gene sets.** `coloc_thr0.5` and
   `sensitivity/thr0.5` hold identical sets in all 162 files, but only 17 are byte-identical:
   `make_gene_sets_coloc.py` writes polars `.unique()` order. Permutation gene sets are written
   sorted, so those can be compared byte for byte.
5. **Your sacct premise is unconfirmed.** The Stage 2 Snakemake log has no `sacct` lines and reports
   squeue as the status source. The 180 s `seconds-between-status-checks` in the profile is a
   sufficient explanation for much of the 4.3 min per-job overhead. Nextflow avoids this either way:
   it detects completion from `.exitcode` files, not from scheduler accounting.
6. **`collect_permutation_ldsc.py` will mislabel single-annotation runs in plan A.** It matches
   `{annot}L2`, but when one set survives LDSC names it `L2_1`. That affects Osteoarthritis (all
   conditions) and the Alzheimer's proxy at 48 h, which would be recorded as dropped. Stage 2 is
   unaffected. `chapter2/vis` handles this with a hard-coded table.
7. **`coloc_strict` does not match its config.** It has 84 AITD anchors, which is the H4 > 0.5 count
   (63 at H4 > 0.8), while `config.yaml` says 0.8.

---

# Part 1: main pipeline

## 1. Inventory

The Snakemake copy in `chapter2/perturb_ldsc_pipeline` is identical to the root
`perturb_ldsc_pipeline`. The sensitivity, pilot and v2 Snakefiles use the **root** copy
(`PIPE_SCRIPTS`), so neither is affected by this work. `chapter2/vis` calls no pipeline script: it
reads gene sets, `.results`, `.log` and `DE_results_strong.csv`.

### scripts/

| File | Description | Status | Reason | Used by |
|---|---|---|---|---|
| extract_de_results.py | h5ad → filtered long DE table | keep (edit) | Filters become arguments; memory fix | rule extract_de_results; output read by chapter2/vis, v2, drug vis |
| make_gene_sets_coloc.py | Regulator/target/anchor sets per study × condition | keep | Works; already has `--min-abs-logfc` | rule make_gene_sets, sensitivity Snakefile, generate_threshold_gene_sets.py |
| run_make_annot.py | make_annot over chr 1–22 | keep | Works | annots.smk, sensitivity, pilot, v2 |
| combine_annots.py | Joins annotations, drops sets < min SNPs | keep | Works | annots.smk, sensitivity, pilot, v2 |
| run_ldsc_l2.py | ldsc --l2 over chr 1–22 | keep | Works | ld_scores.smk, filter_annots.smk, sensitivity, pilot, v2 |
| run_munge_sumstats.py | Format-detecting munge wrapper | keep | Needed for new traits | sumstats.smk |
| filter_annots.py | Second, LD-score-based column filter | drop | Off (`filter_low_snp_annots: false`); `combine_annots` already filters | filter_annots.smk only |
| make_gene_sets_programs.py | Zhu programme sets | drop | Confirmed: rule only | gene_sets_programs.smk |
| collect_program_results.py | Programme results table | drop | Confirmed: its output is read only by the two programme vis scripts | gene_sets_programs.smk |
| make_gene_sets_l2g.py | L2G sets | drop | Confirmed: rule only | gene_sets_l2g.smk |
| build_l2g_study_mapping.py | L2G study mapping | drop | Confirmed: no references anywhere | none |
| select_imds_from_opentargets.py | L2G bootstrap step 1 | drop | L2G only; not in your list | README only |
| compile_studies_from_opentargets.py | L2G bootstrap step 2 | drop | L2G only; not in your list | README, make_gene_sets_l2g.py docstring |
| download_opentargets_sumstats.py | L2G bootstrap step 3 | drop | L2G only; not in your list | README only |
| make_gene_sets_tenk10k.py | Gene sets from the TenK10K coloc tree | not sure | Not in the Snakemake DAG and not `coloc_dat`, but it writes the `disease_groups.csv` taken by two vis scripts | Manual; approved_drug_enrichment_all.py, program_enrichment.py |

The programme and L2G drops are safe. `gene_set_analysis/6_regulator_program_structure.py` and
`perturb_programme_heritability/` read `data/programs_zhu/` directly, not these outputs. Leave
`data/opentargets/` in place: the drug scripts and sensitivity config need
`evidence_clinical_precedence` and `disease`.

### rules/

| File | Description | Status | Reason | Used by |
|---|---|---|---|---|
| common.smk | Paths, study list, reference selection | keep (as `refs.nf`) | Only `_ref_key` logic and lookups needed | All rules; copied into sensitivity, pilot, v2 |
| extract_de.smk | Runs extraction | keep | Now passes all filters | Snakefile |
| gene_sets_coloc.smk | Runs coloc gene sets | keep | H4 only | Snakefile |
| sumstats.smk | Runs munge | keep | | Snakefile |
| annots.smk | make_annot + combine | keep (merged) | Into `LDSC_CHAIN` | Snakefile |
| ld_scores.smk | ldsc --l2 | keep (merged) | Into `LDSC_CHAIN` | Snakefile |
| heritability.smk | ldsc --h2 | keep (merged) | Into `LDSC_CHAIN` | Snakefile |
| filter_annots.smk | Optional second filter | drop | Off | Snakefile (conditional) |
| gene_sets_programs.smk | Programme backend | drop | Confirmed | Snakefile |
| gene_sets_l2g.smk | L2G gene sets | drop | Confirmed | Snakefile |

### vis/

No rule calls these, so all are `drop` (not ported, left in place). All except where noted take
paths as arguments, so they work with the new outputs.

| File | Description | Status | Superseded by chapter2/vis? | Used by |
|---|---|---|---|---|
| parse_ldsc_results.py | `.results` → tidy CSV | drop | Yes: heritability/01_collect_results.py | Three plot scripts below |
| plot_ldsc_results.py | z heatmap + forest | drop | Yes: Fig. 4, S5 | Reads parse output |
| plot_enrichment_panels.py | Four-panel heatmap | drop | Yes: Fig. 4 | Reads parse output |
| plot_enrichment_pvalue_heatmap.py | Enrichment/p heatmap | drop | Yes: Fig. 4 | Reads parse output |
| enrichment_vs_geneset_size.py | Enrichment vs set size | drop | Yes: Fig. S6 | Reads parse output |
| gene_count_barplot.py | Set-size bars | drop | Yes: 02_gene_set_sizes_barplot.py, Fig. S4 | Gene sets |
| study_overlap.py | Pairwise overlap heatmaps | drop | Yes: 04_jaccard_heatmaps.py | Gene sets |
| permutation_test_overlap.py | Overlap permutation test | drop | Partly: 04 uses an analytic expectation | Gene sets; defaults to old `DE_results_sig.csv` |
| gene_vs_snp_overlap.py | Gene vs SNP overlap | drop | No | Gene sets |
| regulator_anchor_matrix.py | Anchor × regulator heatmap | drop | No | Gene sets, DE |
| regulator_anchor_network.py | Anchor–regulator network | drop | No | Gene sets, DE |
| approved_drug_enrichment.py | Drug-target Fisher test | drop | No | Gene sets, DE |
| approved_drug_enrichment_all.py | Per-disease drug enrichment | drop | No | Gene sets, DE; constants shared by name with others |
| approved_drug_enrichment_drilldown.py | Genes behind a hit | drop | No | Gene sets, DE |
| approved_drug_enrichment_network.py | Drug-hit network | drop | No | Gene sets, DE |
| approved_drug_enrichment_permutation.py | Drug-target permutation (old matching) | drop | No; v2 copy in null_v2 | Pilot and sensitivity vis |
| permutation_heatmap.py | Permutation overlay heatmap | drop | No | `permutation_results.csv` |
| plot_drug_enrichment_null.py | Drug null histograms | drop | No | `permutation_results.csv` |
| disease_specificity_matrix.py | Set × indication test | drop | No | Gene sets, DE |
| edge_specificity_exploratory.py | Signed-edge exploration | drop | No | Gene sets, DE |
| novel_candidates.py | Candidate ranking | drop | No | Gene sets, DE |
| program_enrichment.py | Programme over-representation | drop | No (programmes out of scope) | Gene sets, `programs_zhu` |
| ldsc_program_enrichment.py | Programme S-LDSC plots | drop | No (programmes out of scope) | `program_enrichment.csv` |
| program_heritability_heatmap.py | Programme heatmap | drop | No (programmes out of scope) | `program_enrichment.csv` |

Old outputs stay on disk, so nothing reading them breaks. Only `chapter2/vis` needs a change, and
only for new variants (section 4).

## 2. Config

| Key | Controls | Needed? |
|---|---|---|
| workspace | Project root | Yes (`data_root`) |
| container | LDSC image | Yes |
| de_input_h5ad | DE input | Yes |
| de_results | DE output; name says "strong" | Replaced by per-variant path |
| gene_sets_dir, annotations_dir, ld_scores_dir, ldsc_results_dir | Output roots | Replaced by `outdir` |
| ld_scores_filtered_dir | Second-filter output | No |
| padj_threshold | Edge padj | Yes |
| gene_set_module | Selects module | No |
| coloc.study_info, gwas_dir, sumstats_dir, coloc_dat, cell_type, h4_threshold | Coloc inputs | Yes |
| coloc.h3_threshold, h4_only | H3 sets | No. `--h4-only` is `store_true` with `default=True`, so H3 cannot be enabled anyway |
| coloc_strong | Unused block | No. Snakefile line 11 rejects this module |
| coloc_strict (min_reg_k, min_tgt_k) | Convergence filter | Keep as `min_reg_k`, `min_tgt_k` (default 1) |
| programs, l2g | Dropped modules | No |
| conditions | Culture conditions | Yes |
| window_size | Gene ± bp | Yes |
| ldsc_dir | LDSC in container | Yes |
| gene_coord, bim, snp_list, baseline_ld, weights, allele_frq | Reference data | Yes |
| filter_low_snp_annots | Second filter switch | No |
| min_snps_per_annot | `combine_annots` threshold | Yes |

**Filters in `extract_de_results.py`.** None of the first four are in `config.yaml`; `extract_de.smk`
lines 11–14 pass only `--padj-threshold`, so they are silent argparse defaults.

| Filter | Applied | Default |
|---|---|---|
| on-target KD category | line 58 | always on |
| ≥ 75 DEGs (`n_total_de_genes`) | lines 62–64 | line 35 |
| ≥ 50 cells (`n_cells_target`) | lines 65–67 | line 43 |
| significant on-target effect | lines 68–70 | always on |
| padj < 0.05 | line 136 | config |
| \|z\| ≥ 2 | line 137 | line 50; has no effect on the strong table (min \|z\| = 2.26) |
| no self-edges | line 138 | always on |

`on_target` will control both always-on knockdown filters (lines 58 and 68–70).

**What assumes the filter:**

- **Filename:** `data/DE_results_strong.csv` is fixed in `config.yaml:15`, null_v2 `config.yaml:5`,
  `chapter2/vis/common.py:241` and `05_perturbation_graph.py:301`.
- **Reference check:** `verify_real_sets` compares with `data/sensitivity/thr0.5/gene_sets`, valid
  only for default filters; `--reference-dir` must become optional.
- **k-NN matching (`matching.py`):** no code change, but its behaviour shifts (finding 2).
  `n_knn` will rise and the Stage 2 k-NN checks must be repeated for any relaxed null.
- **Chapter text:** Fig. 2.1 ("few are strong knockdowns"), the Fig. 2.2 n values (894, 970, 934)
  and the `INDEX.md` variant line describe strong perturbations only.
- **Not affected:** `build_screen_expression.py` reads all contrasts.

## 3. Current coloc DAG

Runner: `slurm/run.sbatch` (22 CPUs, one allocation per study) or `slurm/config.yaml`. Resources
below are from `slurm/config.yaml`. All steps use `containers/ldsc39_env.sif`.

| # | Step | Inputs | Outputs | Script | CPU / mem / time |
|---|---|---|---|---|---|
| 1 | extract_de_results | `GWCD4i.DE_stats.h5ad` | `DE_results_strong.csv` | extract_de_results.py | 4 / 16 GB / 60 min |
| 2 | make_gene_sets | DE table, `coloc_dat.tsv` | `{study}_{cond}_{annot}.txt` (3 × 54) | make_gene_sets_coloc.py | 1 / 8 GB / 30 min |
| 3 | munge_sumstats | `GWAS/{study}.tsv.gz`, study info | `sumstats/{study}.sumstats.gz` | run_munge_sumstats.py | 1 / 8 GB / 30 min |
| 4 | make_annot (×3) | Gene set, gene coords, bim | `{cond}_{annot}.{chr}.annot.gz` | run_make_annot.py | 8 / 8 GB / 30 min |
| 5 | combine_annots | Step 4 | `{cond}_combined.{chr}.annot` | combine_annots.py | 1 / 4 GB / 15 min |
| 6 | calculate_ld_scores | Step 5, plink, HM3 list | `.l2.ldscore.gz`, `.M`, `.M_5_50` | run_ldsc_l2.py | 8 / 16 GB / 60 min |
| 7 | calculate_heritability | Steps 3 and 6, baselineLD v2.2, weights, frq | `{cond}_combined.{results,log,cov}` | ldsc.py --h2 | 1 / 8 GB / 30 min |

`--max-parallel` defaults to 22 on 8 allocated CPUs in steps 4 and 6; the v2 Snakefile corrects this.

Current outcome for thr0.5: **49 of 54 runs have `.results`.** SCZ × 3 fail with "0 SNPs remain"
(log only; the sumstats do not merge with the reference). ADFH Rest and 8 h produce no files because
`combine_annots` refuses when every set is under 1,000 SNPs. These must stay distinguishable:
`01_collect_results.py` reads log-only as "run failed" and no files as "not run".

| Variant | Location | How produced |
|---|---|---|
| coloc | `data/gene_sets/coloc`, `data/ldsc_results/coloc` | Main pipeline, H4 > 0.8, no logFC filter |
| coloc_strict | `data/{gene_sets,ldsc_results}/coloc_strict` | Main pipeline, k = 3; H4 > 0.5 on disk (finding 7), 15 result folders |
| coloc_strong | none | Config block only; cannot run |
| coloc0.5 | Cluster only, per NOTES | H4 > 0.5; no config |
| sensitivity thr{X} | `data/sensitivity/thr{X}/` | `generate_threshold_gene_sets.py` (10 thresholds) for sets; `perturb_ldsc_sensitivity/Snakefile` (coarse grid: 0, 0.25, 0.5, 1, 1.5, 2) for S-LDSC |
| coloc_thr0.5, coloc_thr1.0 | `data/ldsc_results/` (local) | Copies of `sensitivity/thr0.5` and `thr1` (their logs show those `--out` paths), not separate runs. `coloc_thr0.5` is the chapter variant |

The sensitivity README command omits `--h4-threshold`, and that script's default is 0.5. All ten
local thr folders do have 63 AITD anchors, so they are correct.

## 4. Proposed design

```
nf_ldsc_pipeline/
  main.nf  nextflow.config  run_driver.sbatch  README.md
  modules/prepare.nf     EXTRACT_DE, MAKE_GENE_SETS, MUNGE_SUMSTATS
  modules/ldsc_chain.nf  LDSC_CHAIN
  modules/refs.nf        refs(study): reference paths for a study
  bin/                   the 6 kept scripts
```

| Process | Does | Parallel over | Resources |
|---|---|---|---|
| EXTRACT_DE | Step 1 | Once per filter set | 4 / 64 GB / 2 h (set after fix) |
| MAKE_GENE_SETS | Step 2; emits a directory | Once | 1 / 8 GB / 30 min |
| MUNGE_SUMSTATS | Step 3; `storeDir data/sumstats`, so existing files are reused | Study | 1 / 8 GB / 30 min |
| LDSC_CHAIN | Steps 4–7 | Study × condition (× perm) | `params.chain_cpus` (8) / 16 GB / 90 min |

**Chromosomes stay inside the chain**, using the existing Popen scripts with
`--max-parallel ${task.cpus}`. Per-chromosome tasks would multiply job counts by 22.

`LDSC_CHAIN` takes `tuple val(meta), path(gene_sets_dir), path(sumstats)`, where `meta` holds study,
condition, perm and publish path. It knows nothing about real versus null. Its script:

```bash
# make_annot x3 -> combine -> --l2 -> --h2, then remove intermediates unless keep_ld
for a in regulators_H4 targets_H4 H4_genes; do run_make_annot.py ...; done
# exit 3 = no set reached min_snps: expected, so no files are published ("not run")
set -o pipefail   # otherwise tee hides a combine_annots failure
combine_annots.py ... | tee combine.txt || exit 3
run_ldsc_l2.py ...
# A failed h2 keeps its .log so chapter2/vis can report "run failed" (SCZ)
ldsc.py --h2 ... --overlap-annot --print-coefficients --print-cov --out ... || true
[ "$keep_ld" = true ] || rm -rf annot ld
```

- **Outputs:** `.log` and `combine.txt` required; `.results` and `.cov` optional; LD scores only when
  `keep_ld` is set.
- **`errorStrategy`:** exit 3 ignore; 137/140/143 retry (max 2); otherwise `finish`.
- **Environment, in the config `env` scope:** `OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=1` and
  `LC_ALL=C.UTF-8`. Your CSD3 shell exports 22 threads, which caused the earlier thrashing. The
  driver also needs UTF-8 so the JVM can handle `SCZ_Pardiñas`.

**Parameters:**

| Parameter | Default | Note |
|---|---|---|
| min_de_genes | 75 | `0` turns it off |
| min_cells | 50 | |
| on_target | true | New flag on the extraction script |
| padj, min_abs_z | 0.05, 2.0 | Previously silent |
| min_abs_logfc | 0.5 | Chapter variant |
| h4_threshold, cell_type | 0.8, CD4 | Always passed: the script's default is 0.5 |
| min_reg_k, min_tgt_k | 1, 1 | |
| window_size | 10000 | |
| min_snps_per_annot | 1000 | |
| studies, conditions | study info, all 3 | Subsets for tests |
| chain_cpus, keep_ld | 8, true | `keep_ld` false for nulls |
| outdir | `data/nf` | |

**Profiles:** `local` (no container, for dry runs on small data) and `slurm` (section 7). Apptainer
uses the two `-B` binds from `run.sbatch`, plus `--fakeroot`.

**Reference paths stay absolute strings, not staged inputs.** There are 22 × several file sets per
chain, and `01_collect_results.py` parses `--frqfile-chr \S*/data/...` from the `.log`, which staged
paths would break.

**Variant naming:** the variant tag is derived from the parameters, so two runs cannot share a
folder:

```
deg75_cells50_ot1_lfc0.5_h40.8_k1-1_w10000_snp1000
```

```
data/nf/<tag>/params.json
data/nf/<tag>/DE_results.csv
data/nf/<tag>/gene_sets/{study}_{cond}_{annot}.txt
data/nf/<tag>/ldsc_results/{study}/{cond}_combined.{results,log,cov}
data/nf/<tag>/ld_scores/{study}/...                 # keep_ld only
```

The layout below the tag matches what `chapter2/vis` reads. To use a new variant there,
`common.variant_config` (regex and two paths) needs a small edit, after this work.

---

# Part 2: permutation wrapper

## 5. Inventory: `perturb_heritability_null_v2`

| File | Description | Status | Reason | Used by |
|---|---|---|---|---|
| **Pipeline code** | | | | |
| scripts/matching.py | Matched sampler, `perm_rng` | keep (byte-identical) | Required | make_permutation_gene_sets, checks, drug v2 |
| scripts/make_permutation_gene_sets.py | One perm × condition's sets | keep (~3 lines) | `--reference-dir` optional | Snakefile |
| scripts/collect_permutation_ldsc.py | Per-perm tables, summaries | keep (~16 lines) | `--out-root`, `--runs` replace the YAML read; `L2_1` fix from `combine.txt` | README; compare_with_pilot imports it |
| scripts/plot_permutation_null.py | Null histograms | keep | Imported by collect | collect |
| scripts/build_screen_expression.py | Per-condition baseMean | keep | Runs once; needs `h5py` (via scanpy in the container; to confirm) | Input to matching |
| Snakefile | Driver | drop (replaced) | `main.nf` | run_driver.sbatch |
| config.yaml | Run settings | drop (replaced) | `nextflow.config`, but still read by the uv scripts below | Snakefile, uv scripts |
| run_driver.sbatch | Driver job | drop (replaced) | New driver | |
| **One-off analysis and figures: stay as the uv project** | | | | |
| scripts/matching_check.py | Stage 1 check | drop | Analysis | matching_check_knn imports it |
| scripts/matching_check_knn.py | Stage 1 re-check | drop | Analysis; source data for fig_s_matching_check | |
| scripts/target_size_check.py | Target-size check | drop | Analysis | |
| scripts/compare_with_pilot.py | Old vs new table | drop | Analysis | |
| scripts/estimate_costs.py | Plan costs | drop | Analysis; update its overhead input after the first run | |
| figures/fig_common.py | Figure loaders | drop | `HERIT` is hard-coded to `data/permutations/heritability`; needs an argument for new variants | Figure scripts |
| figures/fig_5_matched_null.py | Fig. 5 | drop | Figure | |
| figures/fig_s_null_histograms.py | Supp. histograms | drop | Figure | |
| figures/fig_s_matching_check.py | Supp. matching figure | drop | Figure | |
| figures/table_s_matched_null.py | Supp. table | drop | Figure | |
| env.sh, pyproject.toml, uv.lock | uv environment | drop | Stay with the uv project | |
| README.md, NOTES_null_v2.md | Documentation | drop | Stay as the record of Stage 1–2 | |
| **Undecided** | | | | |
| scripts/approved_drug_enrichment_permutation_v2.py | Drug permutation with v2 matching, never run | not sure | Imports `matching`; analysis or pipeline? | None yet |

Inputs not built by any step: `gene_coord_protein_coding*.txt` come from
`4_create_gene_coord_file.py --protein-coding-only`, a one-off. They stay as fixed inputs.

## 6. Wrapper design

```
nf_ldsc_permutation/
  main.nf  nextflow.config  run_driver.sbatch  README.md
  modules/permutation.nf   BUILD_EXPRESSION, PERM_GENE_SETS, COLLECT
  bin/                     the 5 kept scripts
```

**Reuse:** `include { LDSC_CHAIN } from '../nf_ldsc_pipeline/modules/ldsc_chain'`, the same for
`EXTRACT_DE` and `refs`, and `includeConfig '../nf_ldsc_pipeline/nextflow.config'`. A project's
`bin/` may not be on `PATH` for a module run from another project, and I could not confirm this in
the docs. The shared module therefore calls scripts through
`python ${moduleDir}/../bin/run_make_annot.py`, which works from either entry point.

**One code path:** both pipelines end in `LDSC_CHAIN`. The only difference is where the gene-set
directory comes from.

```
EXTRACT_DE ─────┐
                ├─ PERM_GENE_SETS (study × cond) ── .collect() ── flatMap to one tuple per perm
BUILD_EXPRESSION┘                                                        │
                                                                    LDSC_CHAIN ── COLLECT
```

- **`PERM_GENE_SETS`:** one job per study × condition. A bash loop calls the existing
  `make_permutation_gene_sets.py` for `real` and `1..n_perm`, so matching and seeds are unchanged:
  `perm_rng(seed, perm, cond)`, i.e. `default_rng([seed + perm, cond_code])`. This avoids one job per
  perm for a step that needs only seconds (not timed separately; roughly 15–35 min per job at
  N = 200).
- **Real-set gate:** `.collect()` holds every chain until all gene-set jobs have finished, and
  `errorStrategy 'terminate'` on this process preserves `verify_real_sets`: a mismatch stops
  everything before LDSC runs. The check applies only when `params.reference_gene_sets` is set, so
  relaxed runs pass none.
- **`EXTRACT_DE` is skippable:** `params.de_results` uses an existing table. The default for
  validation is `data/DE_results_strong.csv`. Otherwise filters are inherited from the main config,
  so edges and matching use the same table.
- **`COLLECT`:** runs the collector once in the container after all chains. It adds `--plot` output
  as today and works on partial sweeps. If the container lacks something, fall back to the current
  `uv run` step.

**Parameters:** `runs` (map of study → conditions, as in `config.yaml`), `n_perm` (200), `seed` (42),
`n_degree_bins` (4), `n_expr_bins` (5), `knn_k` (20), `reference_gene_sets`, plus all edge filters
inherited from the main pipeline.

**Cleanup:** `rm -rf annot ld` as the last line of `LDSC_CHAIN`, plus `scratch true` on slurm. All
~70 MB of intermediates stay on node-local disk and never reach RDS. Because it is one task,
`-resume` keeps working, which is why cleanup is simple here and awkward with step processes.

**Outputs:** `data/permutations/<tag>/{study}/perm{P}/{gene_sets,ldsc_results}/`, the same shape as
Stage 2, so the collector and figures need no restructuring. Summaries go to `<tag>/{study}_{cond}/`
and `<tag>/summary_all.csv`. Stage 2 in `data/permutations/heritability/` is never written to.

## 7. Scheduling

**Job granularity: one SLURM job per perm × condition chain.** Reasons:

- Per-step would mean ~76,000 jobs at plan A.
- Batching multiple chains into one task loses per-chain retry and `-resume`, and it is not needed:
  submission is under 1 job per second even at peak.
- Stage 2 chains (5.9–21.3 min, median queue wait 0.3 min) suit this size.

Nextflow job arrays are not proposed. They add restrictions I could not confirm from the docs, and
the sbatch load does not need them.

| `executor` setting | Value | Why |
|---|---|---|
| queueSize | 40 | Double Stage 2; confirm the account limit (320 CPUs) |
| pollInterval | 30 sec | Reads `.exitcode` from disk, so no scheduler load |
| queueStatInterval | 3 min | One `squeue`; the CSD3 rule is ≥ 2 min between calls |
| submitRateLimit | `60/1min` | Matches `max-jobs-per-second: 1` |
| exitReadTimeout | 5 min | Tolerates RDS file-attribute lag (default 270 s) |

The driver stays a 1-CPU SLURM job as in `run_driver.sbatch` (running from the login node hung in
Stage 2): 8 GB, `NXF_OPTS='-Xms1g -Xmx6g'`, up to 7 days. It needs `NXF_HOME`, `NXF_WORK` and the
Apptainer cache on RDS, because `/home` is full. Resubmit with `-resume` to continue.

Estimates, assuming 9.2 min per slot (8.63 min chain + ~0.5 min pickup). Queue wait at 40 concurrent
jobs is unmeasured.

| Run | Chains | CPU-h allocated | Snakemake at 20 | Nextflow at 20 | Nextflow at 40 |
|---|---|---|---|---|---|
| Stage 2 repeat | 1,005 | 1,156 | 10.9 h (measured) | 7.7 h | 3.9 h |
| Plan B, N = 200 | 2,412 | 2,775 | 26 h | 18.5 h | 9.2 h |
| Plan A, N = 200 | 10,854 | 12,489 | 117 h | 83 h | 42 h |

Concurrency, not the cheaper polling, is the larger gain, and it costs no extra CPU-hours. Allocated
CPU-h are about twice CPU used. Test `chain_cpus = 4` on 20 chains: adopt it if wall time is under
about 1.7× the 8-CPU time, as billing would then roughly halve.

**Fallback if the queue or admins object:** `-profile node` runs the same pipeline with the local
executor inside one 76-core allocation, 19 chains × 4 CPUs (the layout already tested for the
pilot). About 6 lines of config and no code change.

---

# Part 3: checks and decisions

## 8. Validation

All on the cluster: the h5ad, sumstats and reference panels are not local, and neither is Stage 2.
Compare with a small script run outside the pipelines; results go to a test `outdir`.

| # | Check | Pass criterion |
|---|---|---|
| 1 | `EXTRACT_DE` at default filters vs `DE_results_strong.csv` | Same rows after sorting. If they differ, stop: the May file may predate the current script |
| 2 | All 162 gene-set files vs `data/sensitivity/thr0.5/gene_sets` | Identical **as sets** (finding 4) |
| 3 | `LDSC_CHAIN` vs `sensitivity/thr0.5/ldsc_results`: AITD Rest, ASTCO Stim48hr, OA Rest (single column) | All `.results` numeric fields within 1e-6 relative |
| 4 | Expected non-results: ADFH Rest, SCZ Rest | ADFH no files; SCZ log only |
| 5 | `01_collect_results.py` on the new output | Same `heritability_long.csv` |
| 6 | `sha256` of `matching.py` vs null_v2 | Identical |
| 7 | All Stage 2 gene sets and `*_matching.tsv` for the 5 runs (1,005 chains; no LDSC) | Byte-identical |
| 8 | Sample: `real` + perms 1–20 for all 5 runs (105 chains) vs Stage 2 `.results` | Within 1e-6 |
| 9 | Full wrapper run: all 1,005 chains, `COLLECT` | `.results` within 1e-6; `summary_all.csv`, `pooled_mean_z.csv` within 1e-6 (string comparison would fail across polars versions) |
| 10 | Relaxed run: `min_de_genes = 0`, default others | Completes; no self-edges; rerun k-NN checks; set sizes reviewed before any null |

Check 9 costs ~1,156 CPU-h (~4 h at 40). If checks 7 and 8 pass, it mainly tests the collector and
full-scale scheduling (question 5).

## 9. Effort and risks

| File | Lines | | File | Lines |
|---|---|---|---|---|
| pipeline `main.nf` | 50 | | wrapper `main.nf` | 50 |
| `modules/prepare.nf` | 70 | | `modules/permutation.nf` | 60 |
| `modules/ldsc_chain.nf` | 55 | | `nextflow.config` | 35 |
| `modules/refs.nf` | 25 | | `run_driver.sbatch` | 15 |
| `nextflow.config` | 95 | | `README.md` | 20 |
| `run_driver.sbatch` | 15 | | | |
| `README.md` | 30 | | | |

No new file exceeds ~100 lines. **New: ~520 lines, plus ~40 lines edited in copied scripts**
(extraction ~20, permutation sets ~3, collector ~16).

**Effort:** about 1 week: 2 days main pipeline, 1 day wrapper, 2–3 days validation, mostly cluster
turnaround.

**Alternative: only parameterise the DEG filter in Snakemake.** Add 3 `extract_de.smk` lines and 3
config keys, make `de_results` and `reference_gene_sets` per-variant, and apply the memory fix. About
half a day, ~15 lines, and existing results stay directly comparable. It keeps three copies of the
backend, the plan A estimate of 117 h, and silent filter defaults.

**Recommendation:** do the rewrite if plan A or plan B is likely, or if you expect repeated filter
variants. If you need only one no-75-DEG run with a null for a few traits soon, do the Snakemake
change first; the extraction fix carries over.

| Risk | Mitigation |
|---|---|
| Extraction runs out of memory unfiltered (finding 3) | Filter with numpy masks before building the long table; keep row order; check 1 |
| `DE_results_strong.csv` not reproducible from the current script | Check 1 first; if it differs, decide which is the reference |
| Nextflow on CSD3 unverified: only the `env_nf` name is known, and `moduleDir` access and `--fakeroot` are untested | Half-day smoke test on one chain before anything else |
| `/home` full, grace period ends 2026-11-04 per NOTES | Set `NXF_HOME` on RDS; clear space first |
| `chapter2/` is entirely untracked in git, including null_v2 | Commit it before starting (your call) |
| Two copies of `matching.py` drift | Check 6 in the README; null_v2 copy stays frozen |
| Inode use: ~15 files per task × 11,000 in `work/` | `nextflow clean` after each collected run; RDS inode quota unknown |
| Node-local scratch too small or slow | Test in smoke run; `scratch false` needs no code change |
| k-NN under relaxed filters behaves differently; the AITD Rest skew in NOTES may change | Check 10 before any null |
| Plan A includes SCZ (fails) and sparse traits | Fix or exclude SCZ first; 603 chains would be wasted |
| Exit 3 ignored for null chains hides sparse draws | `n_perm_completed` in the summary reports them |

## 10. Questions for you

1. **Is one `LDSC_CHAIN` process acceptable instead of one process per step?** I assumed yes. It costs
   only step-level resume, at most one chain's re-run. If you need separate processes, use
   `-profile node`, with a cleanup process and resume based on existing `.results`.
2. **Should a no-75-DEG run also change on-target filtering?** I assumed `min_de_genes = 0` with
   `min_cells = 50` and both on-target filters kept.
3. **How large is the next null: plan A, plan B, or a few traits?** This decides rewrite versus
   Snakemake first.
4. **What is the account's concurrent CPU limit on icelake?** `queueSize` 40 needs 320 CPUs.
5. **Is check 9 needed in full, or do checks 7 and 8 plus the collector run on existing Stage 2
   `.results` suffice?** Full is planned, as you asked; it costs ~1,156 CPU-h.
6. **Fix SCZ before any rerun?** Rerunning without a fix reproduces the failure. Fixing it adds
   three fits and changes the Fig. 4 FDR families.
7. **Does `make_gene_sets_tenk10k.py` belong?** I would leave it out, as it reads a different coloc
   source.
8. **Is `approved_drug_enrichment_permutation_v2.py` pipeline or analysis?** I would keep it in the uv
   project until it has been run once.
