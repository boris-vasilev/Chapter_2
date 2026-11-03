# nf_ldsc_permutation

Matched-null S-LDSC permutation (the Nextflow port of `permutation_analysis/old_snakemake/`).
For each study × condition in `params.runs` it builds the real coloc gene sets and `n_perm`
matched pseudo-anchor sets, runs each through `LDSC_CHAIN` from `../nf_ldsc_pipeline/`, and
collects the null summaries. Design: Part 2 of `../nf_ldsc_pipeline/PLAN.md`.

```
main.nf                  workflow; includes EXTRACT_DE, LDSC_CHAIN, refs from ../nf_ldsc_pipeline
modules/permutation.nf   BUILD_EXPRESSION, PERM_GENE_SETS (one job per study × condition), COLLECT
nextflow.config          permutation parameters; includes ../nf_ldsc_pipeline/nextflow.config
bin/                     matching.py (byte-identical to Stage 2) and the scripts that use it
```

Outputs go to `data/permutations/<tag>/`, where the tag is the main pipeline's filter variant
plus the matching options. The layout is the same as Stage 2:
`{study}/perm{P}/{gene_sets,ldsc_results}/`, `{study}_{cond}/` summaries, and `summary_all.csv`.
`data/permutations/heritability/` (Stage 2) is never written to.

Run on CSD3 from this folder with `sbatch run_driver.sbatch [options]`. `--studies` and
`--conditions` select a subset of `runs`. To add runs, edit `runs` in `nextflow.config`
(`-params-file` merges maps, so it cannot remove entries). Validation steps are in
`../nf_ldsc_pipeline/CSD3_RUNBOOK.md`.
