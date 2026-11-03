# Matched-null permutation: analysis and figures

Analysis around the S-LDSC matched-anchor null, formerly `perturb_heritability_null_v2/`.
The pipeline itself (gene sets, LDSC chains, collector) is now `../nf_ldsc_permutation/`. Its
`bin/matching.py` is the one copy of the matched sampler, and the scripts here import it from
there. Fixes and measured facts are in `NOTES_null_v2.md`.

## Layout
```
scripts/
  matching_check.py              Stage 1: old vs fixed matching, no LDSC
  matching_check_knn.py          Stage 1 re-check with k-NN matching
  approved_drug_enrichment_permutation_v2.py  drug-target permutation, v2 matching (not run)
  archive/                       target_size_check, compare_with_pilot, estimate_costs (one-offs)
figures/                         Fig. 5 and supplementary figures (read data/permutations/heritability)
old_snakemake/                   Snakefile, config.yaml, run_driver.sbatch: kept until the
                                 Nextflow port is validated. They do not run from here as is,
                                 because the scripts moved; run them from commit 79cfbc8 instead.
env.sh, pyproject.toml, uv.lock  uv environment for everything here
```

## Run
```bash
cd permutation_analysis
source env.sh
uv run python scripts/matching_check.py --out-dir ../data/permutations/matching_check
```
