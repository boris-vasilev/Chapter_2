# S-LDSC matched-anchor permutation null, v2

Corrected version of `perturb_ldsc_permutation/` (the pilot). Same question:
is the S-LDSC enrichment of coloc-anchored regulator/target sets specific to
the coloc anchors, or does expanding *any* matched anchor set through the CD4
Perturb-seq graph give the same enrichment? Fixes and measured facts are in
`NOTES_null_v2.md`.

## Layout
```
config.yaml                  runs (study -> conditions), n_perm, matching, paths
Snakefile                    verify_real_sets (local) + perm_chain (1 slurm job per perm x cond)
env.sh                       UV_CACHE_DIR / UV_PYTHON_INSTALL_DIR on RDS (/home is full)
pyproject.toml, uv.lock      uv env for everything outside the LDSC container
scripts/
  matching.py                    shared matched sampler (container-safe: polars+numpy)
  make_permutation_gene_sets.py  one perm x cond's gene sets (container)
  build_screen_expression.py     per-condition screen baseMean table (uv, once)
  matching_check.py              Stage 1: OLD vs FIXED matching, no LDSC (uv)
  collect_permutation_ldsc.py    per-run summaries + pooled tests (uv)
  plot_permutation_null.py       null histograms (imported by collect)
  compare_with_pilot.py          old_vs_new_comparison.csv (uv)
  estimate_costs.py              plan_costs.csv (uv)
  approved_drug_enrichment_permutation_v2.py  drug-target permutation, v2 matching (not run)
```

Outputs go to `../data/permutations/` (see NOTES for the tree).

## Run
```bash
cd perturb_heritability_null_v2
source env.sh
uv run python scripts/build_screen_expression.py \
    --h5ad ../data/GWCD4i.DE_stats.h5ad \
    --out ../data/permutations/inputs/screen_expression.tsv      # once
uv run python scripts/matching_check.py --out-dir ../data/permutations/matching_check

# Driver (uses snakemake 9.20 from the env_nf conda env + slurm executor plugin)
tmux new -s null_v2 'snakemake --profile ../perturb_ldsc_pipeline/slurm 2>&1 | tee ../logs/null_v2_snakemake.log'

# Any time, over whatever has completed:
uv run python scripts/collect_permutation_ldsc.py --config config.yaml
uv run python scripts/compare_with_pilot.py --config config.yaml
```
