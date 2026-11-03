# Running on CSD3

Compares one study × condition with the existing `coloc_thr0.5` results
(`data/sensitivity/thr0.5/`). Default parameters are that variant.

## Setup

Clone onto RDS: `/home` is full, and the container only sees the RDS binds.

```bash
cd /rds/project/rds-csoP2nj6Y6Y/biv22
git clone -b nf-rewrite git@github.com:boris-vasilev/perturb_heritability_model.git nf_rewrite
cd nf_rewrite/chapter2/nf_ldsc_pipeline

export NXF_HOME=/rds/project/rds-csoP2nj6Y6Y/biv22/.nextflow
nextflow -version        # needs >= 24.04; Java 17+
```

Cluster paths and SLURM settings are all in the `csd3` profile of `nextflow.config`.

## 1. One study × condition

Reuses the existing DE table, so this needs one gene-set job and one chain (about 10 min).

```bash
sbatch run_driver.sbatch \
    --de_results data/DE_results_strong.csv \
    --studies AITD_Saevarsdottir_32581359_1 --conditions Rest \
    --outdir $PWD/test_out
```

If the chain fails, read `.command.err` in the work directory named in `nf_driver_<job>.log`.
The likeliest first-run problems are `--fakeroot`, bin paths inside the container, and
node-local scratch (remove `scratch = true` to rule that out).

## 2. Compare with the existing run

```bash
WS=/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_heritability_model
RUN=$(ls -d test_out/*/)
compare() { apptainer exec --fakeroot -B /rds/project/rds-csoP2nj6Y6Y \
    $WS/containers/ldsc39_env.sif python bin/compare_to_reference.py "$@"; }

compare --gene-sets $RUN/gene_sets $WS/data/sensitivity/thr0.5/gene_sets
compare --results $RUN/ldsc_results $WS/data/sensitivity/thr0.5/ldsc_results
```

Each prints `OK: no differences` or the differing file and column, and exits 1 on a difference.
Gene sets are compared as sets, because their line order varies; `.results` numbers within
1e-6 relative.

## 3. Other planned checks

Rerun section 1 with each of these, then section 2 (or `ls` for the non-results).

| Run | Expected |
|---|---|
| `--studies ASTCO_Ferreira_30929738_1 --conditions Stim48hr` | Matches |
| `--studies OA_McDonald_36411363_3 --conditions Rest` | Matches (single-column fit) |
| `--studies ADFH_Jansen_30617256_1 --conditions Rest` | No files in `ldsc_results`; run ends OK with one ignored task |
| `--studies SCZ_Pardiñas_29483656_1 --conditions Rest` | `.log` only, ending in "0 SNPs remain" |

## 4. Extraction

Omit `--de_results` to rerun extraction (up to 2 h), then compare with the existing table:

```bash
compare --de $RUN/DE_results.csv $WS/data/DE_results_strong.csv
```

If this differs, stop: the May table may predate the current script.

To run without the DEG filter, add `--min_de_genes 0`. Outputs go to a separately named folder.

## 5. Permutation wrapper (`../nf_ldsc_permutation`)

Checks it against Stage 2 (`data/permutations/heritability/`, read-only) in three steps, cheapest
first. Run from `nf_ldsc_permutation/`, with `WS` as in section 2. The options below reproduce
Stage 2. `--min_de_genes 75` only names the output folder after the table's filters.

```bash
cd ../nf_ldsc_permutation
SIF=$WS/containers/ldsc39_env.sif
py() { apptainer exec --fakeroot -B /rds/project/rds-csoP2nj6Y6Y $SIF python "$@"; }
STAGE2="--de_results data/DE_results_strong.csv --min_de_genes 75 \
        --reference_gene_sets data/sensitivity/thr0.5/gene_sets"
REF=$WS/data/permutations/heritability
```

**5a. Gene sets only** for the 5 Stage 2 runs: 5 jobs, no LDSC. The real sets are checked against
`thr0.5` inside the jobs. Then every gene-set file must match Stage 2 byte for byte.

```bash
sbatch run_driver.sbatch $STAGE2 --gene_sets_only
# once the driver job has finished:
RUN=$(ls -d $WS/data/permutations/deg75_*/)
cd $RUN && find . -path '*/gene_sets/*' | wc -l       # expect 4020 (1,005 chains x 4)
n=0; for d in */perm*/gene_sets; do diff -rq "$d" "$REF/$d" || n=$((n+1)); done; echo "$n differing"
cd -
```

**5b. Collector on the existing Stage 2 `.results`**. Summary numbers must be within 1e-6 of
`summary_all.csv`. Stage 2 has no `.combine.txt` files, so the new `L2_1` handling cannot change
anything here. COLLECT needs matplotlib in the container. If it is missing, run the collector with
`uv run --project ../permutation_analysis python` instead.

```bash
RUNS="AITD_Saevarsdottir_32581359_1:Rest,Stim8hr,Stim48hr ASTCO_Ferreira_30929738_1:Rest,Stim48hr"
py bin/collect_permutation_ldsc.py --root $REF --out-dir $PWD/collect_check --runs $RUNS
py -c "
import sys, polars as pl; sys.path.insert(0, '../nf_ldsc_pipeline/bin')
from compare_to_reference import frames_differ
for f in ['summary_all.csv', 'pooled_mean_z.csv', 'matching_log_summary.csv']:
    a, b = (pl.read_csv(d + '/' + f, infer_schema_length=None) for d in ['collect_check', '$REF'])
    print(f, frames_differ(a, b, 1e-6) or 'OK')"
```

**5c. One real chain**: AITD Rest, perm `real` only, compared with its Stage 2 `.results`.

```bash
sbatch run_driver.sbatch $STAGE2 --studies AITD_Saevarsdottir_32581359_1 --conditions Rest \
    --n_perm 0 --perm_root data/permutations/nf_test
# once the driver job has finished:
S=AITD_Saevarsdottir_32581359_1/permreal/ldsc_results
py ../nf_ldsc_pipeline/bin/compare_to_reference.py \
    --results $(ls -d $WS/data/permutations/nf_test/deg75_*/)$S $REF/$S
```

If all three pass, a sample of perms with LDSC (PLAN.md section 8, check 8) is next.
