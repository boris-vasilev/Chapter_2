#!/bin/bash
# Snakemake driver for the v2 permutation null, run as its own slurm job.
# Submit from perturb_heritability_null_v2/:   sbatch run_driver.sbatch
# Chain jobs are submitted by the slurm executor (existing profile); this job
# only schedules them. Resubmitting resumes: finished .results are skipped.
#SBATCH -J null_v2_driver
#SBATCH -A MRC-BSU2-SL2-CPU
#SBATCH -p icelake
#SBATCH -c 1
#SBATCH --mem=4000
#SBATCH -t 20:00:00
#SBATCH -o /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_heritability_model/data/permutations/heritability/driver_%j.log

set -euo pipefail
cd /rds/project/rds-csoP2nj6Y6Y/biv22/perturb_heritability_model/perturb_heritability_null_v2

# --rerun-triggers mtime: only missing/out-of-date .results are (re)run, so
# completed chains (incl. the kept AITD Rest real chain) are never redone
# because of a param/code-hash difference.
~/.local/bin/snakemake --profile ../perturb_ldsc_pipeline/slurm \
    --rerun-triggers mtime --rerun-incomplete
