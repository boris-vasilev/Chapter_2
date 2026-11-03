#!/bin/bash
#SBATCH --job-name=ldsc_pipeline
#SBATCH -A MRC-BSU2-SL2-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=22
#SBATCH --time=1:30:00
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
#
# Runs the whole Snakemake pipeline inside a single allocation.
#
# Works both as:
#   sbatch slurm/run.sbatch [snakemake-args...]      # batch job
#   bash   slurm/run.sbatch [snakemake-args...]      # inside an salloc'd shell
#
# Examples:
#   bash slurm/run.sbatch -n                          # dry run
#   bash slurm/run.sbatch --until calculate_ld_scores # stop early
#   bash slurm/run.sbatch data/ldsc_results/UC_DeLange_28067908_1/Rest_combined.results
#

source ~/.bash_profile
set -eo pipefail

# cwd — resolve symlinks so paths match Path.resolve() inside the Snakefile.
WORKDIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$(realpath "$WORKDIR")"
mkdir -p logs

# Cores — use slurm allocation if available, else fall back.
CORES="${SLURM_CPUS_PER_TASK:-8}"

# Pick $MODULE if set, so submit_all.sh and interactive runs can pick
# the gene-set module without editing config.yaml.
CONFIG_OVERRIDES=()
if [[ -n "${MODULE:-}" ]]; then
    CONFIG_OVERRIDES+=(--config "gene_set_module=$MODULE")
    echo "Snakemake config override: gene_set_module=$MODULE"
fi

snakemake \
    "${CONFIG_OVERRIDES[@]}" \
    --cores "$CORES" \
    --software-deployment-method apptainer \
    --apptainer-args "--fakeroot \
        -B /rds/project/rds-csoP2nj6Y6Y:/rds/project/rds-csoP2nj6Y6Y \
        -B /home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y:/home/biv22/rds/rds-mrc-bsu-csoP2nj6Y6Y" \
    --keep-going \
    --rerun-incomplete \
    --latency-wait 60 \
    "$@"
