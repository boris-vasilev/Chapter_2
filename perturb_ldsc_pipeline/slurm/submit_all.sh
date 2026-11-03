#!/bin/bash
# Submit one sbatch job per study listed in the active module's
# GWAS_study_info.csv. Each job runs the full pipeline for that
# study × 3 conditions.
#
# Usage:
#   bash slurm/submit_all.sh                 # submit everything (coloc by default)
#   MODULE=l2g bash slurm/submit_all.sh      # switch to L2G
#   bash slurm/submit_all.sh --dry-run       # print what would be submitted
#
# Override the study list explicitly if needed:
#   STUDY_INFO=/path/to/custom.csv MODULE=l2g bash slurm/submit_all.sh

set -eo pipefail

DRY=0
[[ "${1:-}" == "--dry-run" ]] && DRY=1

WORKSPACE="/rds/project/rds-csoP2nj6Y6Y/biv22/perturb_heritability_model"
CONDITIONS=(Rest Stim8hr Stim48hr)
CONFIG_YAML="$WORKSPACE/perturb_ldsc_pipeline/config.yaml"

# Module — default to whatever config.yaml says, so this script tracks
# the configured pipeline without needing MODULE=… on every invocation.
# Env override still wins: MODULE=l2g bash slurm/submit_all.sh.
CONFIG_MODULE=$(awk -F'[:#]' '$1=="gene_set_module" {
                                  gsub(/[" \047]/, "", $2); print $2; exit
                              }' "$CONFIG_YAML")
MODULE="${MODULE:-${CONFIG_MODULE:-coloc}}"

# Per-module default study list. Override with $STUDY_INFO if needed.
# Anything coloc-derived (coloc, coloc_strict, coloc_strong, …) uses the
# same curated study list; only L2G has its own.
case "$MODULE" in
    l2g)   STUDY_INFO_DEFAULT="$WORKSPACE/data/opentargets/GWAS_study_info.csv" ;;
    *)     STUDY_INFO_DEFAULT="$WORKSPACE/data/GWAS_study_info.csv" ;;
esac
STUDY_INFO="${STUDY_INFO:-$STUDY_INFO_DEFAULT}"

if [[ ! -f "$STUDY_INFO" ]]; then
    echo "ERROR: study info not found at $STUDY_INFO" >&2
    exit 1
fi

echo "Module:     $MODULE"
echo "Study info: $STUDY_INFO"

# Read study column (first col), skip header.
while IFS=, read -r study _; do
    [[ -z "$study" || "$study" == "study" ]] && continue

    targets=()
    for cond in "${CONDITIONS[@]}"; do
        targets+=("$WORKSPACE/data/ldsc_results/$MODULE/$study/${cond}_combined.results")
    done

    if (( DRY )); then
        echo "[dry-run] sbatch --export=ALL,MODULE=$MODULE --job-name=ldsc_${study} slurm/run.sbatch ${targets[*]}"
    else
        sbatch --export="ALL,MODULE=$MODULE" \
               --job-name="ldsc_${study}" \
               slurm/run.sbatch "${targets[@]}"
        sleep 1   # one sbatch per second
    fi
done < "$STUDY_INFO"
