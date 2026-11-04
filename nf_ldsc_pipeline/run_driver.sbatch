#!/bin/bash
# Nextflow driver as its own SLURM job; it only submits and tracks chain jobs.
# Run from nf_ldsc_pipeline/:  sbatch run_driver.sbatch [nextflow run options]
#SBATCH -J nf_ldsc_driver
#SBATCH -A MRC-BSU2-SL2-CPU
#SBATCH -p icelake
#SBATCH -c 1
#SBATCH --mem=8000
#SBATCH -t 24:00:00
#SBATCH -o nf_driver_%j.log

set +u
source ~/.bash_profile
set -euo pipefail

# /home is full, so Nextflow's own files go on RDS.
export NXF_HOME=/rds/project/rds-csoP2nj6Y6Y/biv22/.nextflow
export NXF_OPTS='-Xms1g -Xmx6g'
export LC_ALL=C.UTF-8 LANG=C.UTF-8    # the JVM must read SCZ_Pardiñas

uv run nextflow run main.nf -profile csd3 -resume "$@"
