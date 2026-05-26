#!/bin/bash
#SBATCH --job-name=mg_mean_clean
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=5-00:00:00
#SBATCH --output=mg_mean_clean_%j.out
#SBATCH --error=mg_mean_clean_%j.err

set -euo pipefail

module purge
module load R/4.5.0

# Run in the directory where this sbatch file is submitted.
cd "${SLURM_SUBMIT_DIR}"

# Default simulation settings. These can be overridden at submission time.
export N_REP=${N_REP:-100}
export N_GENERATIONS=${N_GENERATIONS:-10}
export N_FEMALE=${N_FEMALE:-100}
export N_MALE=${N_MALE:-100}
export N_SNP=${N_SNP:-1000}
export N_QTL=${N_QTL:-200}
export N_FEMALE_SELECT=${N_FEMALE_SELECT:-12}
export N_MALE_SELECT=${N_MALE_SELECT:-6}
export GS_GM_PRESELECT_FEMALE=${GS_GM_PRESELECT_FEMALE:-36}
export GS_GM_PRESELECT_MALE=${GS_GM_PRESELECT_MALE:-18}
export N_RESTART=${N_RESTART:-2}
export N_OFFSPRING_PER_CROSS=${N_OFFSPRING_PER_CROSS:-20}
export H2=${H2:-0.50}
export MU=${MU:-50}
export SCENARIO=${SCENARIO:-S3_strong}
export N_WORKERS=${N_WORKERS:-8}
export OUT_DIR=${OUT_DIR:-multigen_expected_mean_large_results}

# Limit BLAS/OpenMP threads inside each R worker.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "Job started at: $(date)"
echo "Node: $(hostname)"
echo "Working directory: $(pwd)"
echo "R script: simulate_100rep_four_strategies_expected_mean_large.R"
echo "N_REP=${N_REP}, N_GENERATIONS=${N_GENERATIONS}, N_WORKERS=${N_WORKERS}, SCENARIO=${SCENARIO}, N_RESTART=${N_RESTART}"
echo "OUT_DIR=${OUT_DIR}"

Rscript simulate_100rep_four_strategies_expected_mean_large.R

echo "Job finished at: $(date)"
