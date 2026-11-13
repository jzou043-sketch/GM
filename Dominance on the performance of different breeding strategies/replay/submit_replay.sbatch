#!/usr/bin/env bash
#SBATCH --job-name=gm0248
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/gm0248_%j.out
#SBATCH --error=logs/gm0248_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
if [[ "${SKIP_MODULE_LOAD:-0}" != 1 ]] && type module >/dev/null 2>&1; then
  module purge
  module load "${R_MODULE:-R/4.5.0}"
fi
export R_LIBS_USER="${GM_R_LIBS:-$HOME/R/gm_alphasimr_lib}"
mkdir -p "$R_LIBS_USER"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export GM_OUT_ROOT="${GM_OUT_ROOT:-$PWD/results_0248_${SLURM_JOB_ID}}"
mkdir -p "$GM_OUT_ROOT"
# Prevent concurrent jobs from writing into one output directory.
exec 9>"$GM_OUT_ROOT/.run.lock"
flock -n 9 || { echo "Another job is using this output folder." >&2; exit 1; }
trap 'echo "FAILED: exit status $?; inspect logs and preserved output: $GM_OUT_ROOT" >&2' ERR
echo "START $(date) host=$(hostname) output=$GM_OUT_ROOT"
Rscript preflight.R
Rscript run_replay.R
Rscript plot_results.R
archive="${GM_OUT_ROOT}_job${SLURM_JOB_ID}.tar.gz"
if [[ -e "$archive" ]]; then echo "Archive already exists: $archive" >&2; exit 1; fi
tar -czf "$archive" -C "$(dirname "$GM_OUT_ROOT")" "$(basename "$GM_OUT_ROOT")"
echo "ALL SUCCESS $(date)"
echo "Results archive: $archive"
