#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/osnb6_%x_%A_%a.out
#SBATCH --error=logs/osnb6_%x_%A_%a.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
if [[ "${SKIP_MODULE_LOAD:-0}" != 1 ]] && type module >/dev/null 2>&1; then
  module purge
  module load "${R_MODULE:-R/4.5.0}"
fi
export R_LIBS_USER="${GM_R_LIBS:-$HOME/R/gm_package_osnb_benchmark_lib}"
mkdir -p "$R_LIBS_USER" "${OUT_ROOT:?}" "$OUT_ROOT/locks"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
stage="${STAGE:?}"; index="${SLURM_ARRAY_TASK_ID:-0}"
exec 9>"$OUT_ROOT/locks/${stage}_${index}.lock"
flock -n 9 || { echo "Stage/task already running" >&2; exit 1; }
trap 'echo "FAILED stage=$STAGE task=${SLURM_ARRAY_TASK_ID:-0} exit=$?; partial files retained" >&2' ERR
echo "START $(date) host=$(hostname) stage=$stage task=$index OUT_ROOT=$OUT_ROOT"
case "$stage" in
  init)
    Rscript preflight.R
    Rscript tests.R
    Rscript driver.R init
    ;;
  pilot|test) Rscript driver.R "$stage" "$index" ;;
  select) Rscript driver.R select ;;
  collect)
    Rscript driver.R collect
    Rscript verify_outputs.R "$OUT_ROOT"
    data_archive="${OUT_ROOT}_simulation_data.tar.gz"
    results_archive="${OUT_ROOT}_results.tar.gz"
    [[ ! -e "$data_archive" && ! -e "$results_archive" ]] || { echo "Archives exist; not overwriting" >&2; exit 1; }
    tar -czf "${data_archive}.partial" -C "$OUT_ROOT" simulation_data
    mv -n "${data_archive}.partial" "$data_archive"
    tar -czf "${results_archive}.partial" -C "$OUT_ROOT" optimizer_results summary configuration validation ALL_COMPLETE.rds
    mv -n "${results_archive}.partial" "$results_archive"
    echo "DATA: $data_archive"
    echo "RESULTS: $results_archive"
    ;;
  *) echo "Invalid stage" >&2; exit 1 ;;
esac
echo "SUCCESS $(date) stage=$stage task=$index"
