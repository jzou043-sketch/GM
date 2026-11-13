#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/gm6_%x_%A_%a.out
#SBATCH --error=logs/gm6_%x_%A_%a.err
set -euo pipefail
cd "${GM_ROOT:-${SLURM_SUBMIT_DIR:?}}"
if [[ "${SKIP_MODULE_LOAD:-0}" != 1 ]] && type module >/dev/null 2>&1; then
 module purge
 module load "${R_MODULE:-R/4.5.0}"
fi
command -v Rscript >/dev/null || { echo 'Rscript unavailable; set R_MODULE'; exit 1; }
export R_LIBS_USER="${GM6_R_LIB:-$HOME/R/gm6_formula_lib}"
mkdir -p "$R_LIBS_USER" "${OUT_ROOT:?}/locks"
for lib in "${PREVIOUS_R_LIB:-$HOME/R/gm3_50k_lib}" "$HOME/R/gm_package_osnb_benchmark_lib"; do
 if [[ -d "$lib" ]]; then export R_LIBS="$lib${R_LIBS:+:$R_LIBS}"; fi
done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
command -v flock >/dev/null || { echo 'Linux flock is required'; exit 1; }
exec 9>"$OUT_ROOT/locks/${STAGE:?}_${SLURM_ARRAY_TASK_ID:-0}.lock"
flock -n 9 || { echo 'Task already running'; exit 1; }
case "$STAGE" in
 init) chmod u+x vendor/AlphaMate_Unix; Rscript setup.R; Rscript driver.R init ;;
 run)
  task="${SLURM_ARRAY_TASK_ID:?}"
  rep=$(awk -F, -v id="$task" 'NR==id+1 {gsub(/"/,"",$1);print $1}' "$OUT_ROOT/task_map.csv")
  [[ "$rep" =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid task/replicate mapping'; exit 1; }
  flock "$OUT_ROOT/locks/base_rep_$rep.lock" Rscript driver.R prepare "$rep"
  Rscript driver.R run "$task"
  ;;
 collect) Rscript driver.R collect ;;
 *) echo 'Unknown stage'; exit 1 ;;
esac
