#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export N_REP="${N_REP:-1000}" N_PILOT="${N_PILOT:-30}" REPS_PER_TASK="${REPS_PER_TASK:-10}"
export SEED_VALUE="${SEED_VALUE:-20261027}"
for value in "$N_REP" "$N_PILOT" "$REPS_PER_TASK" "$SEED_VALUE"; do
  [[ "$value" =~ ^[1-9][0-9]*$ ]] || { echo "Positive integers required" >&2; exit 1; }
done
if [[ -n "${EVAL_BUDGET:-}${MAX_CONFIGS:-}${TUNING_CONFIGS:-}" ]]; then
  echo "NOTE: old EVAL_BUDGET/MAX_CONFIGS/TUNING_CONFIGS are NOT used by this package."
fi
unset EVAL_BUDGET MAX_CONFIGS TUNING_CONFIGS
export OUT_ROOT="${OUT_ROOT:-$PWD/package_osnb_run_$(date +%Y%m%d_%H%M%S)_$$}"
if [[ -e "$OUT_ROOT/configuration/config.rds" ]]; then
  echo "Use a new OUT_ROOT for submit_all.sh; resume individual stages separately." >&2
  exit 1
fi
mkdir -p logs "$OUT_ROOT"
tasks=$(( (N_REP + REPS_PER_TASK - 1) / REPS_PER_TASK ))
part="${PARTITION:-epyc}"
setup=$(sbatch --parsable --partition="$part" --job-name=osnb6_setup --time="${SETUP_TIME:-0-04:00:00}" --export=ALL,STAGE=init worker.sbatch)
setup=${setup%%;*}
pilot=$(sbatch --parsable --partition="$part" --job-name=osnb6_pilot --time="${PILOT_TIME:-3-00:00:00}" \
  --array="1-${N_PILOT}%${PILOT_CONCURRENCY:-4}" --dependency="afterok:$setup" --export=ALL,STAGE=pilot worker.sbatch)
pilot=${pilot%%;*}
select=$(sbatch --parsable --partition="$part" --job-name=osnb6_select --time=01:00:00 \
  --dependency="afterok:$pilot" --export=ALL,STAGE=select worker.sbatch)
select=${select%%;*}
test=$(sbatch --parsable --partition="$part" --job-name=osnb6_test --time="${JOB_TIME:-7-00:00:00}" \
  --array="1-${tasks}%${ARRAY_CONCURRENCY:-10}" --dependency="afterok:$select" --export=ALL,STAGE=test worker.sbatch)
test=${test%%;*}
collect=$(sbatch --parsable --partition="$part" --job-name=osnb6_collect --time=04:00:00 \
  --dependency="afterok:$test" --export=ALL,STAGE=collect worker.sbatch)
collect=${collect%%;*}
printf 'OUTPUT=%s\nSETUP=%s\nPILOT=%s\nSELECT=%s\nTEST=%s\nCOLLECT=%s\n' \
  "$OUT_ROOT" "$setup" "$pilot" "$select" "$test" "$collect" | tee "logs/submission_${setup}.txt"
echo "Formal populations: $N_REP; separate pilot populations: $N_PILOT; SIX methods, no HC."
echo "Original OSNB unchanged: pool=200 (actual 100), n_restart=2, best search."
echo "Other methods receive its measured request allowance; no OSNB truncation."
echo 'SA: stats::optim(method="SANN"); Tabu: tabuSearch::tabuSearch. No custom SA/Tabu search loops.'
