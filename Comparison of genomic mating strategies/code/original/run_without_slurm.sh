#!/usr/bin/env bash
# Alternative for a Linux server without sbatch. No recursive/background daemon.
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export GM_ROOT="$PWD"
export GM_PROFILE="${GM_PROFILE:-inbreeding_restored}"
[[ "$GM_PROFILE" == inbreeding_restored ]] || { echo 'Unset old GM_PROFILE'; exit 1; }
export OUT_ROOT="${OUT_ROOT:-$PWD/results_inbreeding_$(date +%Y%m%d_%H%M%S)_$$}"
[[ "$OUT_ROOT" = /* ]] || { echo 'OUT_ROOT must be absolute'; exit 1; }
if [[ -d "$OUT_ROOT" && -n "$(ls -A "$OUT_ROOT")" ]]; then
 echo 'Use a new output directory for this driver'; exit 1
fi
mkdir -p logs "$OUT_ROOT"
export SKIP_MODULE_LOAD=1
concurrency="${ARRAY_CONCURRENCY:-4}"
[[ "$concurrency" =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid concurrency'; exit 1; }
count=600
if [[ "${SMOKE:-0}" == 1 ]]; then count=1; fi
printf 'OUTPUT=%s\nPROFILE=%s\n' "$OUT_ROOT" "$GM_PROFILE" | tee "logs/direct_$(date +%Y%m%d_%H%M%S).txt"
env STAGE=init bash worker.sbatch
seq 1 "$count" | xargs -P "$concurrency" -I '{}' env STAGE=run SLURM_ARRAY_TASK_ID='{}' bash worker.sbatch
env STAGE=collect bash worker.sbatch
