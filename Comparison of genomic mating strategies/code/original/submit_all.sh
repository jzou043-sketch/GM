#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export GM_ROOT="$PWD"
export GM_PROFILE="${GM_PROFILE:-inbreeding_restored}"
[[ "$GM_PROFILE" == inbreeding_restored ]] || { echo 'Invalid GM_PROFILE'; exit 1; }
export OUT_ROOT="${OUT_ROOT:-$PWD/results_inbreeding_$(date +%Y%m%d_%H%M%S)_$$}"
[[ "$OUT_ROOT" = /* ]] || { echo 'OUT_ROOT must be absolute'; exit 1; }
if [[ -d "$OUT_ROOT" && -n "$(ls -A "$OUT_ROOT")" ]]; then
 echo 'Output directory is not empty. Unset OUT_ROOT for a fresh experiment; use resume_failed.sh to resume.'
 exit 1
fi
mkdir -p logs "$OUT_ROOT"
part="${PARTITION:-intel}";concurrency="${ARRAY_CONCURRENCY:-10}"
[[ "$concurrency" =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid concurrency'; exit 1; }
count=600
if [[ "${SMOKE:-0}" == 1 ]]; then count=1; fi
setup=$(sbatch --parsable --partition="$part" --job-name=gm6i_setup --time="${SETUP_TIME:-1-00:00:00}" --export=ALL,STAGE=init worker.sbatch)
setup=${setup%%;*}
run=$(sbatch --parsable --partition="$part" --job-name=gm6i_run --time="${JOB_TIME:-2-00:00:00}" --array="1-$count%$concurrency" --dependency="afterok:$setup" --export=ALL,STAGE=run worker.sbatch)
run=${run%%;*}
collect=$(sbatch --parsable --partition="$part" --job-name=gm6i_collect --time=04:00:00 --dependency="afterok:$run" --export=ALL,STAGE=collect worker.sbatch)
collect=${collect%%;*}
printf 'OUTPUT=%s\nSETUP=%s\nRUN=%s\nCOLLECT=%s\nPROFILE=%s\n' "$OUT_ROOT" "$setup" "$run" "$collect" "$GM_PROFILE" | tee "logs/submission_$setup.txt"
echo "Submitted $count scale/replicate tasks; each runs 3 methods sequentially."
echo 'Production: 100 maximum populations, 6 nested scales each, 1800 method results.'
