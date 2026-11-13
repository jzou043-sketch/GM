#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p logs
export GM_OUT_ROOT="${GM_OUT_ROOT:-$PWD/results_0248_$(date +%Y%m%d_%H%M%S)_$$}"
job=$(sbatch --parsable --export=ALL --partition="${PARTITION:-epyc}" \
  --time="${JOB_TIME:-7-00:00:00}" submit_replay.sbatch)
job=${job%%;*}
printf 'Job ID: %s\nOutput: %s\nLog: logs/gm0248_%s.out\n' "$job" "$GM_OUT_ROOT" "$job" | tee "logs/submitted_${job}.txt"
echo "One job: four scenarios sequentially, four strategies per scenario."
