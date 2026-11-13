#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export GM_ROOT="$PWD"
if [[ -z "${OUT_ROOT:-}" ]]; then
 latest=$(ls -t logs/submission_*.txt 2>/dev/null | head -n1 || true)
 [[ -n "$latest" ]] || { echo 'Set OUT_ROOT to the experiment to resume'; exit 1; }
 export OUT_ROOT=$(sed -n 's/^OUTPUT=//p' "$latest")
fi
[[ -f "$OUT_ROOT/configuration.rds" ]] || { echo 'Setup did not complete; fix setup and start in a new output directory'; exit 1; }
command -v Rscript >/dev/null || { echo 'Load your R module first'; exit 1; }
ids=$(Rscript -e 'o<-Sys.getenv("OUT_ROOT");t<-read.csv(file.path(o,"task_map.csv"));ok<-file.exists(file.path(o,"runs",t$Scale,sprintf("rep_%03d",t$Replicate),"COMPLETE.rds"));cat(paste(t$Task[!ok],collapse=","))')
part="${PARTITION:-intel}";concurrency="${ARRAY_CONCURRENCY:-10}"
[[ "$concurrency" =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid concurrency'; exit 1; }
mkdir -p logs
if [[ -z "$ids" ]]; then
 collect=$(sbatch --parsable --partition="$part" --job-name=gm6i_collect --time=04:00:00 --export=ALL,STAGE=collect worker.sbatch)
 printf 'OUTPUT=%s\nCOLLECT=%s\n' "$OUT_ROOT" "${collect%%;*}" | tee "logs/submission_collect_${collect%%;*}.txt"
else
 run=$(sbatch --parsable --partition="$part" --job-name=gm6i_resume --time="${JOB_TIME:-2-00:00:00}" --array="$ids%$concurrency" --export=ALL,STAGE=run worker.sbatch)
 run=${run%%;*}
 collect=$(sbatch --parsable --partition="$part" --job-name=gm6i_collect --time=04:00:00 --dependency="afterok:$run" --export=ALL,STAGE=collect worker.sbatch)
 printf 'OUTPUT=%s\nRUN=%s\nCOLLECT=%s\n' "$OUT_ROOT" "$run" "${collect%%;*}" | tee "logs/submission_resume_$run.txt"
fi
