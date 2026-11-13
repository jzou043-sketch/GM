#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if command -v squeue >/dev/null 2>&1; then squeue -u "$USER"; fi
latest=$(ls -t logs/submission_*.txt logs/direct_*.txt 2>/dev/null | head -n1 || true)
if [[ -n "$latest" ]]; then
 ids=$(sed -n -E 's/^(SETUP|RUN|COLLECT)=([0-9]+)$/\2/p' "$latest" | paste -sd, -)
 if [[ -n "$ids" ]] && command -v sacct >/dev/null 2>&1; then
  sacct -j "$ids" --format=JobID,JobName%20,State,ExitCode,Elapsed,MaxRSS
 fi
 out=$(sed -n 's/^OUTPUT=//p' "$latest")
 echo "OUTPUT=$out"
 for scale in N100 N200 N300 N400 N500 N600 smoke; do
  if [[ -d "$out/runs/$scale" ]]; then
   printf '%s completed: ' "$scale"
   find "$out/runs/$scale" -name COMPLETE.rds | wc -l
  fi
 done
 [[ ! -f "$out/ALL_COMPLETE.rds" ]] || echo 'All results collected: summary/ and figures/ are ready.'
fi
