#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
squeue -u "$USER"
latest=$(ls -t logs/submission_*.txt 2>/dev/null | head -n 1 || true)
if [[ $# -gt 0 ]]; then
  [[ "$1" =~ ^[0-9]+(_[0-9]+)?$ ]] || { echo 'Supply a real numeric job ID.' >&2; exit 1; }
  sacct -j "$1" --format=JobID,JobName%20,State,ExitCode,Elapsed,MaxRSS
elif [[ -n "$latest" ]]; then
  ids=$(sed -n -E 's/^(SETUP|PILOT|SELECT|TEST|COLLECT)=([0-9]+)$/\2/p' "$latest" | paste -sd, -)
  if [[ -n "$ids" ]]; then sacct -j "$ids" --format=JobID,JobName%20,State,ExitCode,Elapsed,MaxRSS; fi
fi
if [[ -z "${OUT_ROOT:-}" && -n "$latest" ]]; then OUT_ROOT=$(sed -n 's/^OUTPUT=//p' "$latest"); fi
if [[ -n "${OUT_ROOT:-}" && -d "$OUT_ROOT" ]]; then
  echo "OUTPUT=$OUT_ROOT"
  for phase in pilot test; do
    printf '%s completed populations: ' "$phase"
    find "$OUT_ROOT/optimizer_results/$phase" -name COMPLETE.rds 2>/dev/null | wc -l || true
  done
  if [[ -f "$OUT_ROOT/ALL_COMPLETE.rds" ]]; then echo 'All formal results collected.'; fi
fi
