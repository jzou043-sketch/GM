#!/usr/bin/env bash
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
squeue -u "$USER"
if [[ $# -gt 0 ]]; then
  sacct -j "$1" --format=JobID,JobName%20,State,ExitCode,Elapsed,MaxRSS
  for f in "logs/gm0248_$1.out" "logs/gm0248_$1.err"; do
    if [[ -f "$f" ]]; then printf '\n%s\n' "$f"; tail -n 40 "$f"; fi
  done
fi
