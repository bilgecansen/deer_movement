#!/usr/bin/env bash
# Iterate every deer with a fitted GAM model in results/ and run
# scripts/gam/run_sims_gam.R on each. By default, skip deer whose simulation
# output already exists; pass --overwrite to reprocess everything.
#
# Failures in Rscript do not halt the loop — we capture the exit code and
# continue. Final summary prints counts and elapsed time.
#
# Usage:
#   bash scripts/gam/run_sims_gam_all.sh              # resumable (default)
#   bash scripts/gam/run_sims_gam_all.sh --overwrite  # reprocess all

shopt -s nullglob

overwrite=false
if [[ "${1:-}" == "--overwrite" ]]; then
  overwrite=true
fi

n_done=0
n_skipped=0
n_failed=0
start=$(date +%s)

for f in results/gam/results_gam_*.rds; do
  base=$(basename "$f" .rds)
  key=${base#results_gam_}
  IFS='_' read -r id season year <<< "$key"

  out_path="sims/gam/sims_gam_${key}.rds"

  if [[ "$overwrite" == false && -f "$out_path" ]]; then
    echo "[skip] $key"
    n_skipped=$((n_skipped + 1))
    continue
  fi

  echo "[run]  $key"
  if Rscript scripts/gam/run_sims_gam.R "$id" "$season" "$year"; then
    n_done=$((n_done + 1))
  else
    rc=$?
    echo "[fail] $key (exit $rc)"
    n_failed=$((n_failed + 1))
  fi
done

elapsed=$(( $(date +%s) - start ))
mins=$(( elapsed / 60 ))
secs=$(( elapsed % 60 ))
echo
echo "Done: $n_done   Skipped: $n_skipped   Failed: $n_failed   Elapsed: ${mins}m ${secs}s"
