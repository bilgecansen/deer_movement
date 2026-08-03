#!/usr/bin/env bash
# Iterate every deer with a fitted GAM model in results/ and run
# scripts/gam/run_logscore_gam.R on each. By default, skip deer whose log-score
# output already exists; pass --overwrite to reprocess everything.
#
# Failures in Rscript do not halt the loop. Final summary prints counts and
# elapsed time.
#
# Usage:
#   bash scripts/gam/run_logscore_gam_all.sh              # resumable (default)
#   bash scripts/gam/run_logscore_gam_all.sh --overwrite  # reprocess all

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

  # results_gam_null_<key>.rds also matches this glob. The null is not a deer of
  # its own — run_logscore_gam.R loads it alongside the numbered models — so skip
  # it silently here rather than invoking R on a bogus key.
  [[ "$base" == results_gam_null_* ]] && continue

  key=${base#results_gam_}
  IFS='_' read -r id season year <<< "$key"

  # One invocation writes both files; only skip when both are already present.
  out_path="filters/gam/logscore_gam_${key}.rds"
  null_out="filters/gam/logscore_gam_null_${key}.rds"

  if [[ "$overwrite" == false && -f "$out_path" && -f "$null_out" ]]; then
    echo "[skip] $key"
    n_skipped=$((n_skipped + 1))
    continue
  fi

  echo "[run]  $key"
  if Rscript scripts/gam/run_logscore_gam.R "$id" "$season" "$year"; then
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
