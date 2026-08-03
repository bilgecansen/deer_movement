#!/usr/bin/env bash
# Iterate every deer with a fitted GAM model in results/ and run
# scripts/gam/run_logscore_gam_test.R on each, scoring the test-year (year+1) track
# under the train-year model. Deer with no test-year wrangled track are
# counted as "no data" and skipped without invoking R.
#
# By default, deer whose test log-score output already exists are skipped;
# pass --overwrite to reprocess everything. Failures in Rscript do not halt
# the loop. Final summary prints counts and elapsed time.
#
# Usage:
#   bash scripts/gam/run_logscore_gam_all_test.sh              # resumable (default)
#   bash scripts/gam/run_logscore_gam_all_test.sh --overwrite  # reprocess all

shopt -s nullglob

overwrite=false
if [[ "${1:-}" == "--overwrite" ]]; then
  overwrite=true
fi

n_done=0
n_skipped=0
n_no_data=0
n_failed=0
start=$(date +%s)

for f in results/gam/results_gam_*.rds; do
  base=$(basename "$f" .rds)

  # results_gam_null_<key>.rds also matches this glob; the null is scored inside
  # run_logscore_gam_test.R, not as a deer of its own. Skip it silently.
  [[ "$base" == results_gam_null_* ]] && continue

  key=${base#results_gam_}
  IFS='_' read -r id season year <<< "$key"

  test_year=$((year + 1))
  test_key="${id}_${season}_${test_year}"
  test_track="data/tracks/data_${test_key}.rds"

  # One invocation writes both files; only skip when both are already present.
  out_path="filters/gam/logscore_gam_test_${key}.rds"
  null_out="filters/gam/logscore_gam_null_test_${key}.rds"

  if [[ ! -f "$test_track" ]]; then
    echo "[no-data] $key (no $test_track)"
    n_no_data=$((n_no_data + 1))
    continue
  fi

  if [[ "$overwrite" == false && -f "$out_path" && -f "$null_out" ]]; then
    echo "[skip]    $key"
    n_skipped=$((n_skipped + 1))
    continue
  fi

  echo "[run]     $key"
  if Rscript scripts/gam/run_logscore_gam_test.R "$id" "$season" "$year"; then
    n_done=$((n_done + 1))
  else
    rc=$?
    echo "[fail]    $key (exit $rc)"
    n_failed=$((n_failed + 1))
  fi
done

elapsed=$(( $(date +%s) - start ))
mins=$(( elapsed / 60 ))
secs=$(( elapsed % 60 ))
echo
echo "Done: $n_done   Skipped: $n_skipped   No data: $n_no_data   Failed: $n_failed   Elapsed: ${mins}m ${secs}s"
