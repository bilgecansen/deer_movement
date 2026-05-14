#!/usr/bin/env bash
# Iterate every deer with a fitted iSSF model in results/ and run
# scripts/run_logscore_test.R on each, scoring the test-year (year+1) track
# under the train-year model. Deer with no test-year wrangled track are
# counted as "no data" and skipped without invoking R.
#
# By default, deer whose test log-score output already exists are skipped;
# pass --overwrite to reprocess everything. Failures in Rscript do not halt
# the loop. Final summary prints counts and elapsed time.
#
# Usage:
#   bash scripts/run_logscore_all_test.sh              # resumable (default)
#   bash scripts/run_logscore_all_test.sh --overwrite  # reprocess all

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

for f in results/results_issf_*.rds; do
  base=$(basename "$f" .rds)
  key=${base#results_issf_}
  IFS='_' read -r id season year <<< "$key"

  test_year=$((year + 1))
  test_key="${id}_${season}_${test_year}"
  test_track="data/tracks/data_${test_key}.rds"
  out_path="filters/logscore_test_${key}.rds"

  if [[ ! -f "$test_track" ]]; then
    echo "[no-data] $key (no $test_track)"
    n_no_data=$((n_no_data + 1))
    continue
  fi

  if [[ "$overwrite" == false && -f "$out_path" ]]; then
    echo "[skip]    $key"
    n_skipped=$((n_skipped + 1))
    continue
  fi

  echo "[run]     $key"
  if Rscript scripts/run_logscore_test.R "$id" "$season" "$year"; then
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
