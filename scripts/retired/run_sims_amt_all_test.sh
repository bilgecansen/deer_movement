#!/usr/bin/env bash
# Iterate every deer with a fitted amt model in results/ and run
# scripts/amt/run_sims_amt_test.R on each — i.e. simulate year+1 paths
# using the
# year-trained model. By default, skip deer whose test-simulation output
# already exists; pass --overwrite to reprocess everything.
#
# run_sims_test.R exits with status 2 when no test-year wrangled track exists
# for that deer; those are counted separately as "no data". All other non-zero
# exits are counted as failures. Failures do not halt the loop. Final summary
# prints counts and elapsed time.
#
# Usage:
#   bash scripts/amt/run_sims_all_test.sh              # resumable (default)
#   bash scripts/amt/run_sims_all_test.sh --overwrite  # reprocess all

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

for f in results/amt/results_amt_*.rds; do
  base=$(basename "$f" .rds)
  key=${base#results_amt_}
  IFS='_' read -r id season year <<< "$key"

  out_path="sims/amt/sims_amt_test_${key}.rds"

  if [[ "$overwrite" == false && -f "$out_path" ]]; then
    echo "[skip]    $key"
    n_skipped=$((n_skipped + 1))
    continue
  fi

  echo "[run]     $key"
  Rscript scripts/amt/run_sims_amt_test.R "$id" "$season" "$year"
  rc=$?
  case $rc in
    0)
      n_done=$((n_done + 1))
      ;;
    2)
      echo "[no-data] $key"
      n_no_data=$((n_no_data + 1))
      ;;
    *)
      echo "[fail]    $key (exit $rc)"
      n_failed=$((n_failed + 1))
      ;;
  esac
done

elapsed=$(( $(date +%s) - start ))
mins=$(( elapsed / 60 ))
secs=$(( elapsed % 60 ))
echo
echo "Done: $n_done   Skipped: $n_skipped   No data: $n_no_data"
echo "Failed: $n_failed   Elapsed: ${mins}m ${secs}s"
