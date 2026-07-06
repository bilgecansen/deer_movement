#!/usr/bin/env bash
# Iterate every deer with test-year simulations in sims/issf/sims_issf_test_*.rds and
# run scripts/issf/run_udoverlap_issf_test.R on each. Driving off the test sims (rather
# than fitted models) means deer with no test-year data are skipped naturally
# — they have no sims_test_*.rds file. By default, skip deer whose test
# UD-overlap output already exists; pass --overwrite to reprocess everything.
#
# Failures in Rscript do not halt the loop. Final summary prints counts and
# elapsed time.
#
# Usage:
#   bash scripts/issf/run_udoverlap_all_test.sh              # resumable (default)
#   bash scripts/issf/run_udoverlap_all_test.sh --overwrite  # reprocess all

shopt -s nullglob

overwrite=false
if [[ "${1:-}" == "--overwrite" ]]; then
  overwrite=true
fi

n_done=0
n_skipped=0
n_failed=0
start=$(date +%s)

for f in sims/issf/sims_issf_test_*.rds; do
  base=$(basename "$f" .rds)
  key=${base#sims_issf_test_}
  IFS='_' read -r id season year <<< "$key"

  out_path="filters/issf/udoverlap_issf_test_${key}.rds"

  if [[ "$overwrite" == false && -f "$out_path" ]]; then
    echo "[skip] $key"
    n_skipped=$((n_skipped + 1))
    continue
  fi

  echo "[run]  $key"
  if Rscript scripts/issf/run_udoverlap_issf_test.R "$id" "$season" "$year"; then
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
