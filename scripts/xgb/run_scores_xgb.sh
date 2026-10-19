#!/usr/bin/env bash
# Run score_sims_xgb.R batch by batch, BATCH deer-years at a time, each batch
# in a fresh R, with a pause between them. ctmm fits (HR overlap) now and
# then hang and never return, which would hold a single long run forever.
# Here a batch that writes no score file for STALL_MINS minutes is stopped,
# workers and all, and the next batch takes up what is left, the stopped
# deer-years included. Stops when every deer-year is scored, or after
# MAX_IDLE batches in a row that score nothing; the last run of the scorer
# writes the table of all of them.
#
# Usage, from the repository root:
#   bash scripts/xgb/run_scores_xgb.sh SET MODELS [BATCH]
#   e.g. bash scripts/xgb/run_scores_xgb.sh train gam 50
# Progress: logs/xgb/scores_<MODELS>_<SET>_gamma99.log

set -u
# Job control: every batch gets its own process group, so stopping it stops
# its workers too
set -m

SET=$1
MODELS=$2
BATCH=${3:-50}
STALL_MINS=10
MAX_IDLE=3
PAUSE_SECONDS=5

if [ "$SET" = "test" ]; then
  sims_dir=sims/xgb_paths_gamma99
else
  sims_dir=sims/xgb_paths_train_gamma99
fi
if [ "$MODELS" = "xgb" ]; then
  out_dir=filters/xgb/scores_${SET}_gamma99
else
  out_dir=filters/xgb/scores_gam_${SET}_gamma99
fi
log=logs/xgb/scores_${MODELS}_${SET}_gamma99.log

n_scored() {
  ls "$out_dir" 2>/dev/null | grep -c "^scores_${MODELS}_.*\.rds$"
}
n_left() {
  comm -23 \
    <(ls "$sims_dir" | sed -n "s/^sims_${MODELS}_\(.*\)\.rds$/\1/p" | sort) \
    <(ls "$out_dir" 2>/dev/null |
        sed -n "s/^scores_${MODELS}_\(.*\)\.rds$/\1/p" | sort) |
    wc -l | tr -d ' '
}
note() {
  echo "$(date '+%m-%d %H:%M') $*" >> "$log"
}
# Running, and not a finished process waiting to be reaped
alive() {
  local state
  state=$(ps -o stat= -p "$1" 2>/dev/null)
  [ -n "$state" ] && [ "${state#Z}" = "$state" ]
}

: > "$log"
note "start: $(n_left) $SET $MODELS deer-years to score, $BATCH a batch"
idle=0
while [ "$(n_left)" -gt 0 ]; do
  before=$(n_scored)
  Rscript scripts/xgb/score_sims_xgb.R "$SET" "$MODELS" "$BATCH" \
    >> "$log" 2>&1 &
  pid=$!
  count=$before
  last=$(date +%s)
  while alive "$pid"; do
    sleep 20
    now=$(n_scored)
    if [ "$now" -gt "$count" ]; then
      count=$now
      last=$(date +%s)
    elif [ $(( $(date +%s) - last )) -ge $(( STALL_MINS * 60 )) ]; then
      note "no new score in $STALL_MINS min: stopping the batch"
      kill -TERM -- -"$pid" 2>/dev/null
      sleep 5
      kill -KILL -- -"$pid" 2>/dev/null
      break
    fi
  done
  wait "$pid" 2>/dev/null
  done_now=$(( $(n_scored) - before ))
  note "batch done: $done_now scored, $(n_left) left"
  if [ "$done_now" -eq 0 ]; then
    idle=$(( idle + 1 ))
    if [ "$idle" -ge "$MAX_IDLE" ]; then
      note "$MAX_IDLE batches in a row scored nothing: giving up"
      exit 1
    fi
  else
    idle=0
  fi
  sleep "$PAUSE_SECONDS"
done
# Nothing left to score: this run only writes the table
Rscript scripts/xgb/score_sims_xgb.R "$SET" "$MODELS" >> "$log" 2>&1
note "finished"
