# Retired scripts

Out-of-sample (`_test`) variants of the amt workflow, parked here while the
out-of-sample approach is reconsidered. They are the versions that existed
before the amt path was restructured to four numbered models scored against the
GAM null, so they still assume the old 1..6 model set and the old log-score
contract.

Do not run them as-is. If the `_test` idea is revived, port them the same way
the in-sample scripts were:

  * four numbered models (see scripts/amt/fit_amt.R)
  * the incoming heading passed to make_start(), first-in-burst steps skipped
  * totals over the step set common to every model, with status counts
  * the GAM null as the reference, scored in the same run

The GAM path has no `_test` variants at all -- they were never written.

## plot_interactions_xgb.R

The interaction plot of the selected xgb model types without deer: the
held-out log score of each pair of variables (xgb_heldout_split_gain()).
Pair scores are only valid at max_depth 2, so the deer types (depth 3)
never had them, and the xgb workflow now fits deer types only. The engine
still computes split gain for the types without deer, which the simulation
test (scripts/xgb_sim/) reads directly. The script expects both families
in results/xgb/ and does not run as-is.
