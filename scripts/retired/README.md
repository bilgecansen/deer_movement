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
