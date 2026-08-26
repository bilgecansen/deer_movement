#!/usr/bin/env bash
# Run every correctness check. Exits non-zero if any check failed, so this can
# gate a batch run or a commit.
#
# Usage: bash scripts/checks/run_checks_all.sh [<id> <season> <year>]
#   The optional deer key is passed to the per-deer kernel checks; without it
#   they use their built-in default deer.

set -u
rc=0

echo "### Tier 6: cross-output audit (gam) ###"
Rscript scripts/checks/check_outputs.R gam || rc=1
Rscript scripts/checks/check_outputs.R amt || rc=1

echo
echo "### Tier 3: analytic reductions + contracts (gam) ###"
Rscript scripts/checks/check_kernel_gam.R "$@" || rc=1

echo
echo "### Tier 2: amt vs GAM on real deer ###"
Rscript scripts/checks/check_differential_gam.R 5 30 || rc=1

echo
if [[ $rc -eq 0 ]]; then
  echo "ALL CHECKS PASSED"
else
  echo "SOME CHECKS FAILED (exit $rc)"
fi
exit $rc
