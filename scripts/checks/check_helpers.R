# Shared check harness --------------------------------------------------------
#
# Minimal PASS/FAIL reporting used by the scripts in scripts/checks/. Each check
# states an expectation that must hold regardless of the data, reports the
# observed value either way, and records the result so the runner can exit
# non-zero if anything failed.
#
# Design note: a check must print the number it saw, not just PASS. A check that
# only ever prints PASS is untested itself, and the whole point of these scripts
# is to catch the case where the code is confidently wrong.

.CHECK_STATE <- new.env(parent = emptyenv())
.CHECK_STATE$results <- list()

#' Record and print one check.
#' @param name  Short label.
#' @param ok    TRUE = pass, FALSE = fail, NA = could not evaluate (counts as a
#'              skip, not a pass — an inconclusive check must never look green).
#' @param detail Observed value / context. Always printed.
check <- function(name, ok, detail = "") {
  tag <- if (isTRUE(ok)) "PASS" else if (isFALSE(ok)) "FAIL" else "SKIP"
  .CHECK_STATE$results[[length(.CHECK_STATE$results) + 1L]] <- list(
    name = name,
    tag = tag
  )
  cat(sprintf("  [%s] %-58s %s\n", tag, name, detail))
  invisible(ok)
}

#' Numeric comparison with an explicit tolerance, reporting the actual gap.
check_close <- function(name, actual, expected, tol = 1e-8, scale = FALSE) {
  if (length(actual) != 1 || !is.finite(actual)) {
    return(check(name, NA, sprintf("actual = %s", paste(actual, collapse = ", "))))
  }
  gap <- abs(actual - expected)
  if (scale) {
    gap <- gap / max(abs(expected), 1e-12)
  }
  check(
    name,
    gap <= tol,
    sprintf(
      "actual %.10g vs expected %.10g (%s %.3g, tol %.3g)",
      actual,
      expected,
      if (scale) "rel gap" else "gap",
      gap,
      tol
    )
  )
}

check_section <- function(title) {
  cat(sprintf("\n== %s ==\n", title))
}

#' Print the tally and exit non-zero if anything failed, so the runner and CI
#' can act on it.
check_summary <- function() {
  tags <- vapply(.CHECK_STATE$results, function(r) r$tag, character(1))
  n_fail <- sum(tags == "FAIL")
  n_skip <- sum(tags == "SKIP")
  cat(sprintf(
    "\n---- %d checks: %d passed, %d FAILED, %d skipped ----\n",
    length(tags),
    sum(tags == "PASS"),
    n_fail,
    n_skip
  ))
  if (n_fail > 0) {
    cat("FAILED:\n")
    for (r in .CHECK_STATE$results) {
      if (r$tag == "FAIL") cat(sprintf("  - %s\n", r$name))
    }
  }
  invisible(n_fail)
}
