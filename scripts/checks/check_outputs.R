#' @description
#' Cross-output audit: scan every per-deer filter output for values that should
#' not be possible, or patterns that indicate a silent bug rather than a real
#' result. Reads only what is already on disk — no refitting, no simulation — so
#' it runs in seconds and can be re-run after every batch.
#'
#' This is the cheapest tier of correctness checking. It cannot tell you the
#' model is right; it tells you when the pipeline produced something the maths
#' says it shouldn't have.
#'
#' Checks, per path (issf / gam):
#'   * bat_uds and svf_agree inside [0, 1] (both are overlap/agreement measures)
#'   * neither pinned at exactly 0 or 1 (a degenerate fit, not a real overlap)
#'   * n_steps identical across models within a deer — delta_logp differences
#'     computed over different step sets are not comparable
#'   * delta_logp not numerically zero for a model that adds terms to the null
#'   * energy_score strictly positive and finite
#'   * no metric value repeated across an implausible number of deer (a caching
#'     or copy/paste bug looks exactly like this)
#'   * the same model set present in every deer's file
#'
#' Usage: Rscript scripts/checks/check_outputs.R [issf|gam|both]

suppressPackageStartupMessages(library(tidyverse))
source("scripts/checks/check_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
paths <- if (length(args) && args[1] != "both") args[1] else c("issf", "gam")

# Tunables --------------------------------------------------------------------
# DELTA_ZERO_TOL: a model that adds terms should change the log likelihood by
#   more than round-off. Anything under this is either a term shrunk fully to
#   zero (legitimate, but worth seeing) or a term that never entered the fit.
# DUP_FRAC: fraction of deer sharing one identical metric value that we treat as
#   suspicious. Continuous metrics should essentially never tie.
DELTA_ZERO_TOL <- 1e-3
DUP_FRAC <- 0.02

read_udoverlap <- function(f) {
  ud <- readRDS(f)
  purrr::imap_dfr(ud, function(x, m) {
    if (length(x) == 1 && is.na(x)) {
      return(tibble(model = m, bat_uds = NA_real_, svf_agree = NA_real_))
    }
    # Accept the pre-rename element name (see filter_models_*.R).
    tibble(
      model = m,
      bat_uds = x$bat_uds,
      svf_agree = if (is.null(x[["svf_agree"]])) {
        x[["svf_score"]]
      } else {
        x[["svf_agree"]]
      }
    )
  }) |>
    mutate(key = gsub("^udoverlap_[a-z]+_(.*)\\.rds$", "\\1", basename(f)))
}

for (path in paths) {
  check_section(sprintf("%s outputs", toupper(path)))
  dir <- file.path("filters", path)
  if (!dir.exists(dir)) {
    check(sprintf("%s: filters dir exists", path), NA, dir)
    next
  }

  # ---- UD overlap + SVF -----------------------------------------------------
  udf <- list.files(
    dir, sprintf("^udoverlap_%s_[^t].*\\.rds$", path), full.names = TRUE
  )
  udf <- udf[!grepl("_test_", basename(udf))]
  if (length(udf)) {
    ud <- purrr::map_dfr(udf, read_udoverlap)
    cat(sprintf("  (%d deer, %d deer-model rows)\n", length(udf), nrow(ud)))

    for (v in c("bat_uds", "svf_agree")) {
      x <- ud[[v]]
      xf <- x[is.finite(x)]
      bad <- sum(xf < 0 | xf > 1)
      check(
        sprintf("%s within [0, 1]", v),
        bad == 0,
        sprintf(
          "%d/%d outside; range [%.4f, %.4f]; NA %d",
          bad, length(xf), min(xf), max(xf), sum(is.na(x))
        )
      )
      n_pin <- sum(xf == 0 | xf == 1)
      check(
        sprintf("%s never pinned at exactly 0 or 1", v),
        n_pin == 0,
        sprintf("%d pinned", n_pin)
      )
      tt <- table(xf)
      worst <- if (length(tt)) max(tt) else 0L
      check(
        sprintf("%s has no over-repeated value", v),
        worst <= max(2, DUP_FRAC * length(xf)),
        sprintf("most repeated value occurs %d times (of %d)",
                worst, length(xf))
      )
    }
  } else {
    check(sprintf("%s: udoverlap files found", path), NA, "none")
  }

  # ---- Log score ------------------------------------------------------------
  lsf <- list.files(
    dir, sprintf("^logscore_%s_.*\\.rds$", path), full.names = TRUE
  )
  lsf <- lsf[!grepl("_test_|_null_", basename(lsf))]
  if (length(lsf)) {
    ls_df <- purrr::map_dfr(lsf, function(f) {
      readRDS(f) |>
        as_tibble() |>
        mutate(key = gsub("^logscore_[a-z]+_(.*)\\.rds$", "\\1", basename(f)))
    })
    cat(sprintf("  (%d deer, %d deer-model rows)\n", length(lsf), nrow(ls_df)))

    # Every model in a deer must be scored on the same steps, or delta_logp is
    # differencing likelihoods computed over different data.
    per_deer <- ls_df |>
      group_by(key) |>
      summarise(n_distinct_steps = n_distinct(n_steps[n_steps > 0]),
                .groups = "drop")
    bad <- sum(per_deer$n_distinct_steps > 1)
    check(
      "n_steps identical across models within a deer",
      bad == 0,
      sprintf("%d/%d deer disagree", bad, nrow(per_deer))
    )

    check(
      "total_logp finite and negative",
      all(ls_df$total_logp[!is.na(ls_df$total_logp)] < 0),
      sprintf(
        "range [%.1f, %.1f]; NA %d",
        min(ls_df$total_logp, na.rm = TRUE),
        max(ls_df$total_logp, na.rm = TRUE),
        sum(is.na(ls_df$total_logp))
      )
    )

    if ("delta_logp" %in% names(ls_df)) {
      # A model that adds terms should move the likelihood by more than
      # round-off. Exclude the reference itself (delta exactly 0 by definition).
      d <- ls_df$delta_logp
      d <- d[is.finite(d) & d != 0]
      near0 <- sum(abs(d) < DELTA_ZERO_TOL)
      check(
        "delta_logp not numerically zero",
        near0 == 0,
        sprintf("%d/%d rows |delta| < %.0e", near0, length(d), DELTA_ZERO_TOL)
      )
    }

    # Model set consistency
    sets <- ls_df |>
      group_by(key) |>
      summarise(s = paste(sort(unique(model)), collapse = ","),
                .groups = "drop")
    check(
      "same model set in every deer",
      n_distinct(sets$s) == 1,
      sprintf("%d distinct sets: %s", n_distinct(sets$s),
              paste(head(unique(sets$s), 3), collapse = " | "))
    )
  } else {
    check(sprintf("%s: logscore files found", path), NA, "none")
  }

  # ---- Energy score ---------------------------------------------------------
  esf <- file.path(dir, sprintf("es_%s.rds", path))
  if (file.exists(esf)) {
    es <- readRDS(esf)
    e <- es$energy_score
    ef <- e[is.finite(e)]
    check(
      "energy_score strictly positive",
      all(ef > 0),
      sprintf("%d/%d <= 0; range [%.1f, %.1f]; NA %d",
              sum(ef <= 0), length(ef), min(ef), max(ef), sum(is.na(e)))
    )
    tt <- table(ef)
    worst <- if (length(tt)) max(tt) else 0L
    check(
      "energy_score has no over-repeated value",
      worst <= max(2, DUP_FRAC * length(ef)),
      sprintf("most repeated value occurs %d times (of %d)", worst, length(ef))
    )
  } else {
    check(sprintf("%s: es file found", path), NA, esf)
  }
}

quit(status = if (check_summary() > 0) 1L else 0L)
