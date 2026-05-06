#' @description
#' Fit deer iSSF models for every deer with a wrangled track file in
#' data/tracks/. For each (id, season, year):
#'   * If the output rds already exists and `overwrite` is FALSE, skip.
#'   * Otherwise read the wrangled data, fit all candidate formulas with
#'     fit_mod() (which itself returns failure-status objects rather than
#'     erroring out per-model), and save the list of results.
#' Errors are caught per deer so a single failure does not halt the loop.
#'
#' Configuration: edit `overwrite` below before running.

# Configuration ---------------------------------------------------------------
overwrite <- TRUE # set to TRUE to refit deer that already have output

# Load packages ---------------------------------------------------------------
library(amt)
library(terra)
library(tidyverse)

# helper functions
source("scripts/helper_functions.R")

# Formulas --------------------------------------------------------------------
formulas <- c(
  "case_ ~ (log(sl_) + cos(ta_)):tod_start_",

  "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + HR_edge_end",

  "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + HR_center_log_end",

  "case_ ~ (log(sl_) + cos(ta_)):tod_start_ +
    (log(sl_) + cos(ta_)):HR_center_log_start",

  "case_ ~ (log(sl_) + cos(ta_)):tod_start_:HR_center_log_start",

  "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + HR_center_log_end +
    east_end + east2_end + north_end + dist_end + wiscland_end + ndvi_end +
    wiscland_end:ndvi_end",

  "case_ ~ log(sl_) + cos(ta_) + east_end + east2_end +
    north_end + dist_end + wiscland_end + ndvi_end + wiscland_end:ndvi_end",

  "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + HR_center_log_end +
    east_end + east2_end + north_end + dist_end + wiscland_end",

  "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + east_end + east2_end +
    north_end + dist_end + wiscland_end",

  "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + HR_center_log_end + dist_end +
    wiscland_end + HR_center_log_end:dist_end + HR_center_log_end:wiscland_end"
)

formulas <- paste(formulas, "+ strata(step_id_)")

# Discover deer from data/tracks/ ---------------------------------------------
track_files <- list.files(
  "data/tracks",
  pattern = "^data_.*\\.rds$",
  full.names = TRUE
)

keys <- gsub("^data_(.*)\\.rds$", "\\1", basename(track_files))

cat(sprintf("Found %d wrangled deer in data/tracks/\n", length(track_files)))

# Output dir
dir.create("results", showWarnings = FALSE)

# Loop ------------------------------------------------------------------------
n_done <- 0L
n_skipped <- 0L
n_failed <- 0L

start_time <- Sys.time()

for (i in seq_along(track_files)) {
  key <- keys[i]
  in_path <- track_files[i]
  out_path <- sprintf("results/results_issf_%s.rds", key)

  # Skip: output already exists and we're not overwriting
  if (!overwrite && file.exists(out_path)) {
    cat(sprintf("[skip-exists]   %s\n", key))
    n_skipped <- n_skipped + 1L
    next
  }

  # Process — wrapped so per-deer errors don't halt the loop. fit_mod itself
  # already returns failure-status objects per model, so the global tryCatch
  # is here for unforeseen errors only (data load issues, missing columns,
  # etc.).
  cat(sprintf("[run]           %s\n", key))
  ok <- tryCatch(
    {
      one_deer <- readRDS(in_path)
      ssf_data <- one_deer$stp.var[[1]]

      results_issf <- purrr::map(seq_along(formulas), function(j) {
        cat(sprintf("  Formula %d\n", j))
        fit_mod(ssf_data, formulas[j])
      })

      saveRDS(results_issf, out_path)
      TRUE
    },
    error = function(e) {
      cat(sprintf("[fail]          %s: %s\n", key, conditionMessage(e)))
      FALSE
    }
  )

  if (ok) {
    n_done <- n_done + 1L
  } else {
    n_failed <- n_failed + 1L
  }
}

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf(
  "Done: %d   Skipped: %d   Failed: %d   Elapsed: %.1f min\n",
  n_done,
  n_skipped,
  n_failed,
  elapsed
))
