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
# make_formulas() returns the 6 candidate models for one deer, named by model
# number ("1","2","3","4","5","6"). Slots 4 & 5 are the resource-selection
# models and depend on season: non-breeding (winter) deer have no valid MODIS
# NDVI (snow / dormancy -> all-or-mostly NA), so the NDVI models (4 & 5) can't
# be fit. For those deer we substitute the landcover-only models and record them
# in slots 4 & 5. Either way the result holds 6 models and slots 4/5 always mean
# "the resource-selection model" (NDVI where it exists, a landcover-only
# substitute where it doesn't). Model 6 is the HR-center × landcover
# interaction. Numbering is contiguous "1".."6" and parallels make_formulas() in
# fit_GAM.R (position = model number).
make_formulas <- function(season) {
  move <- "case_ ~ (sl_):tod_start_ + log(sl_) + cos(ta_)"

  f1 <- move # 1 movement only
  f2 <- paste(move, "+ HR_edge_end") # 2 + HR edge
  f3 <- paste(move, "+ HR_center_end") # 3 + HR center
  f6 <- paste(
    move,
    "+ HR_center_end + wiscland_end + HR_center_end:wiscland_end"
  ) # 6 + HR-center x LC

  if (season == "nb") {
    # 4 landcover (with HR-center)
    f4 <- paste(move, "+ HR_center_end + wiscland_end")
    f5 <- paste(move, "+ wiscland_end") # 5 landcover (no HR-center)
  } else {
    f4 <- paste(
      move,
      "+ HR_center_end + wiscland_end + ndvi_end + wiscland_end:ndvi_end"
    ) # 4 NDVI (with HR-center)
    f5 <- paste(
      move,
      "+ wiscland_end + ndvi_end + wiscland_end:ndvi_end"
    ) # 5 NDVI (no HR-center)
  }

  f <- paste(c(f1, f2, f3, f4, f5, f6), "+ strata(step_id_)")
  stats::setNames(f, c("1", "2", "3", "4", "5", "6"))
}

# Discover deer from data/tracks/ ---------------------------------------------
track_files <- list.files(
  "data/tracks",
  pattern = "^data_.*\\.rds$",
  full.names = TRUE
)

keys <- gsub("^data_(.*)\\.rds$", "\\1", basename(track_files))

cat(sprintf("Found %d wrangled deer in data/tracks/\n", length(track_files)))

# Output dir
dir.create("results/issf", showWarnings = FALSE, recursive = TRUE)

# Loop ------------------------------------------------------------------------
n_done <- 0L
n_skipped <- 0L
n_failed <- 0L

start_time <- Sys.time()

for (i in seq_along(track_files)) {
  key <- keys[i]
  in_path <- track_files[i]
  out_path <- sprintf("results/issf/results_issf_%s.rds", key)

  # Skip: output already exists and we're not overwriting
  if (!overwrite && file.exists(out_path)) {
    cat(sprintf("[skip-exists]   %s\n", key))
    n_skipped <- n_skipped + 1L
    next
  }

  # Process — wrapped so per-deer errors don't halt the loop. fit_mod itself
  # already returns failure-status objects per model, so the global tryCatch
  # is here for unforeseen errors only (data load issues, missing
  # columns, etc.).
  cat(sprintf("[run]           %s\n", key))
  ok <- tryCatch(
    {
      one_deer <- readRDS(in_path)
      ssf_data <- one_deer$stp.var[[1]]

      # Season picks the slot-4/5 resource models (NDVI vs landcover; see
      # make_formulas). results_issf is named by model number ("1".."6").
      season <- strsplit(key, "_")[[1]][2]
      formulas <- make_formulas(season)

      results_issf <- purrr::map(formulas, function(f) fit_mod(ssf_data, f))

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
