#' @description
#' Fit the amt (conditional-logistic) movement models for every deer with a
#' wrangled track file in data/tracks/. For each (id, season, year):
#'   * If the output rds already exists and `overwrite` is FALSE, skip.
#'   * Otherwise read the wrangled data, fit all candidate formulas with
#'     fit_mod() (which itself returns failure-status objects rather than
#'     erroring out per-model), and save the list of results.
#' Errors are caught per deer so a single failure does not halt the loop.
#'
#' NO NULL MODEL IS FIT HERE. The amt models are compared against the GAM null
#' (movement + s(HR_center_end)) from results/gam/results_gam_null_<key>.rds, so
#' fit_GAM.R must have been run before any amt log score. The consequence worth
#' remembering: the GAM null carries a cyclic time-of-day smooth, a more
#' flexible movement block than the day/night factor these models use, so
#' delta_logp mixes "does habitat help" with "my movement block is less flexible
#' than the reference's".
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
# make_formulas() returns the 4 candidate models for one deer, named by model
# number ("1".."4"). The set mirrors make_formulas() in fit_GAM.R slot for slot,
# with linear / factor terms where the GAM uses penalised smooths:
#
#   #  amt                                        GAM
#   1  movement only                              movement only
#   2  + HR_center + resource                     + s(HR_center) + resource
#   3  + resource                                 + resource
#   4  + HR_center x landcover                    + s(HR_center) + hc_fs
#
# The resource block is season-dependent, exactly as in the GAM path: winter
# (nb) deer have no usable MODIS NDVI (snow / dormancy -> all-or-mostly NA), so
# slots 2 and 3 fall back to landcover alone. Either way slots 2/3 always mean
# "the resource-selection model".
#
# Retired relative to the old 1..6 set: the HR-edge model (old 2), and the
# HR-center-only model (old 3), which is now the GAM null and is not refit here.
# Old 1, 4, 5, 6 become new 1, 2, 3, 4 -- the same remapping the GAM path used.
make_formulas <- function(season) {
  move <- "case_ ~ (sl_):tod_start_ + log(sl_) + cos(ta_)"

  f1 <- move # 1 movement only
  # 4 HR-center x landcover
  f4 <- paste(
    move,
    "+ HR_center_end + wiscland_end + HR_center_end:wiscland_end"
  )

  if (season == "nb") {
    # 2 landcover, with HR-center
    f2 <- paste(move, "+ HR_center_end + wiscland_end")
    # 3 landcover, no HR-center
    f3 <- paste(move, "+ wiscland_end")
  } else {
    # 2 NDVI x landcover, with HR-center
    f2 <- paste(
      move,
      "+ HR_center_end + wiscland_end + ndvi_end + wiscland_end:ndvi_end"
    )
    # 3 NDVI x landcover, no HR-center
    f3 <- paste(
      move,
      "+ wiscland_end + ndvi_end + wiscland_end:ndvi_end"
    )
  }

  f <- paste(c(f1, f2, f3, f4), "+ strata(step_id_)")
  stats::setNames(f, c("1", "2", "3", "4"))
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
dir.create("results/amt", showWarnings = FALSE, recursive = TRUE)

# Loop ------------------------------------------------------------------------
n_done <- 0L
n_skipped <- 0L
n_failed <- 0L

start_time <- Sys.time()

for (i in seq_along(track_files)) {
  key <- keys[i]
  in_path <- track_files[i]
  out_path <- sprintf("results/amt/results_amt_%s.rds", key)

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
      # make_formulas). results_amt is named by model number ("1".."6").
      season <- strsplit(key, "_")[[1]][2]
      formulas <- make_formulas(season)

      results_amt <- purrr::map(formulas, function(f) fit_mod(ssf_data, f))

      saveRDS(results_amt, out_path)
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
