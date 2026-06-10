#' @description
#' Fit deer movement GAMs (SSFs with penalised smooths) for every deer with a
#' wrangled track file in data/tracks/. This is the GAM analogue of fit_issf.R.
#'
#' Following Klappstein et al. (2024, Methods Ecol Evol), an (integrated) SSF is
#' fit as a stratified Cox proportional-hazards model in mgcv:
#'   gam(cbind(times, stratum) ~ <smooth/linear terms>,
#'       data = ..., family = cox.ph, weights = obs)
#' where `times` is a constant event time, `stratum` groups each observed step
#' with its random points, and `obs` is the used/available (1/0) indicator
#' passed as the cox.ph censoring weight. The reshape is done by
#' prepare_gam_data(); the per-model fit by fit_gam_mod() (both in
#' helper_functions.R).
#'
#' For each (id, season, year):
#'   * If the output rds already exists and `overwrite` is FALSE, skip.
#'   * Otherwise read the wrangled data, reshape the chosen design's covariate
#'     table (`gam_input`) into mgcv cox.ph form, fit all candidate formulas
#'     with fit_gam_mod() (which returns failure-status objects rather than
#'     erroring per-model), and save the list to results/results_gam_<key>.rds.
#' Errors are caught per deer so a single failure does not halt the loop.
#'
#' Random points: two designs are wrangled per deer.
#'   * "stp.var.nonp" — random steps sampled uniformly over a disc of radius
#'     ~ max step length (wrangle_deer_mvt.R, model = "nonp"). Required for
#'     *non-parametric* movement smooths s(sl_) / s(ta_) (Klappstein et al.
#'     2024, Section 3.1); also valid for parametric kernels.
#'   * "stp.var" — gamma / von Mises (iSSF-style). Usable for GAMs with a
#'     *parametric* movement kernel, and more efficient there.
#' Pick via `gam_input` below.
#'
#' Configuration: edit `overwrite` / `gam_input` below before running.

# Configuration ---------------------------------------------------------------
overwrite <- TRUE # set to TRUE to refit deer that already have output

# Which per-deer design to fit:
#   "stp.var"      — gamma / von Mises (parametric movement kernel; default)
#   "stp.var.nonp" — uniform-disc (non-parametric movement kernel)
gam_input <- "stp.var"

# Load packages ---------------------------------------------------------------
library(mgcv)
library(tidyverse)

# helper functions
source("scripts/helper_functions.R")

# Formulas --------------------------------------------------------------------
# RHS only — fit_gam_mod() prepends the "cbind(times, stratum) ~" Cox-PH
# response. Each model is the GAM translation of the like-numbered iSSF formula
# in fit_issf.R (iSSF #5, the three-way interaction, is intentionally omitted).
#
# MOVE = parametric gamma / von Mises movement kernel with a zebra-style cyclic
# spline letting the gamma scale (the step-length coefficient) vary smoothly
# with time of day (Klappstein et al. 2024). The gamma design (stp.var) supplies
# the tentative rate, so log(sl_) + cos(ta_) suffice. Landcover interactions use
# hierarchical "fs" smooths (one shrunk curve per class) rather than independent
# by= smooths, which are unidentifiable when a deer visits only a few of the
# 10 landcover classes.
move <- "log(sl_) + cos(ta_) + s(tod_, bs = 'cc', k = 5, by = sl_)"

formulas <- c(
  # 1  movement only
  move,
  # 2  + distance to HR edge
  paste(move, "+ s(HR_edge_end)"),
  # 3  + distance to HR center
  paste(move, "+ s(HR_center_end)"),
  # 4  movement varies with distance to HR center at the step start
  paste(
    move,
    "+ s(HR_center_start, by = sl_) + s(HR_center_start, by = cos(ta_))"
  ),
  # 6  + HR center + NDVI-by-landcover (hierarchical)
  paste(
    move,
    "+ s(HR_center_end) + s(ndvi_end, wiscland_end, bs = 'fs', k = 5)"
  ),
  # 7  + NDVI-by-landcover (hierarchical), no HR
  paste(move, "+ s(ndvi_end, wiscland_end, bs = 'fs', k = 5)"),
  # 8  + HR center + landcover
  paste(move, "+ s(HR_center_end) + wiscland_end"),
  # 9  + landcover
  paste(move, "+ wiscland_end"),
  # 10 + HR-center-by-landcover (hierarchical)
  paste(move, "+ s(HR_center_end, wiscland_end, bs = 'fs', k = 5)")
)

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
  out_path <- sprintf("results/results_gam_%s.rds", key)

  # Skip: output already exists and we're not overwriting
  if (!overwrite && file.exists(out_path)) {
    cat(sprintf("[skip-exists]   %s\n", key))
    n_skipped <- n_skipped + 1L
    next
  }

  # Process — wrapped so per-deer errors don't halt the loop. fit_gam_mod
  # itself already returns failure-status objects per model, so the global
  # tryCatch is here for unforeseen errors only (data load, missing columns).
  cat(sprintf("[run]           %s\n", key))
  ok <- tryCatch(
    {
      one_deer <- readRDS(in_path)
      gam_data <- prepare_gam_data(one_deer[[gam_input]][[1]])

      results_gam <- purrr::map(seq_along(formulas), function(j) {
        cat(sprintf("  Formula %d\n", j))
        fit_gam_mod(gam_data, formulas[j])
      })

      saveRDS(results_gam, out_path)
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
