#' @description
#' Deer movement data wrangling — full track per deer.
#' Iterates sequentially over every deer on the canonical list
#' (canonical_deer(), the same list create_hr_rasters.R uses) and writes one
#' rds per deer to data/tracks/data_<id>_<season>_<year>.rds. Each file
#' carries the deer's `set`, "train" or "test".
#'
#' Behavior per deer:
#'   * If the HR rasters are missing for that (id, season, year), skip with a
#'     warning. extract_step_variables can't run without them.
#'   * If the output rds already exists and `overwrite` is FALSE, skip.
#'   * Otherwise run make_random_pt_extraction -> extract_step_variables for
#'     each design and save the result. The raster window is the deer's own
#'     (deer_crop_buffer), wide enough for its longest step. Errors are
#'     caught per deer so a single failure does not halt the loop.
#'
#' Configuration: edit `overwrite` below before running.

# Configuration ---------------------------------------------------------------
overwrite <- TRUE # set to TRUE to reprocess deer that already have output

# Load packages ---------------------------------------------------------------
library(amt)
library(terra)
library(tidyverse)
library(sf)

# helper functions
source("scripts/helper_functions.R")

# Load shared inputs ----------------------------------------------------------
start_time <- Sys.time()

# landscape data — landcover, NDVI and LANDFIRE are annual, so all three are
# loaded per deer inside the loop (see load_landcover / make_water_mask /
# load_ndvi / load_landfire). Elevation, northness and eastness are static,
# so they are loaded once here. No env_old.
topo <- load_topo()

# Deer movement data
sw_deer_tracks <- readRDS('library/SW_filtered_deer.RData')

# The canonical deer list, shared with create_hr_rasters.R so the HR rasters
# line up by deer key
deer_mvt <- canonical_deer(sw_deer_tracks)

cat(sprintf("Total qualifying deer: %d\n", nrow(deer_mvt)))

# Output dir
dir.create("data/tracks", recursive = TRUE, showWarnings = FALSE)

# Loop -------------------------------------------------------------------------
n_done <- 0L
n_skipped <- 0L
n_failed <- 0L

for (i in seq_len(nrow(deer_mvt))) {
  one_deer <- deer_mvt[i, ]
  key <- sprintf("%s_%s_%d", one_deer$id, one_deer$season, one_deer$year)
  out_path <- sprintf("data/tracks/data_%s.rds", key)
  hr_path <- sprintf("data/HR/HRbin_%s.tif", key)

  # Skip: HR raster missing -> can't run extract_step_variables
  if (!file.exists(hr_path)) {
    cat(sprintf("[skip-missing-HR] %s\n", key))
    n_skipped <- n_skipped + 1L
    next
  }

  # Skip: output already exists and we're not overwriting
  if (!overwrite && file.exists(out_path)) {
    cat(sprintf("[skip-exists]      %s\n", key))
    n_skipped <- n_skipped + 1L
    next
  }

  # Process — wrapped so per-deer errors don't halt the loop
  cat(sprintf("[run]              %s\n", key))
  ok <- tryCatch(
    {
      # Season-specific annual landcover (+ open-water mask), NDVI stack and
      # LANDFIRE layers
      landcover <- load_landcover(one_deer$year, one_deer$season)
      water <- make_water_mask(landcover)
      ndvi <- load_ndvi(one_deer$year)
      landfire <- load_landfire(one_deer$year, one_deer$season)

      # Gamma / von Mises design (parametric): used by amt and by parametric
      # GAMs -> stp.random / stp.var. 25 random points (Klappstein et al. 2024).
      one_deer <- make_random_pt_extraction(
        data = one_deer,
        n_pts = 25,
        water = water,
        model = "amt",
        stp_col = "stp",
        output_col = "stp.random"
      )
      one_deer <- extract_step_variables(
        data = one_deer,
        env = landcover,
        ndvi = ndvi,
        landfire = landfire,
        topo = topo,
        random_col = "stp.random",
        output_col = "stp.var"
      )

      # Non-parametric design, the one the xgb path reads: 100 random points
      # uniform in distance out to the deer's longest step with a turning
      # angle, every one counting the same (see nonp_random_points) ->
      # stp.random.nonp / stp.var.nonp
      one_deer <- make_random_pt_extraction(
        data = one_deer,
        n_pts = 100,
        water = water,
        model = "nonp",
        stp_col = "stp",
        output_col = "stp.random.nonp"
      )
      one_deer <- extract_step_variables(
        data = one_deer,
        env = landcover,
        ndvi = ndvi,
        landfire = landfire,
        topo = topo,
        random_col = "stp.random.nonp",
        output_col = "stp.var.nonp"
      )
      saveRDS(one_deer, out_path)
      TRUE
    },
    error = function(e) {
      cat(sprintf("[fail]             %s: %s\n", key, conditionMessage(e)))
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
