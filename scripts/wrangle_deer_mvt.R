#' @description
#' Deer movement data wrangling — full track per deer (no train/test split).
#' Iterates sequentially over every deer matching the create_hr_rasters.R filter
#' and writes one rds per deer to data/tracks/data_<id>_<season>_<year>.rds.
#'
#' Behavior per deer:
#'   * If the HR rasters are missing for that (id, season, year), skip with a
#'     warning. extract_step_variables can't run without them.
#'   * If the output rds already exists and `overwrite` is FALSE, skip.
#'   * Otherwise run make_random_pt_extraction -> extract_step_variables and
#'     save the result. Errors are caught per deer so a single failure does
#'     not halt the loop.
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

# landscape data — landcover and NDVI are now annual, so both are loaded per
# deer inside the loop (see load_landcover / make_water_mask / load_ndvi).
# No env_old.

# Deer movement data
sw_deer_tracks <- readRDS('library/SW_filtered_deer.RData')

# Filter — must match create_hr_rasters.R so the HR rasters line up by deer key
deer_mvt <- sw_deer_tracks %>%
  dplyr::filter(keep == T & excursion == F & unstable_hr_center == F) %>%
  dplyr::filter(year %in% 2017:2021)

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
      # Season-specific annual landcover (+ open-water mask) and NDVI stack
      landcover <- load_landcover(one_deer$year, one_deer$season)
      water <- make_water_mask(landcover)
      ndvi <- load_ndvi(one_deer$year)

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
        random_col = "stp.random",
        output_col = "stp.var"
      )

      # Non-parametric design: uniform-disc random steps (required for
      # non-parametric movement smooths s(sl_)/s(ta_); also usable for any GAM)
      # -> stp.random.nonp / stp.var.nonp
      one_deer <- make_random_pt_extraction(
        data = one_deer,
        n_pts = 50,
        water = water,
        model = "nonp",
        stp_col = "stp",
        output_col = "stp.random.nonp"
      )
      one_deer <- extract_step_variables(
        data = one_deer,
        env = landcover,
        ndvi = ndvi,
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
