#' @description
#' Deer movement data wrangling — full track per deer (no train/test split).
#' Iterates sequentially over every deer matching the create_hr_rasters.R
#' filter and writes one rds per deer to data/tracks/data_<id>_<season>_<year>.rds.
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
overwrite <- FALSE   # set to TRUE to reprocess deer that already have output

# Load packages ---------------------------------------------------------------
library(amt)
library(terra)
library(tidyverse)
library(sf)

# helper functions
source("scripts/helper_functions.R")

# Load shared inputs ----------------------------------------------------------
start_time <- Sys.time()

# landscape data
env_raster <- terra::rast("data/env/wiscland/wiscland2_binary.tif")
env_old <- terra::rast("Example_code/Env_2017.tif")
env_raster <- terra::crop(env_raster, env_old) %>%
  terra::resample(env_old)

env_raster$ele <- env_old$ele
env_raster$east <- env_old$eastness
env_raster$east2 <- env_old$eastness2
env_raster$north <- env_old$northness
env_raster$dist <- env_old$fe.dist

water_binary <- terra::ifel(env_raster$wiscland == "water", 1, 0)
names(water_binary) <- "Water"

# NDVI data
ndvi_rasters <- list(
  ndvi_2017 = terra::rast('data/NDVI_year/NDVI_2017.tif'),
  ndvi_2018 = terra::rast('data/NDVI_year/NDVI_2018.tif'),
  ndvi_2019 = terra::rast('data/NDVI_year/NDVI_2019.tif'),
  ndvi_2020 = terra::rast('data/NDVI_year/NDVI_2020.tif'),
  ndvi_2021 = terra::rast('data/NDVI_year/NDVI_2021.tif')
)

# Deer movement data
sw_deer_tracks <- readRDS('data/SW_filtered_deer.RData')

# Filter — must match create_hr_rasters.R so the HR rasters line up by deer key
deer_mvt <- sw_deer_tracks %>%
  dplyr::filter(keep == T & excursion == F & unstable_hr_center == F) %>%
  dplyr::filter(year %in% 2017:2021)

cat(sprintf("Total qualifying deer: %d\n", nrow(deer_mvt)))

# Output dir
dir.create("data/tracks", recursive = TRUE, showWarnings = FALSE)

# Loop --------------------------------------------------------------------------
n_done    <- 0L
n_skipped <- 0L
n_failed  <- 0L

for (i in seq_len(nrow(deer_mvt))) {
  one_deer <- deer_mvt[i, ]
  key      <- sprintf("%s_%s_%d", one_deer$id, one_deer$season, one_deer$year)
  out_path <- sprintf("data/tracks/data_%s.rds", key)
  hr_path  <- sprintf("data/HR/HRbin_%s.tif",  key)

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
      one_deer <- make_random_pt_extraction(
        data       = one_deer,
        n_pts      = 10,
        water      = water_binary,
        stp_col    = "stp",
        output_col = "stp.random"
      )
      one_deer <- extract_step_variables(
        data       = one_deer,
        env        = env_raster,
        ndvi_list  = ndvi_rasters,
        random_col = "stp.random",
        output_col = "stp.var"
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
