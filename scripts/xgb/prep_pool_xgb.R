#' @description
#' Build the pooled step table the xgboost importance run reads. One season
#' and year at a time; every qualifying deer goes into one table.
#'
#' Deer are kept when they have at least MIN_STEPS_PER_EDF observed steps per
#' effective degree of freedom of their GAM null
#' (results/gam/results_gam_null_<key>.rds), so the pool is not dominated by
#' deer whose own models are barely identifiable.
#'
#' Random points are regenerated here rather than read from data/tracks/,
#' because this path uses the uniform-disc design with N_PTS points. Uniform
#' availability keeps the movement term meaning the same thing for every deer
#' in the pool, which the per-deer fitted gamma does not.
#'
#' Columns written: the movement block (sl_, tod_day, cos_ta), the habitat
#' variables in HAB_VARS, landcover as integer codes for xgboost's
#' categorical splits, and the FAMD axes. FAMD stays NA off forest — the fit
#' handles it with a switch (see helpers_xgb.R); nothing is filled in.
#'
#' Output: data/xgb/pooled_<season>_<year>.rds, a list of
#'   * pooled — the step table, sorted by stratum with the observed step first
#'   * filter — steps, edf and steps-per-edf for every deer considered
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
SEASON <- "fa"
YEAR <- 2021L
# Random points per observed step, drawn uniformly over the disc
N_PTS <- 100L
# Keep a deer when observed steps / edf of its GAM null is at least this
MIN_STEPS_PER_EDF <- 10
# Habitat variables; landcover enters as one categorical column
HAB_VARS <- c(
  "HR_center_end",
  "ndvi_end",
  "landcover",
  "forest_edge_end",
  "elevation_end",
  "northness_end",
  "eastness_end"
)
# FAMD axes; NA off forest, so they are fit through the switch
FAMD_VARS <- sprintf("famd%d_end", 1:5)
overwrite <- TRUE # set to FALSE to keep an existing pooled file

# Load packages ---------------------------------------------------------------
library(amt)
library(terra)
library(tidyverse)
library(sf)
library(mgcv)

# helper functions
source("scripts/helper_functions.R")

key <- sprintf("%s_%d", SEASON, YEAR)
out_path <- sprintf("data/xgb/pooled_%s.rds", key)
dir.create("data/xgb", showWarnings = FALSE, recursive = TRUE)

if (!overwrite && file.exists(out_path)) {
  stop(sprintf("%s exists and overwrite is FALSE", out_path))
}

# Which deer qualify ----------------------------------------------------------
null_files <- list.files(
  "results/gam",
  pattern = sprintf("^results_gam_null_.*_%s[.]rds$", key),
  full.names = TRUE
)
if (!length(null_files)) {
  stop(sprintf("No GAM null models for %s in results/gam/", key))
}

filter_tab <- purrr::map_dfr(null_files, function(f) {
  gfit <- readRDS(f)$gam
  deer_key <- sub("^results_gam_null_(.*)[.]rds$", "\\1", basename(f))
  if (!inherits(gfit, "gam")) {
    return(tibble(key = deer_key, n_steps = NA_integer_, edf = NA_real_))
  }
  tibble(
    key = deer_key,
    # The Cox-PH response is cbind(times, stratum); column 2 is the stratum,
    # so its distinct values count the steps the null was fit on.
    n_steps = length(unique(gfit$model[[1]][, 2])),
    edf = sum(gfit$edf)
  )
}) |>
  mutate(
    steps_per_edf = n_steps / edf,
    keep = !is.na(steps_per_edf) & steps_per_edf >= MIN_STEPS_PER_EDF
  )

keys <- filter_tab$key[filter_tab$keep]
cat(sprintf(
  "%s: keeping %d of %d deer (>= %g steps per null edf), %d observed steps\n",
  key, length(keys), nrow(filter_tab), MIN_STEPS_PER_EDF,
  sum(filter_tab$n_steps[filter_tab$keep])
))

# Shared inputs ---------------------------------------------------------------
tracks <- readRDS("library/SW_filtered_deer.RData") |>
  dplyr::filter(keep == TRUE, excursion == FALSE,
                unstable_hr_center == FALSE) |>
  dplyr::filter(year %in% 2017:2021)
tracks$key <- sprintf("%s_%s_%d", tracks$id, tracks$season, tracks$year)

landcover <- load_landcover(YEAR, SEASON)
water <- make_water_mask(landcover)
ndvi <- load_ndvi(YEAR)
landfire <- load_landfire(YEAR, SEASON)
# Elevation, northness and eastness; the same layers for every deer and year.
# extract_step_variables() needs them, and they reach the model only if
# HAB_VARS names them.
topo <- load_topo()

# One row per (deer, step, point) ---------------------------------------------
set.seed(YEAR)
start_time <- Sys.time()

pooled <- build_xgb_pool(
  keys = keys,
  tracks = tracks,
  rasters = list(landcover = landcover, water = water, ndvi = ndvi,
                 landfire = landfire, topo = topo),
  n_pts = N_PTS,
  hab_vars = HAB_VARS,
  famd_vars = FAMD_VARS
)

on <- !is.na(pooled[[FAMD_VARS[1]]])
stopifnot(all(vapply(FAMD_VARS,
                     function(v) identical(!is.na(pooled[[v]]), on),
                     logical(1))))

saveRDS(list(pooled = pooled, filter = filter_tab), out_path)

cat(sprintf(
  "deer %d | steps %d | rows %d of %d kept | %.0f s\n",
  dplyr::n_distinct(pooled$deer), max(pooled$stratum), nrow(pooled),
  (N_PTS + 1) * sum(filter_tab$n_steps[filter_tab$keep]),
  as.numeric(difftime(Sys.time(), start_time, units = "secs"))
))
cat(sprintf(
  "FAMD present on %.1f%% of rows and %.1f%% of observed steps\n",
  100 * mean(on), 100 * mean(on[pooled$case_ == 1])
))
cat(sprintf("-> %s\n", out_path))
