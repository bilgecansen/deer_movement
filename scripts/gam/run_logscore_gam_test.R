#' @description
#' Out-of-sample one-step-ahead log score for GAMs: score the test-year observed
#' track under a GAM trained on (id, season, year) using environmental inputs
#' (HR rasters, NDVI) from year + 1. GAM analogue of run_logscore_test.R.
#'
#' Usage: Rscript run_logscore_gam_test.R <id> <season> <year>
#'   id     — deer ID
#'   season — season string (e.g. "fa", "nb")
#'   year   — training year (test year is year + 1)
#'
#' Assumes results/gam/results_gam_<train_key>.rds, the TRAIN-year wrangled track
#' (data_<train_key>.rds, for the tentative gamma / von Mises distributions the
#' GAM corrections are relative to), and the TEST-year wrangled track all exist;
#' the bash wrapper gates on those so this script needs no "no test data" guard.

# Parse command line arguments -------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Usage: Rscript run_logscore_gam_test.R <id> <season> <year>\nExample: Rscript run_logscore_gam_test.R 5004 fa 2017"
  )
}

id <- args[1]
season <- args[2]
year <- as.integer(args[3])
test_year <- year + 1L

train_key <- sprintf("%s_%s_%d", id, season, year)
test_key <- sprintf("%s_%s_%d", id, season, test_year)

cat(sprintf("Training key: %s\n", train_key))
cat(sprintf("Test key:     %s\n", test_key))

# Load packages ----------------------------------------------------------------
library(mgcv)
library(amt)
library(terra)
library(tidyverse)
library(foreach)
library(sf)
library(furrr)
library(parallel)

# helper functions
source("scripts/helper_functions.R")

# Load data --------------------------------------------------------------------
start_time <- Sys.time()

# Deer movement data — single-row per-deer file for the TEST year
deer_mvt <- readRDS(sprintf("data/tracks/data_%s.rds", test_key))

# landscape data — TEST-year season-specific landcover (categorical band + per-
# class binary indicator layers). env_old terrain covariates are dropped.
env_raster <- load_landcover(test_year, season)

# NDVI data — TEST year
ndvi_year <- load_ndvi(test_year)

# GAM models — fitted on TRAIN year; a named list keyed by model number
results_gam <- readRDS(sprintf("results/gam/results_gam_%s.rds", train_key))

# Tentative movement distributions for the GAM kernel (gamma / von Mises). Use
# the TRAIN-year stp.var attributes — the same proposal the GAM was fit against,
# so the movement-kernel compensation matches the fitted corrections. (The iSSF
# test runner gets the equivalent from iss_i$sl_/iss_i$ta_ baked into the fitted
# object; an mgcv GAM does not carry them, so we read the train track here.)
train_mvt <- readRDS(sprintf("data/tracks/data_%s.rds", train_key))
sl_distr <- attr(train_mvt$stp.var[[1]], "sl_")
ta_distr <- attr(train_mvt$stp.var[[1]], "ta_")
rm(train_mvt)

# Score every fitted model — failed fits ($gam == "Error") return NULL below and
# are recorded as NA. Older results_gam files are 8-element positional (no
# names); fall back to positional numbering so the script runs on both formats.
if (is.null(names(results_gam))) {
  names(results_gam) <- as.character(seq_along(results_gam))
}
models_to_run <- names(results_gam)

# Pre-crop rasters for this single deer ---------------------------------------
crop_extent <- sf::st_buffer(
  sf::st_as_sf(
    deer_mvt$stp[[1]],
    coords = c('x1_', 'y1_'),
    crs = 6610
  ),
  CROP_BUFFER_M
)

env_cropped <- terra::crop(env_raster, crop_extent)

# HR rasters — load TEST-year rasters (test_year keys the file names).
# Same convention as run_logscore_gam.R: no HR_bin; HR_edge in metres; HR_center
# is transformed via log1p to match the formulas.
env_cropped$HR_edge <- load_hr_edge_raster(id, season, test_year, env_cropped)
env_cropped$HR_center <- load_hr_center_raster(
  id,
  season,
  test_year,
  env_cropped
)
env_cropped$HR_center_log <- log1p(env_cropped$HR_center)

deer_input <- list(
  crop_env = terra::wrap(env_cropped),
  crop_ndvi = terra::wrap(terra::crop(ndvi_year, crop_extent)),
  stp = deer_mvt$stp[[1]]
)

# Collect fitted GAMs ---------------------------------------------------------
# No coefficient surgery (unlike the iSSF path): the fitted mgcv cox.ph object is
# used directly by redistribution_kernel_gam. Failed fits become NULL.
cat("Collecting fitted GAMs...\n")

model_gams <- purrr::map(results_gam, function(r) {
  g <- r$gam
  if (is.character(g)) NULL else g
})

# Free large objects
rm(env_raster, ndvi_year, results_gam, deer_mvt)
gc()

# Run log score across all models in parallel ---------------------------------
cat("Computing one-step log scores...\n")

future::plan(
  multisession,
  workers = min(length(models_to_run), parallel::detectCores() - 1)
)

results_logscore <- furrr::future_map(
  models_to_run,
  function(m) {
    cat("  Model:", m, "\n")

    model_gam <- model_gams[[m]]

    if (is.null(model_gam)) {
      return(data.frame(
        model = as.integer(m),
        total_logp = NA_real_,
        n_steps = 0L
      ))
    }

    env_local <- terra::unwrap(deer_input$crop_env)
    ndvi_local <- terra::unwrap(deer_input$crop_ndvi)

    ll_df <- onestep_logscore_gam(
      stp_data = deer_input$stp,
      env_test = env_local,
      ndvi_test = ndvi_local,
      gam_train = model_gam,
      sl_distr = sl_distr,
      ta_distr = ta_distr
    )

    data.frame(
      model = as.integer(m),
      total_logp = sum(ll_df$logp, na.rm = TRUE),
      n_steps = sum(!is.na(ll_df$logp))
    )
  },
  .options = furrr::furrr_options(
    packages = c("mgcv", "amt", "terra", "sf", "dplyr", "lubridate", "foreach"),
    stdout = TRUE,
    seed = TRUE
  )
)

future::plan(sequential)
gc()

# Combine and compute delta_logp relative to model 2 --------------------------
results <- dplyr::bind_rows(results_logscore)

null_logp <- results$total_logp[results$model == 2]

results <- results %>%
  dplyr::mutate(delta_logp = total_logp - null_logp)

cat("Results:\n")
print(results)

# Save -------------------------------------------------------------------------
dir.create("filters/gam", showWarnings = FALSE, recursive = TRUE)
saveRDS(results, sprintf("filters/gam/logscore_gam_test_%s.rds", train_key))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf(
  "Deer %s (test on %s) completed in %.1f minutes\n",
  train_key,
  test_key,
  elapsed
))
