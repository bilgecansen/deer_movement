#' @description
#' Out-of-sample one-step-ahead log score for GAMs: score the test-year observed
#' track under a GAM trained on (id, season, year) using environmental inputs
#' (HR rasters, NDVI) from year + 1. GAM analogue of run_logscore_test.R.
#'
#' Like the in-sample runner, this scores the null model alongside the numbered
#' ones in a single pass and writes them to separate files:
#'   filters/gam/logscore_gam_test_<train_key>.rds       numbered (+ delta_logp)
#'   filters/gam/logscore_gam_null_test_<train_key>.rds  the null alone
#'
#' Usage: Rscript run_logscore_gam_test.R <id> <season> <year>
#'   id     — deer ID
#'   season — season string (e.g. "fa", "nb")
#'   year   — training year (test year is year + 1)
#'
#' Assumes results/gam/results_gam_<train_key>.rds and
#' results/gam/results_gam_null_<train_key>.rds, the TRAIN-year wrangled track
#' (data_<train_key>.rds, for the tentative gamma / von Mises distributions the
#' GAM corrections are relative to), and the TEST-year wrangled track all exist;
#' the bash wrapper gates on those so this script needs no "no test data" guard.

# Parse command line arguments -------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste0(
      "Usage: Rscript run_logscore_gam_test.R <id> <season> <year>\n",
      "Example: Rscript run_logscore_gam_test.R 5004 fa 2017"
    )
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

# GAM models — fitted on TRAIN year; a named list keyed by model number, plus
# the separately-saved null model (a single fit_gam_mod() return).
results_gam <- readRDS(sprintf("results/gam/results_gam_%s.rds", train_key))
results_gam_null <- readRDS(sprintf(
  "results/gam/results_gam_null_%s.rds",
  train_key
))

# Tentative movement distributions for the GAM kernel (gamma / von Mises). Use
# the TRAIN-year stp.var attributes — the same proposal the GAM was fit against,
# so the movement-kernel compensation matches the fitted corrections. (The amt
# test runner gets the equivalent from iss_i$sl_/iss_i$ta_ baked into the fitted
# object; an mgcv GAM does not carry them, so we read the train track here.)
train_mvt <- readRDS(sprintf("data/tracks/data_%s.rds", train_key))
sl_distr <- attr(train_mvt$stp.var[[1]], "sl_")
ta_distr <- attr(train_mvt$stp.var[[1]], "ta_")
rm(train_mvt)

# Score every fitted model — failed fits ($gam == "Error") return NULL below and
# are recorded as NA. Pre-naming-contract results files are positional (no
# names); fall back to positional numbering so the script still runs on them.
if (is.null(names(results_gam))) {
  names(results_gam) <- as.character(seq_along(results_gam))
}

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
# No coefficient surgery (unlike the amt path): the fitted mgcv cox.ph object
# is used directly by redistribution_kernel_gam. Failed fits become NULL.
cat("Collecting fitted GAMs...\n")

model_gams <- purrr::map(results_gam, function(r) {
  g <- r$gam
  if (is.character(g)) NULL else g
})

null_gam <- if (is.character(results_gam_null$gam)) {
  NULL
} else {
  results_gam_null$gam
}

# One scoring pass over the null plus every numbered model; see
# run_logscore_gam.R for why they share a run. The null is carried as a named
# unit ("null") so it can never be mistaken for a model number.
score_units <- c(list(null = null_gam), model_gams)
units_to_run <- names(score_units)

# Free large objects
rm(env_raster, ndvi_year, results_gam, results_gam_null, deer_mvt)
gc()

# Run log score across all models in parallel ---------------------------------
cat("Computing one-step log scores...\n")

future::plan(
  multisession,
  workers = min(length(units_to_run), parallel::detectCores() - 1)
)

results_logscore <- furrr::future_map(
  units_to_run,
  function(m) {
    cat("  Model:", m, "\n")

    model_gam <- score_units[[m]]

    if (is.null(model_gam)) {
      return(data.frame(
        model = m,
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
      model = m,
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

# Split null from numbered, and compute delta_logp vs. the null ---------------
# See run_logscore_gam.R: delta_logp lives on the numbered rows because it is a
# property of the comparison; the null's own scores are saved untouched.
results_all <- dplyr::bind_rows(results_logscore)

null_results <- results_all %>% dplyr::filter(model == "null")
null_logp <- null_results$total_logp

results <- results_all %>%
  dplyr::filter(model != "null") %>%
  dplyr::mutate(model = as.integer(model)) %>%
  dplyr::arrange(model) %>%
  dplyr::mutate(delta_logp = total_logp - null_logp)

cat("Null:\n")
print(null_results)
cat("Numbered models:\n")
print(results)

# Save -------------------------------------------------------------------------
dir.create("filters/gam", showWarnings = FALSE, recursive = TRUE)
saveRDS(results, sprintf("filters/gam/logscore_gam_test_%s.rds", train_key))
saveRDS(
  null_results,
  sprintf("filters/gam/logscore_gam_null_test_%s.rds", train_key)
)

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf(
  "Deer %s (test on %s) completed in %.1f minutes\n",
  train_key,
  test_key,
  elapsed
))
