#' @description
#' Out-of-sample simulation: simulate movement using a model trained on
#' (id, season, year) but with environmental inputs and observed step
#' structure from (id, season, year + 1).
#'
#' Usage: Rscript run_sims_test.R <id> <season> <year>
#'   id     — deer ID
#'   season — season string (e.g. "fa", "nb")
#'   year   — training year (test year is year + 1)
#'
#' Exits with status 2 if no test-year wrangled track exists for this deer.

# Parse command line arguments -------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Usage: Rscript run_sims_test.R <id> <season> <year>\nExample:
      Rscript run_sims_test.R 5000 fa 2017"
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

# Guard: bail out cleanly if there is no test-year track for this deer --------
test_track_path <- sprintf("data/tracks/data_%s.rds", test_key)
if (!file.exists(test_track_path)) {
  message(sprintf("no test data available for %s", test_key))
  quit(status = 2)
}

# Load packages ----------------------------------------------------------------
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
deer_mvt <- readRDS(test_track_path)

# landscape data — TEST-year season-specific landcover (categorical band + per-
# class binary indicator layers). env_old terrain covariates are dropped.
env_raster <- load_landcover(test_year, season)

# NDVI data — TEST year
ndvi_year <- load_ndvi(test_year)

# amt models — fitted on TRAIN year
results_amt <- readRDS(sprintf("results/amt/results_amt_%s.rds", train_key))

# Simulate every model that was fitted — null/failed models are filtered out
# automatically inside the precompute step (returns NULL for those).
n_models <- length(results_amt)
models_to_sim <- seq_len(n_models)
n_sim <- 10

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
# Same convention as run_sims.R: no HR_bin; HR_edge in metres; HR_center is
# transformed via log1p to match the formulas.
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

# Precompute simulation models ------------------------------------------------
cat("Precomputing simulation models...\n")

model_sims <- purrr::map(models_to_sim, function(m) {
  iss_i <- results_amt[[m]]$iss

  if (is.character(iss_i)) {
    return(NULL)
  }

  # Test-year stp.var supplies the model.matrix template for coefficient-name
  # matching. Factor levels are set explicitly in extract_step_variables, so
  # train- and test-year stp.var are interchangeable here.
  train_i <- deer_mvt$stp.var[[1]]

  coefs <- iss_i$model$coefficients
  names(coefs) <- rename_landcover_coefs(names(coefs))
  coefs <- coefs[!is.na(coefs)]

  dummy_sim <- amt::make_issf_model(coefs = coefs)
  mm_names <- colnames(model.matrix(
    amt:::ssf_formula(dummy_sim$model$formula),
    data = train_i
  ))

  for (idx in seq_along(coefs)) {
    if (grepl(":", names(coefs)[idx]) && !(names(coefs)[idx] %in% mm_names)) {
      parts <- strsplit(names(coefs)[idx], ":")[[1]]
      perms <- combinat::permn(parts)
      for (p in perms) {
        candidate <- paste(p, collapse = ":")
        if (candidate %in% mm_names) {
          names(coefs)[idx] <- candidate
          break
        }
      }
    }
  }

  amt::make_issf_model(
    coefs = coefs,
    sl = iss_i$sl_,
    ta = iss_i$ta_
  )
})
names(model_sims) <- as.character(models_to_sim)

# Free large objects
rm(env_raster, ndvi_year, results_amt, deer_mvt)
gc()

# Simulate across models in parallel ------------------------------------------
cat("Simulating movement...\n")

future::plan(
  multisession,
  workers = min(length(models_to_sim), parallel::detectCores() - 1)
)

results_sim <- furrr::future_map(
  models_to_sim,
  function(m) {
    cat("Simulating model:", m, "\n")

    model_sim <- model_sims[[as.character(m)]]

    if (is.null(model_sim)) {
      return(NA)
    }

    env_local <- terra::unwrap(deer_input$crop_env)
    ndvi_local <- terra::unwrap(deer_input$crop_ndvi)

    foreach(h = 1:n_sim, .combine = "rbind") %do%
      {
        res <- simulate_movement(
          stp_data = deer_input$stp,
          env_test = env_local,
          ndvi_test = ndvi_local,
          model = model_sim,
          method = "amt"
        )
        res$nsim <- h
        res
      }
  },
  .options = furrr::furrr_options(
    packages = c("amt", "terra", "sf", "dplyr", "lubridate", "foreach"),
    stdout = FALSE,
    seed = TRUE
  )
)

future::plan(sequential)
gc()

names(results_sim) <- as.character(models_to_sim)

# Save -------------------------------------------------------------------------
dir.create("sims/amt", showWarnings = FALSE, recursive = TRUE)
saveRDS(results_sim, sprintf("sims/amt/sims_amt_test_%s.rds", train_key))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf(
  "Deer %s (test on %s) completed in %.1f minutes\n",
  train_key,
  test_key,
  elapsed
))
