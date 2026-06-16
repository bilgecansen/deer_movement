#' @description
#' Simulate deer movement from all fitted models for a single deer.
#'
#' Usage: Rscript run_sims.R <id> <season> <year>
#'   id     — deer ID
#'   season — season string (e.g. "fa", "nb")
#'   year   — year (integer)

# Parse command line arguments -------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Usage: Rscript run_sims.R <id> <season> <year>\nExample: 
      Rscript run_sims.R 5000 fa 2017"
  )
}

id <- args[1]
season <- args[2]
year <- as.integer(args[3])

key <- sprintf("%s_%s_%d", id, season, year)
cat(sprintf("Running deer %s\n", key))

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

# Deer movement data — single-row per-deer file
deer_mvt <- readRDS(sprintf("data/tracks/data_%s.rds", key))

# landscape data — season-specific annual landcover (categorical band + per-
# class binary indicator layers). env_old terrain covariates are dropped.
env_raster <- load_landcover(year, season)

# NDVI data
ndvi_year <- load_ndvi(year)

# issf models
results_issf <- readRDS(sprintf("results/results_issf_%s.rds", key))

# Simulate every model that was fitted — null/failed models are filtered out
# automatically inside the precompute step (returns NULL for those).
n_models <- length(results_issf)
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

# HR rasters — no HR_bin (no model uses it). HR_edge stays in metres; for
# HR_center we also build the log1p transformed layer the formulas reference.
env_cropped$HR_edge <- load_hr_edge_raster(id, season, year, env_cropped)
env_cropped$HR_center <- load_hr_center_raster(id, season, year, env_cropped)
env_cropped$HR_center_log <- log1p(env_cropped$HR_center)

deer_input <- list(
  crop_env = terra::wrap(env_cropped),
  crop_ndvi = terra::wrap(terra::crop(ndvi_year, crop_extent)),
  stp = deer_mvt$stp[[1]]
)

# Precompute simulation models ------------------------------------------------
cat("Precomputing simulation models...\n")

model_sims <- purrr::map(models_to_sim, function(m) {
  iss_i <- results_issf[[m]]$iss

  if (is.character(iss_i)) {
    return(NULL)
  }

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
rm(env_raster, ndvi_year, results_issf, deer_mvt)
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
          method = "issf"
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
dir.create("sims", showWarnings = FALSE)
saveRDS(results_sim, sprintf("sims/sims_%s.rds", key))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("Deer %s completed in %.1f minutes\n", key, elapsed))
