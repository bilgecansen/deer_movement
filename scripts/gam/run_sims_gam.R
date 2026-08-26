#' @description
#' Simulate deer movement from the numbered GAM models for a single deer. GAM
#' analogue of scripts/amt/run_sims.R.
#'
#' Only the NUMBERED models are simulated. The null model (fit_GAM.R writes it
#' to results/gam/results_gam_null_<key>.rds) is deliberately left out: it is
#' scored by log score for comparison, not simulated. Everything built on
#' simulations — UD overlap, SVF, energy score, and hence the four filter gates
#' — is therefore numbered-models-only by construction.
#'
#' Usage: Rscript scripts/gam/run_sims_gam.R <id> <season> <year>
#'   id     — deer ID
#'   season — season string (e.g. "fa", "nb")
#'   year   — year (integer)

# Parse command line arguments -------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Usage: Rscript scripts/gam/run_sims_gam.R <id> <season> <year>\nExample:
      Rscript scripts/gam/run_sims_gam.R 5004 fa 2017"
  )
}

id <- args[1]
season <- args[2]
year <- as.integer(args[3])

key <- sprintf("%s_%s_%d", id, season, year)
cat(sprintf("Running deer %s\n", key))

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

# Deer movement data — single-row per-deer file
deer_mvt <- readRDS(sprintf("data/tracks/data_%s.rds", key))

# landscape data — season-specific annual landcover (categorical band + per-
# class binary indicator layers). env_old terrain covariates are dropped.
env_raster <- load_landcover(year, season)

# NDVI data
ndvi_year <- load_ndvi(year)

# GAM models — a named list keyed by model number ("1".."4")
results_gam <- readRDS(sprintf("results/gam/results_gam_%s.rds", key))

# Guard against being handed the null file. results_gam_null_<key>.rds holds a
# single fit_gam_mod() return (top-level $gam), not a list of them, so it would
# otherwise be silently misread as a one-model set named after its own fields.
if (!is.null(results_gam[["gam"]])) {
  stop(
    "Loaded a single fitted model, not the numbered model list. ",
    "The null model is not simulated; pass a numbered-model key."
  )
}

# Tentative movement distributions for the GAM kernel (gamma / von Mises). The
# GAMs are fit on the parametric stp.var design, so these ride along on stp.var
# as attributes and supply the movement-kernel compensation at simulation time.
sl_distr <- attr(deer_mvt$stp.var[[1]], "sl_")
ta_distr <- attr(deer_mvt$stp.var[[1]], "ta_")

# Simulate every numbered model — failed fits ($gam == "Error") return NA below.
# Pre-naming-contract results files are positional (no names); fall back to
# positional numbering so the script still runs on them.
if (is.null(names(results_gam))) {
  names(results_gam) <- as.character(seq_along(results_gam))
}
models_to_sim <- names(results_gam)
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

# Collect fitted GAMs ---------------------------------------------------------
# No coefficient surgery (unlike the amt path): the fitted mgcv cox.ph object
# is used directly by redistribution_kernel_gam. Failed fits become NULL.
cat("Collecting fitted GAMs...\n")

model_gams <- purrr::map(results_gam, function(r) {
  g <- r$gam
  if (is.character(g)) NULL else g
})

# Free large objects
rm(env_raster, ndvi_year, results_gam, deer_mvt)
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

    model_gam <- model_gams[[m]]

    if (is.null(model_gam)) {
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
          model = model_gam,
          method = "gam",
          sl_distr = sl_distr,
          ta_distr = ta_distr
        )
        res$nsim <- h
        res
      }
  },
  .options = furrr::furrr_options(
    packages = c("mgcv", "amt", "terra", "sf", "dplyr", "lubridate", "foreach"),
    stdout = FALSE,
    seed = TRUE
  )
)

future::plan(sequential)
gc()

names(results_sim) <- models_to_sim

# Save -------------------------------------------------------------------------
dir.create("sims/gam", showWarnings = FALSE, recursive = TRUE)
saveRDS(results_sim, sprintf("sims/gam/sims_gam_%s.rds", key))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("Deer %s completed in %.1f minutes\n", key, elapsed))
