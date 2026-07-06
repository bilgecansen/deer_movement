#' @description
#' One-step-ahead log score for a single deer, computed for every fitted GAM
#' model in results/gam/results_gam_<key>.rds. GAM analogue of run_logscore.R.
#'
#' Usage: Rscript run_logscore_gam.R <id> <season> <year>
#'   id     — deer ID
#'   season — season string (e.g. "fa", "nb")
#'   year   — year (integer)

# Parse command line arguments -------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Usage: Rscript run_logscore_gam.R <id> <season> <year>\nExample: Rscript run_logscore_gam.R 5004 fa 2017"
  )
}

id <- args[1]
season <- args[2]
year <- as.integer(args[3])

key <- sprintf("%s_%s_%d", id, season, year)
cat(sprintf("Running GAM log score for deer %s\n", key))

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

# GAM models — a named list keyed by model number ("1".."6")
results_gam <- readRDS(sprintf("results/gam/results_gam_%s.rds", key))

# Tentative movement distributions for the GAM kernel (gamma / von Mises). The
# GAMs are fit on the parametric stp.var design, so these ride along on stp.var
# as attributes and supply the movement-kernel compensation at scoring time.
sl_distr <- attr(deer_mvt$stp.var[[1]], "sl_")
ta_distr <- attr(deer_mvt$stp.var[[1]], "ta_")

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
saveRDS(results, sprintf("filters/gam/logscore_gam_%s.rds", key))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("Deer %s completed in %.1f minutes\n", key, elapsed))
