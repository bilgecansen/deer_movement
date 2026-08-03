#' @description
#' One-step-ahead log score for a single deer, computed for every numbered GAM
#' model in results/gam/results_gam_<key>.rds AND for the null model in
#' results/gam/results_gam_null_<key>.rds. GAM analogue of run_logscore.R.
#'
#' Unlike the simulation path, the null IS scored here — that is the whole point
#' of the log score, which is what the numbered models are compared against.
#' Null and numbered results are written to separate files:
#'   filters/gam/logscore_gam_<key>.rds       numbered models (+ delta_logp)
#'   filters/gam/logscore_gam_null_<key>.rds  the null model alone
#'
#' Both are scored in the SAME run, over the same cropped rasters and the same
#' observed steps. That is deliberate: delta_logp is only meaningful as a
#' like-for-like difference, and splitting the null into its own invocation would
#' repeat all the raster work while letting the two sides drift apart.
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

# GAM models — a named list keyed by model number ("1".."4"), plus the null
# model, which fit_GAM.R saves separately as a single fit_gam_mod() return.
results_gam <- readRDS(sprintf("results/gam/results_gam_%s.rds", key))
results_gam_null <- readRDS(sprintf(
  "results/gam/results_gam_null_%s.rds",
  key
))

# Tentative movement distributions for the GAM kernel (gamma / von Mises). The
# GAMs are fit on the parametric stp.var design, so these ride along on stp.var
# as attributes and supply the movement-kernel compensation at scoring time.
sl_distr <- attr(deer_mvt$stp.var[[1]], "sl_")
ta_distr <- attr(deer_mvt$stp.var[[1]], "ta_")

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

null_gam <- if (is.character(results_gam_null$gam)) NULL else results_gam_null$gam

# One scoring pass over the null plus every numbered model. The null is carried
# as a named unit ("null") rather than a model number so it can never be mistaken
# for one; the results are split apart again before saving.
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
# delta_logp stays on the numbered rows because it is a property of the
# comparison, not of the null: it is only defined against the null scored on the
# same steps in this same run. The null's own scores are saved untouched in
# their own file.
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
saveRDS(results, sprintf("filters/gam/logscore_gam_%s.rds", key))
saveRDS(null_results, sprintf("filters/gam/logscore_gam_null_%s.rds", key))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("Deer %s completed in %.1f minutes\n", key, elapsed))
