#' @description
#' Out-of-sample one-step-ahead log score: score the test-year observed track
#' under a model trained on (id, season, year) using environmental inputs
#' (HR rasters, NDVI) from year + 1.
#'
#' Usage: Rscript run_logscore_test.R <id> <season> <year>
#'   id     — deer ID
#'   season — season string (e.g. "fa", "nb")
#'   year   — training year (test year is year + 1)
#'
#' Assumes both results/issf/results_issf_<train_key>.rds and the test-year
#' wrangled track exist; the bash wrapper gates on those so this script does
#' not need a "no test data" guard.

# Parse command line arguments -------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste0(
      "Usage: Rscript run_logscore_test.R <id> <season> <year>\n",
      "Example: Rscript run_logscore_test.R 5000 fa 2017"
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

# issf models — fitted on TRAIN year
results_issf <- readRDS(sprintf("results/issf/results_issf_%s.rds", train_key))

# Score every fitted model — null/failed models filtered out automatically
# inside the precompute step (returns NULL for those).
n_models <- length(results_issf)
models_to_run <- seq_len(n_models)

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
# Same convention as run_logscore.R: no HR_bin; HR_edge in metres; HR_center
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

# Precompute simulation models ------------------------------------------------
cat("Precomputing simulation models...\n")

model_sims <- purrr::map(models_to_run, function(m) {
  iss_i <- results_issf[[m]]$iss

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
names(model_sims) <- as.character(models_to_run)

# Free large objects
rm(env_raster, ndvi_year, results_issf, deer_mvt)
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

    model_sim <- model_sims[[as.character(m)]]

    if (is.null(model_sim)) {
      return(data.frame(
        model = m,
        total_logp = NA_real_,
        n_steps = 0L
      ))
    }

    env_local <- terra::unwrap(deer_input$crop_env)
    ndvi_local <- terra::unwrap(deer_input$crop_ndvi)

    ll_df <- onestep_logscore(
      stp_data = deer_input$stp,
      env_test = env_local,
      ndvi_test = ndvi_local,
      issf_train = model_sim
    )

    data.frame(
      model = m,
      total_logp = sum(ll_df$logp, na.rm = TRUE),
      n_steps = sum(!is.na(ll_df$logp))
    )
  },
  .options = furrr::furrr_options(
    packages = c("amt", "terra", "sf", "dplyr", "lubridate", "foreach"),
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
dir.create("filters/issf", showWarnings = FALSE, recursive = TRUE)
saveRDS(results, sprintf("filters/issf/logscore_issf_test_%s.rds", train_key))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf(
  "Deer %s (test on %s) completed in %.1f minutes\n",
  train_key,
  test_key,
  elapsed
))
