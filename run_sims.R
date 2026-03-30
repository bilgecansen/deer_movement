#' @description
#' Simulate deer movement from models

# Parse command line arguments -------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop(
    "Usage: Rscript run_sims.R <row_number>\nExample: Rscript run_sims.R 5"
  )
}

row_no <- as.integer(args[1])

cat(sprintf("Running deer %d\n", row_no))

# Load packages ----------------------------------------------------------------
library(amt)
library(terra)
library(tidyverse)
library(foreach)
library(sf)
library(furrr)
library(parallel)

# helper functions
source("helper_functions.R")

# Load data --------------------------------------------------------------------
start_time <- Sys.time()

# Deer movement data
deer_mvt <- readRDS("data_deer_1_119.rds")

# Pick deers to simulate
deer_mvt <- deer_mvt %>%
  dplyr::slice(row_no)

# landscape data
env_raster <- terra::rast("env/wiscland/wiscland2_binary.tif")
env_old <- terra::rast("Example_code/Env_2017.tif")
env_raster <- terra::crop(env_raster, env_old) %>%
  terra::resample(env_old)

env_raster$ele <- env_old$ele
env_raster$east <- env_old$eastness
env_raster$east2 <- env_old$eastness2
env_raster$north <- env_old$northness
env_raster$dist <- env_old$fe.dist

# NDVI data
ndvi_year <- terra::rast(paste(
  paste("NDVI", deer_mvt$year, sep = "_"),
  ".tif",
  sep = ""
))


# issf models
results_issf <- readRDS(sprintf("results_issf_1_119.rds"))

# Simulate movement ------------------------------------------------------------

n_models <- length(results_issf)
n_deer <- nrow(deer_mvt)
n_sim <- 10

# Pre-crop rasters for this single deer
crop_extent <- sf::st_buffer(
  sf::st_as_sf(
    deer_mvt$stp_test[[1]],
    coords = c('x1_', 'y1_'),
    crs = 6610
  ),
  5000
)

deer_input <- list(
  crop_env = terra::wrap(terra::crop(env_raster, crop_extent)),
  crop_ndvi = terra::wrap(terra::crop(ndvi_year, crop_extent)),
  stp_test = deer_mvt$stp_test[[1]],
  x_median = deer_mvt$x_median,
  y_median = deer_mvt$y_median
)

# Precompute simulation models for all models for this deer
cat("Precomputing simulation models...\n")

model_sims <- purrr::map(1:n_models, function(m) {
  coeff_i <- results_issf[[m]]$coeff[[1]]

  if (length(coeff_i) == 1) {
    return(NULL)
  }

  iss_i <- results_issf[[m]]$iss[[1]]
  train_i <- deer_mvt$stp.var.train[[1]]

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

# Free large objects
rm(env_raster, ndvi_year, env_old, results_issf, deer_mvt)
gc()

# Simulate across models in parallel
cat("Simulating movement...\n")

future::plan(multisession, workers = parallel::detectCores() - 1)

results_sim <- furrr::future_map(
  1:n_models,
  function(m) {
    cat("Simulating model:", m, "\n")

    if (is.null(model_sims[[m]])) {
      return(NA)
    }

    env_local <- terra::unwrap(deer_input$crop_env)
    ndvi_local <- terra::unwrap(deer_input$crop_ndvi)

    foreach(h = 1:n_sim, .combine = "rbind") %do%
      {
        res <- simulate_movement(
          stp_data = deer_input$stp_test,
          x_median = deer_input$x_median,
          y_median = deer_input$y_median,
          env_test = env_local,
          ndvi_test = ndvi_local,
          issf_train = model_sims[[m]]
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

names(results_sim) <- 1:n_models

saveRDS(results_sim, sprintf("results/results_sim_%d.rds", row_no))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf(
  "Row %d completed in %.1f minutes\n",
  row_no,
  elapsed
))
