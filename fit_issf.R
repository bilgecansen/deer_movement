#' @description
#' Deer movement analysis code
#' 1. Fit iSSF models
#' 2. Simulate movement from models
#' 3. Estimate and compare utilization distributions
#' 4. Calculate Energy Scores

# Parse command line arguments -------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage: Rscript fit_issf.R <start_row> <end_row>\nExample: Rscript fit_issf.R 1 30"
  )
}

row_start <- as.integer(args[1])
row_end <- as.integer(args[2])

cat(sprintf("Running rows %d to %d\n", row_start, row_end))

# Load packages ----------------------------------------------------------------
library(amt)
library(terra)
library(tidyverse)
library(foreach)
library(sf)
library(ctmm)
library(scoringRules)
library(furrr)
library(parallel)

# helper functions
source("helper_functions.R")

# Load data --------------------------------------------------------------------

start_time <- Sys.time()

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

water_binary <- ifel(env_raster$wiscland == "water", 1, 0)
names(water_binary) <- "Water"

# NDVI data
ndvi_rasters <- list(
  ndvi_2017 = terra::rast('NDVI_2017.tif'),
  ndvi_2018 = terra::rast('NDVI_2018.tif')
)

# Deer movement data
deer_mvt <- readRDS("data_deer_1_119.rds")

# Separate training and test splits
deer_mvt <- deer_mvt %>%
  dplyr::mutate(
    stp_train = purrr::map(stp, ~ .x[1:ceiling(nrow(.x) * 0.7), ]),
    stp_test = purrr::map(stp, ~ .x[(ceiling(nrow(.x) * 0.7) + 1):nrow(.x), ])
  ) %>%
  dplyr::slice(row_start:row_end)

# Main workflow ----------------------------------------------------------------

# Step 1: Fit models
cat("Fitting models...\n")

formulas <-
  c(
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + HR_end",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + HR_end:days",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + (log(sl_) + cos(ta_)):HR_start",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_:HR_start",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_:HR_start + HR_end:days",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + (log(sl_) + cos(ta_)):HR_start + 
      HR_end:days",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + HR_end + east_end + east2_end + 
      north_end + dist_end + wiscland_end + wiscland_end:ndvi_end",
    "case_ ~ log(sl_) + cos(ta_) + east_end + east2_end + 
      north_end + dist_end + wiscland_end + wiscland_end:ndvi_end"
  )

formulas <- paste(formulas, "+ strata(step_id_)")

formula_df <- data.frame(
  formula = formulas,
  name = 1:length(formulas)
)

## Fit models for each formula
results <- list()

results$models <- purrr::map(1:nrow(formula_df), function(i) {
  cat("Fitting formula", i, "\n")
  fit_mods(deer_mvt, formula_df[i, 1])
})

# Step 2: Simulate movement
cat("Simulating movement...\n")

n_models <- nrow(formula_df)
n_deer <- nrow(deer_mvt)
n_sim <- 10

future::plan(multisession, workers = parallel::detectCores() - 1)

results$sim <- purrr::map(1:n_models, function(m) {
  cat("Simulating model:", formula_df$name[m], "\n")

  run_simulations(
    deer_mvt,
    results$models[[m]],
    env_raster,
    ndvi_rasters,
    n_sim = n_sim,
    n_deer = nrow(deer_mvt)
  )
})
names(results$sim) <- formula_df$name

# Step 3: Estimate UD overlap
cat("Estimating overlap of UDs and CTMMs...\n")

results$ud <- purrr::map(1:n_models, function(m) {
  cat("UD and CTMM overlap for model:", formula_df$name[m], "\n")

  sim_m <- results$sim[[m]]
  stp_test <- deer_mvt$stp_test

  ud_m <- furrr::future_map(
    1:n_deer,
    function(i) {
      overlap_ud(stp_test[[i]], sim_m[[i]], n_sim = n_sim)
    },
    .options = furrr_options(
      packages = c("amt", "terra", "sf", "tidyverse", "ctmm"),
      stdout = FALSE,
      seed = T
    )
  )

  gc()
  ud_m
})

names(results$ud) <- formula_df$name

# Step 4: Calculate Energy Scores
cat("Calculating Energy Scores...\n")

results$es <- foreach(m = 1:n_models, .combine = "rbind") %do%
  {
    cat("Energy Score for model:", formula_df$name[m], "\n")

    scores <- purrr::map_dbl(1:n_deer, function(i) {
      sim_i <- results$sim[[m]][[i]]

      # Skip if simulation failed
      if (length(sim_i) == 1 && is.na(sim_i)) {
        return(NA_real_)
      }

      calc_energy_score(
        obs = deer_mvt$stp_test[[i]],
        sim = sim_i
      )
    })

    data.frame(
      model = formula_df$name[m],
      deer = deer_mvt$id,
      energy_score = scores
    )
  }

results$es <- results$es %>%
  dplyr::group_by(deer) %>%
  dplyr::mutate(energy_skill = 1 - (energy_score / energy_score[2]))

saveRDS(results, sprintf("results_issf_%d_%d.rds", row_start, row_end))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf(
  "Rows %d-%d completed in %.1f minutes\n",
  row_start,
  row_end,
  elapsed
))
