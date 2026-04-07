#' @description
#' Calculate movement prediction skill for a single deer
#' 1. Estimate and compare utilization distributions
#' 2. Calculate Energy Scores

# Parse command line arguments -------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop(
    "Usage: Rscript calc_pred_skill.R <row_no>\nExample: Rscript calc_pred_skill.R 1"
  )
}

row_no <- as.integer(args[1])
cat(sprintf("Running deer %d\n", row_no))

# Load packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(ctmm)
  library(scoringRules)
  library(furrr)
  library(parallel)
})

# helper functions
source("helper_functions.R")

# Load data --------------------------------------------------------------------

start_time <- Sys.time()

# Observed deer data
deer_mvt <- readRDS("data_deer_1_119.rds") %>%
  dplyr::slice(row_no)

# Simulated paths for this deer
results_sim <- readRDS(sprintf("results/results_sim_%d.rds", row_no))

n_models <- length(results_sim)
n_sim <- 10
results_skill <- list()

# Step 1: Estimate UD overlap (parallelized across models)
future::plan(multisession, workers = parallel::detectCores() - 1)

results_skill$ud <- suppressMessages(suppressWarnings(
  furrr::future_map(
    1:n_models,
    function(m) {
      cat("  UD model:", m, "\n")
      sim_m <- results_sim[[m]]

      if (length(sim_m) == 1 && is.na(sim_m)) {
        return(list(bat_uds = NA, bat_ctmm = NA))
      }

      overlap_ud(deer_mvt$stp_test[[1]], sim_m, n_sim = n_sim)
    },
    .options = furrr_options(
      packages = c("sf", "tidyverse", "ctmm"),
      stdout = TRUE,
      seed = TRUE
    )
  )
))
names(results_skill$ud) <- 1:n_models

future::plan(sequential)
gc()

# Step 2: Calculate Energy Scores
results_skill$es <- suppressMessages(suppressWarnings(
  purrr::map_dfr(1:n_models, function(m) {
    sim_m <- results_sim[[m]]

    if (length(sim_m) == 1 && is.na(sim_m)) {
      return(data.frame(
        model = m,
        deer = deer_mvt$id,
        energy_score = NA_real_
      ))
    }

    data.frame(
      model = m,
      deer = deer_mvt$id,
      energy_score = calc_energy_score(
        obs = deer_mvt$stp_test[[1]],
        sim = sim_m
      )
    )
  })
))

saveRDS(results_skill, sprintf("results/results_skill_%d.rds", row_no))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("Deer %d completed in %.1f minutes\n", row_no, elapsed))
