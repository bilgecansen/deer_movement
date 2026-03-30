#' @description
#' Calculate movement prediction skill
#' 1. Estimate and compare utilization distributions
#' 2. Calculate Energy Scores

# Parse command line arguments -------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage: Rscript calc_pred_skill.R <start_row> <end_row>\nExample: Rscript calc_pred_skill.R 1 30"
  )
}

row_start <- as.integer(args[1])
row_end <- as.integer(args[2])

cat(sprintf("Running rows %d to %d\n", row_start, row_end))

# Load packages ----------------------------------------------------------------
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

# Observed and simulated deers
deer_mvt <- readRDS("data_deer_1_119.rds") %>%
  dplyr::slice(row_start:row_end)

results_sim <- readRDS(sprintf("results_sim_%d_%d.rds", row_start, row_end))

n_models <- length(results_sim)
n_deer <- nrow(deer_mvt)
n_sim <- 10
results_skill <- list()

# Step 1: Estimate UD overlap
cat("Estimating overlap of UDs and CTMMs...\n")

results_skill$ud <- purrr::map(1:n_models, function(m) {
  cat("UD and CTMM overlap for model:", m, "\n")

  sim_m <- results_sim[[m]]
  stp_test <- deer_mvt$stp_test

  future::plan(multisession, workers = parallel::detectCores() - 1)

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

  # Reset workers to free memory
  future::plan(sequential)
  gc()

  ud_m
})

names(results_skill$ud) <- 1:n_models

# Step 2: Calculate Energy Scores
cat("Calculating Energy Scores...\n")

results_skill$es <- foreach(m = 1:n_models, .combine = "rbind") %do%
  {
    cat("Energy Score for model:", m, "\n")

    scores <- purrr::map_dbl(1:n_deer, function(i) {
      sim_i <- results_sim[[m]][[i]]

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
      model = m,
      deer = deer_mvt$id,
      energy_score = scores
    )
  }

results_skill$es <- results_skill$es %>%
  dplyr::group_by(deer) %>%
  dplyr::mutate(energy_skill = 1 - (energy_score / energy_score[2]))

saveRDS(results_skill, sprintf("results_skill_%d_%d.rds", row_start, row_end))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf(
  "Rows %d-%d completed in %.1f minutes\n",
  row_start,
  row_end,
  elapsed
))
