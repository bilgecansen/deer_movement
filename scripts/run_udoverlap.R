#' @description
#' Compute UD overlap (Bhattacharyya) and SVF agreement between observed and
#' simulated paths for every model that was simulated for a single deer.
#'
#' Usage: Rscript run_udoverlap.R <id> <season> <year>
#'   id     — deer ID
#'   season — season string (e.g. "fa", "nb")
#'   year   — year (integer)

# Parse command line arguments -------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Usage: Rscript run_udoverlap.R <id> <season> <year>\nExample: Rscript run_udoverlap.R 5000 fa 2017"
  )
}

id <- args[1]
season <- args[2]
year <- as.integer(args[3])

key <- sprintf("%s_%s_%d", id, season, year)
cat(sprintf("Running UD overlap for deer %s\n", key))

# Load packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(ctmm)
  library(furrr)
  library(parallel)
})

# helper functions
source("scripts/helper_functions.R")

# Load data --------------------------------------------------------------------
start_time <- Sys.time()

# Observed deer data — single-row per-deer file
deer_mvt <- readRDS(sprintf("data/tracks/data_%s.rds", key))

# Simulated paths for this deer — named list keyed by the model numbers that
# were actually simulated (run_sims.R simulates every fitted model now).
results_sim <- readRDS(sprintf("sims/sims_%s.rds", key))

# Model numbers actually present in results_sim
model_nums <- as.integer(names(results_sim))
n_sim <- 10

# Compute UD overlap + SVF score per model in parallel ------------------------
cat("Computing UD overlap and SVF score per model...\n")

future::plan(
  multisession,
  workers = min(length(model_nums), parallel::detectCores() - 1)
)

results_ud <- suppressMessages(suppressWarnings(
  furrr::future_map(
    model_nums,
    function(m) {
      cat("  Model:", m, "\n")
      sim_m <- results_sim[[as.character(m)]]

      if (length(sim_m) == 1 && is.na(sim_m)) {
        return(list(bat_uds = NA_real_, svf_score = NA_real_))
      }

      overlap_ud(deer_mvt$stp[[1]], sim_m, n_sim = n_sim)
    },
    .options = furrr_options(
      packages = c("sf", "tidyverse", "ctmm"),
      stdout = TRUE,
      seed = TRUE
    )
  )
))
names(results_ud) <- as.character(model_nums)

future::plan(sequential)
gc()

# Save -------------------------------------------------------------------------
dir.create("filters", showWarnings = FALSE)
saveRDS(results_ud, sprintf("filters/udoverlap_%s.rds", key))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("Deer %s completed in %.1f minutes\n", key, elapsed))
