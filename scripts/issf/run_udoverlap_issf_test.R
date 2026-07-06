#' @description
#' Out-of-sample UD overlap + SVF agreement: compare test-year simulated paths
#' (from a model trained on year, simulated for year+1 by run_sims_test.R)
#' against the test-year observed track.
#'
#' Usage: Rscript run_udoverlap_test.R <id> <season> <year>
#'   id     — deer ID
#'   season — season string (e.g. "fa", "nb")
#'   year   — training year (test year is year + 1)
#'
#' Assumes sims/issf/sims_issf_test_<train_key>.rds and the test-year wrangled track
#' both exist; the bash wrapper iterates over the test sims directly so this
#' script does not need a "no test data" guard.

# Parse command line arguments -------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Usage: Rscript run_udoverlap_test.R <id> <season> <year>\n
    Example: Rscript run_udoverlap_test.R 5000 fa 2017"
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

# Observed deer data — TEST year
deer_mvt <- readRDS(sprintf("data/tracks/data_%s.rds", test_key))

# Simulated paths — test simulations from the train-year model
results_sim <- readRDS(sprintf("sims/issf/sims_issf_test_%s.rds", train_key))

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

      tryCatch(
        overlap_ud(deer_mvt$stp[[1]], sim_m, n_sim = n_sim),
        error = function(e) {
          message(sprintf(
            "    Model %d failed: %s",
            m,
            conditionMessage(e)
          ))
          list(bat_uds = NA_real_, svf_score = NA_real_)
        }
      )
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
dir.create("filters/issf", showWarnings = FALSE, recursive = TRUE)
saveRDS(results_ud, sprintf("filters/issf/udoverlap_issf_test_%s.rds", train_key))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf(
  "Deer %s (test on %s) completed in %.1f minutes\n",
  train_key,
  test_key,
  elapsed
))
