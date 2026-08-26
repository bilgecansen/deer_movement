#' @description
#' Compute UD overlap (Bhattacharyya) and SVF agreement between observed and
#' GAM-simulated paths for every model that was simulated for a single deer.
#' GAM analogue of scripts/issf/run_udoverlap.R: it reads sims_gam_<key> and
#' writes udoverlap_gam_<key>. The metric itself (overlap_ud) is model-agnostic,
#' so only the input/output filenames differ from the iSSF version.
#'
#' Numbered models only — no filtering needed here: sims_gam_<key> contains only
#' the numbered models (run_sims_gam.R does not simulate the null), so whatever
#' is in the file is by definition the right set.
#'
#' Usage: Rscript scripts/gam/run_udoverlap_gam.R <id> <season> <year>
#'   id     — deer ID
#'   season — season string (e.g. "fa", "nb")
#'   year   — year (integer)

# Parse command line arguments -------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste0(
      "Usage: Rscript scripts/gam/run_udoverlap_gam.R ",
      "<id> <season> <year>\n",
      "Example: Rscript scripts/gam/run_udoverlap_gam.R 5004 fa 2017"
    )
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

# GAM-simulated paths for this deer — named list keyed by the model numbers that
# were actually simulated (run_sims_gam.R simulates every numbered model; the
# null is not simulated and so never appears here).
results_sim <- readRDS(sprintf("sims/gam/sims_gam_%s.rds", key))

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
        return(list(bat_uds = NA_real_, svf_agree = NA_real_))
      }

      # If overlap_ud errors (commonly due to a degenerate ctmm fit on one
      # or more simulated paths -> dimension-mismatch inside ctmm internals
      # like `R$par %*% tJ`), record NA for both metrics for this model
      # rather than letting the worker crash and halt the whole batch.
      tryCatch(
        overlap_ud(deer_mvt$stp[[1]], sim_m, n_sim = n_sim),
        error = function(e) {
          message(sprintf("  Model %d failed: %s", m, conditionMessage(e)))
          list(bat_uds = NA_real_, svf_agree = NA_real_)
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
dir.create("filters/gam", showWarnings = FALSE, recursive = TRUE)
saveRDS(results_ud, sprintf("filters/gam/udoverlap_gam_%s.rds", key))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("Deer %s completed in %.1f minutes\n", key, elapsed))
