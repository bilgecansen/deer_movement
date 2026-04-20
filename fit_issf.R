#' @description
#' Fit deer iSSF models for a single deer.
#'
#' Usage: Rscript fit_issf.R <row_number>
#'   row_number — deer index (row in data_deer_1_119.rds)

# Parse command line arguments -------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop(
    "Usage: Rscript fit_issf.R <row_number>\nExample: Rscript fit_issf.R 5"
  )
}

row_no <- as.integer(args[1])
cat(sprintf("Fitting models for deer %d\n", row_no))

# Load packages ----------------------------------------------------------------
library(amt)
library(terra)
library(tidyverse)

# helper functions
source("helper_functions.R")

# Load data --------------------------------------------------------------------
start_time <- Sys.time()

# Deer movement data
deer_mvt <- readRDS("data_deer_1_119.rds") %>%
  dplyr::slice(row_no)

# Formulas ---------------------------------------------------------------------
formulas <-
  c(
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + HR_end",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + HR_end + HR_end:days",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + (log(sl_) + cos(ta_)):HR_start",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_:HR_start",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_:HR_start + HR_end + HR_end:days",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + (log(sl_) + cos(ta_)):HR_start +
      HR_end + HR_end:days",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + HR_end + east_end + east2_end +
      north_end + dist_end + wiscland_end + ndvi_end + wiscland_end:ndvi_end",
    "case_ ~ log(sl_) + cos(ta_) + east_end + east2_end +
      north_end + dist_end + wiscland_end + ndvi_end + wiscland_end:ndvi_end",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + HR_end + east_end + east2_end +
      north_end + dist_end + wiscland_end"
  )

formulas <- paste(formulas, "+ strata(step_id_)")

# Fit models -------------------------------------------------------------------
cat("Fitting models...\n")

ssf_data <- deer_mvt$stp.var.train[[1]]

results_issf <- purrr::map(seq_along(formulas), function(i) {
  cat("  Formula", i, "\n")
  fit_mod(ssf_data, formulas[i])
})

# Save results -----------------------------------------------------------------
saveRDS(results_issf, sprintf("results/results_issf_%d.rds", row_no))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("Deer %d completed in %.1f minutes\n", row_no, elapsed))
