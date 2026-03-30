#' @description
#' Fit deer iSSF models

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

# Deer movement data
deer_mvt <- readRDS("data_deer_1_119.rds")

# Pick deers to fit models
deer_mvt <- deer_mvt %>%
  dplyr::slice(row_start:row_end)

# Fit issf models --------------------------------------------------------------

# Set formulas
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

# Fit models for each formula
results_issf <- purrr::map(1:nrow(formula_df), function(i) {
  cat("Fitting formula", i, "\n")
  fit_mods(deer_mvt, formula_df[i, 1])
})

saveRDS(results_issf, sprintf("results_issf_%d_%d.rds", row_start, row_end))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf(
  "Rows %d-%d completed in %.1f minutes\n",
  row_start,
  row_end,
  elapsed
))
