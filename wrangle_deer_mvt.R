#' @description
#' Deer movement data wrangling
#' 1. Generate random steps
#' 2. Extract environmental variables for all steps

# Parse command line arguments -------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage: Rscript wrangle_deer_mvt.R <start_row> <end_row>\nExample: Rscript wrangle_deer_mvt.R 1 30"
  )
}

row_start <- as.integer(args[1])
row_end <- as.integer(args[2])

cat(sprintf("Running rows %d to %d\n", row_start, row_end))

# Load packages ----------------------------------------------------------------
library(amt)
library(terra)
library(tidyverse)
library(sf)
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

water_binary <- terra::ifel(env_raster$wiscland == "water", 1, 0)
names(water_binary) <- "Water"

# NDVI data
ndvi_rasters <- list(
  ndvi_2017 = terra::rast('NDVI_2017.tif'),
  ndvi_2018 = terra::rast('NDVI_2018.tif')
)

# Deer movement data
raw_data <- readRDS('Example_code/SW_filtered_deer.RData')

# Prepare sample data
deer_mvt <- raw_data %>%
  dplyr::filter(keep == T & excursion == F & unstable_hr_center == F) %>%
  dplyr::filter(year %in% c(2017, 2018))

# Separate training and test splits
deer_mvt <- deer_mvt %>%
  dplyr::mutate(
    stp_train = purrr::map(stp, ~ .x[1:ceiling(nrow(.x) * 0.7), ]),
    stp_test = purrr::map(stp, ~ .x[(ceiling(nrow(.x) * 0.7) + 1):nrow(.x), ])
  ) %>%
  dplyr::slice(row_start:row_end)


# Main workflow ----------------------------------------------------------------

future::plan(multisession, workers = parallel::detectCores() - 1)

# Step 1: Generate random steps
cat("Generating random steps for train...\n")
deer_mvt <- make_random_pt_extraction(
  data = deer_mvt,
  n_pts = 10,
  water = water_binary,
  stp_col = "stp_train",
  output_col = "stp.random.train"
)

cat("Generating random steps for test...\n")
deer_mvt <- make_random_pt_extraction(
  data = deer_mvt,
  n_pts = 10,
  water = water_binary,
  stp_col = "stp_test",
  output_col = "stp.random.test"
)

# Step 2: Extract environmental variables
cat("Extracting environmental variables for train...\n")
deer_mvt <- extract_step_variables(
  data = deer_mvt,
  env = env_raster,
  ndvi_list = ndvi_rasters,
  random_col = "stp.random.train",
  output_col = "stp.var.train"
)

cat("Extracting environmental variables for test...\n")
deer_mvt <- extract_step_variables(
  data = deer_mvt,
  env = env_raster,
  ndvi_list = ndvi_rasters,
  random_col = "stp.random.test",
  output_col = "stp.var.test"
)

saveRDS(deer_mvt, sprintf("data_deer_%d_%d.rds", row_start, row_end))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf(
  "Rows %d-%d completed in %.1f minutes\n",
  row_start,
  row_end,
  elapsed
))
