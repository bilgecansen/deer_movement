#' @description
#' Merge raw monthly HLS NDVI tiles into one 12-layer raster per year, aligned
#' to the working env grid. Pre-aligning here means downstream scripts
#' (simulate_movement, onestep_loglik, extract_step_variables) don't need to
#' resample NDVI again at runtime.

library(terra)
library(tidyverse)
library(foreach)

# warp_to_template() lives in helper_functions.R
source("scripts/helper_functions.R")

# Spatial template — final NDVI grid matches this so downstream scripts can
# stack NDVI directly into the env raster without further resampling.
template <- terra::rast("data/env/wiscland/wiscland2.tif")


merge_NDVI <- function(
  year,
  template,
  outfile,
  ndvi_folder = "data/NDVI_month/"
) {
  # Load and merge raw monthly tiles
  idx <- stringr::str_detect(list.files(ndvi_folder), as.character(year))

  ndvi_files <- paste(
    ndvi_folder,
    list.files(ndvi_folder)[idx],
    sep = ""
  )

  z <- map(ndvi_files, terra::rast)

  z_all <- do.call(terra::mosaic, z)

  names(z_all) <- c(
    "Jan",
    "Feb",
    "March",
    "Apr",
    "May",
    "June",
    "July",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec"
  )

  # Mid-month timestamps, matching load_ndvi(). Stamping the 1st makes amt's
  # nearest-in-time layer lookup at fit time resolve to the FOLLOWING month for
  # any step past mid-month, which then disagrees with the calendar-month lookup
  # used at simulation / scoring time. See ndvi_layer_times() for the full note.
  time(z_all) <- ndvi_layer_times(year)

  # Warp onto the working grid AND write the output in one shot. Layer names
  # and time stamps travel via the .aux.xml sidecar that warp_to_template
  # copies from the temp input to the final output.
  warp_to_template(
    r = z_all,
    template = template,
    outfile = outfile,
    method = "near",
    datatype = "FLT4S"
  )
}


for (i in 2016:2021) {
  cat(sprintf("Merging NDVI for %d\n", i))

  outfile <- sprintf("data/NDVI_year/NDVI_%d.tif", i)

  merge_NDVI(year = i, template = template, outfile = outfile)
}
