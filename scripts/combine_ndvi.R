library(terra)
library(tidyverse)
library(foreach)

merge_NDVI <- function(ndvi_folder = "HLS_NDVI_monthly/", year) {
  # Load and merge NDVI rasters
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

  start_month <- paste(year, "-01-01", sep = "")
  months <- seq(as.Date(start_month), by = "month", length.out = 12)

  time(z_all) <- months
  z_all
}

for (i in 2017:2018) {
  r <- merge_NDVI(year = i)

  writeRaster(
    r,
    filename = paste("NDVI", paste(i, ".tif", sep = ""), sep = "_"),
    overwrite = T
  )
}
