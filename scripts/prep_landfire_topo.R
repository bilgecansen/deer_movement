#' @description
#' Download LANDFIRE 2020 elevation and aspect for the study area and put them
#' on the landcover grid, as library/landfire/elevation.tif and
#' library/landfire/aspect.tif.
#'
#' Source: LANDFIRE Product Service (LFPS) layers LF2020_Elev (metres) and
#' LF2020_Asp (whole degrees clockwise from true north, 0-359, with -1 marking
#' flat cells), 30 m, derived by LANDFIRE from the USGS 1 arc-second DEM. They
#' are fetched with rlandfire::landfireAPIv2 for the landcover grid's WGS84
#' bounding box plus 0.05 degrees, in LANDFIRE's default Albers projection,
#' and warped once onto the landcover grid:
#'   * elevation — bilinear, float
#'   * aspect    — nearest cell, because degrees wrap at 0/360 and -1 is a flag
#'
#' The LFPS API requires an email address. Pass it in the environment, e.g.
#'   LANDFIRE_EMAIL=you@example.edu Rscript scripts/prep_landfire_topo.R
#' The downloaded zip is kept in library/landfire/raw/ and reused if present.

# Load packages ---------------------------------------------------------------
library(terra)
library(rlandfire)

# helper functions
source("scripts/helper_functions.R")

email <- Sys.getenv("LANDFIRE_EMAIL")
if (!nzchar(email)) {
  stop(
    "Set LANDFIRE_EMAIL to the address for the LFPS request, e.g. ",
    "LANDFIRE_EMAIL=you@example.edu Rscript scripts/prep_landfire_topo.R"
  )
}

out_dir <- "library/landfire"
raw_dir <- file.path(out_dir, "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
zip_path <- file.path(raw_dir, "lf2020_topo.zip")

lc_files <- list.files(
  "library/landcover",
  "^landcover_\\d{4}\\.tif$",
  full.names = TRUE
)
template <- terra::rast(lc_files[1])[[1]]

# Download --------------------------------------------------------------------
if (!file.exists(zip_path)) {
  bbox <- terra::as.polygons(terra::ext(template), crs = terra::crs(template))
  bb <- as.vector(terra::ext(terra::project(bbox, "EPSG:4326")))
  aoi <- round(c(bb[1] - 0.05, bb[3] - 0.05, bb[2] + 0.05, bb[4] + 0.05), 4)
  resp <- rlandfire::landfireAPIv2(
    c("LF2020_Elev", "LF2020_Asp"),
    aoi,
    email,
    path = zip_path,
    verbose = FALSE
  )
  if (!identical(resp$status, "Succeeded")) {
    stop("LFPS request did not succeed: ", resp$status)
  }
}

files <- utils::unzip(zip_path, exdir = raw_dir)
src <- terra::rast(grep("[.]tif$", files, value = TRUE))

# Warp onto the landcover grid ------------------------------------------------
elevation <- warp_to_template(
  src[["LF2020_Elev_CONUS"]],
  template,
  file.path(out_dir, "elevation.tif"),
  method = "bilinear",
  datatype = "FLT4S"
)
aspect <- warp_to_template(
  src[["LF2020_Asp_CONUS"]],
  template,
  file.path(out_dir, "aspect.tif"),
  method = "near",
  datatype = "INT2S"
)

# gdalwarp fills any cell the download does not cover with 0, a valid-looking
# elevation and a due-north aspect, so make sure the download covered the
# whole grid: elevation here is never below ~190 m
elev_min <- terra::global(elevation, "min", na.rm = TRUE)[1, 1]
if (elev_min < 100) {
  stop(sprintf(
    "elevation.tif has a minimum of %g m: the download does not cover the grid",
    elev_min
  ))
}
cat(sprintf(
  "Wrote %s and %s (elevation %g-%g m; %.1f%% of cells flat)\n",
  file.path(out_dir, "elevation.tif"),
  file.path(out_dir, "aspect.tif"),
  elev_min,
  terra::global(elevation, "max", na.rm = TRUE)[1, 1],
  100 * terra::global(aspect == -1, "mean", na.rm = TRUE)[1, 1]
))
