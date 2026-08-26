# Landscape / raster loading ---------------------------------------------------
#
# Readers that turn the annual rasters in library/ and the per-deer home-range
# rasters in data/HR into the SpatRasters the rest of the pipeline consumes.
# Loading only -- no warping (see helpers_warp.R) and no covariate extraction
# (see helpers_track.R).
#
# Part of the helper library split out of scripts/helper_functions.R, which
# now sources every file in this folder. Scripts keep sourcing that one
# aggregator, so nothing here needs to be sourced directly.

#' Load a season-specific annual landcover stack
#'
#' Reads library/landcover/landcover_<year>.tif, selects the band matching the
#' deer's `season`, and returns a SpatRaster carrying:
#'   * `wiscland` — the categorical landcover band (factor; drives the iSSF
#'     `wiscland_*` design columns at fit time)
#'   * one binary indicator layer per non-reference predictor class, each named
#'     exactly like the class so the redistribution kernel can read
#'     `<class>_end` covariates at simulation / scoring time.
#'
#' `open_water` is the water/exclusion class and `forest` is the reference
#' level; neither gets an indicator layer.
#'
#' @param year   Year (integer); selects landcover_<year>.tif
#' @param season Deer season code: one of "br", "nb", "fa", "pf"
#' @param folder Folder of annual landcover rasters
#' @param ref    Reference (intercept) class — emitted with no indicator layer
#' @return SpatRaster: `wiscland` + one binary layer per non-reference class
load_landcover <- function(
  year,
  season,
  folder = "library/landcover",
  ref = "forest"
) {
  band_for <- c(
    br = "breeding",
    nb = "non_breeding",
    fa = "fawning",
    pf = "post_fawning"
  )
  if (!season %in% names(band_for)) {
    stop(sprintf("Unknown season '%s' (expected one of br/nb/fa/pf)", season))
  }

  fname <- sprintf("landcover_%d.tif", year)
  lc <- terra::rast(file.path(folder, fname))[[band_for[[season]]]]
  names(lc) <- "wiscland"

  # One binary indicator per predictor class (everything except the reference
  # level and the water class). terra is lazy here, so these per-class layers
  # are only materialised over the window callers later crop to.
  pred_levels <- setdiff(LANDCOVER_LEVELS, ref)
  bins <- lapply(pred_levels, function(lv) {
    b <- terra::ifel(lc == lv, 1, 0)
    names(b) <- lv
    b
  })

  do.call(c, c(list(lc), bins))
}

#' Binary open-water mask from a landcover stack (1 = water, 0 = land)
#'
#' Built from the `open_water` class of the categorical `wiscland` band. Used to
#' drop random/available steps that land on water in make_random_pt_extraction.
#'
#' @param landcover SpatRaster from load_landcover (must carry `wiscland`)
#' @return SpatRaster with a single layer named "Water"
make_water_mask <- function(landcover) {
  w <- terra::ifel(landcover[["wiscland"]] == "open_water", 1, 0)
  names(w) <- "Water"
  w
}

#' Load a year's NDVI as a 12-layer monthly stack
#'
#' Reads the twelve library/ndvi/ndvi_<year>_<MM>.tif single-band files for
#' `year` and stacks them in calendar order (layer i = month i), with per-layer
#' `time()` set to the first of each month. Drop-in replacement for the old
#' data/NDVI_year/NDVI_<year>.tif 12-band rasters:
#'   * indexable by month integer (`ndvi[[mo]]`) for simulate_movement /
#'     onestep_logscore / plot_kernel_grid, and
#'   * time-stamped for amt::extract_covariates_var_time in
#'     extract_step_variables.
#'
#' @param year   Year (integer)
#' @param folder Folder of monthly NDVI rasters
#' @return SpatRaster with 12 layers (Jan..Dec), time() set
load_ndvi <- function(year, folder = "library/ndvi") {
  files <- file.path(folder, sprintf("ndvi_%d_%02d.tif", year, 1:12))
  missing <- !file.exists(files)
  if (any(missing)) {
    stop(sprintf(
      "Missing NDVI files for %d: %s",
      year,
      paste(basename(files[missing]), collapse = ", ")
    ))
  }
  r <- terra::rast(files)
  names(r) <- month.abb
  terra::time(r) <- ndvi_layer_times(year)
  r
}

#' Representative timestamp for each monthly NDVI composite: the MIDPOINT of the
#' month, not its first day.
#'
#' This matters because the two sides of the pipeline select an NDVI layer by
#' different rules, and they have to agree:
#'   * fit time — extract_step_variables() uses amt::extract_covariates_var_time
#'     with when = "any", i.e. the layer NEAREST IN TIME to the step;
#'   * kernel time — simulate_movement() / onestep_logscore*() index the stack
#'     by
#'     lubridate::month(t1_), i.e. the step's calendar month.
#'
#' Stamped on the 1st, "nearest in time" resolves to the FOLLOWING month for any
#' step past mid-month: a step on 20 May sits 12 days from 1 Jun but 19 days
#' from 1 May. Models were therefore fit on one NDVI layer and simulated /
#' scored on another for roughly half of all steps (confirmed: 0% mismatch on
#' days 1-15, 100% on days 16+; see scripts/checks/check_kernel_gam.R,
#' check K5).
#'
#' Stamping the midpoint makes the two rules agree on 364-365 days a year,
#' against 195 under the old stamps. Measured exhaustively over 2016-2022, the
#' residual is always and only Jan 31 and Mar 01: months differ in length, so
#' the nearest-layer boundary cannot sit exactly on the month boundary either
#' side of February. Both days fall in winter, where the models use landcover
#' rather than NDVI, so no NDVI-using model is affected. check K5b asserts the
#' residual stays confined to those two days.
#'
#' @param year Year (integer)
#' @return POSIXct vector of length 12, the midpoint of each month.
ndvi_layer_times <- function(year) {
  starts <- as.POSIXct(
    sprintf("%d-%02d-01 00:00:00", year, 1:12),
    tz = "UTC"
  )
  nexts <- c(
    starts[-1],
    as.POSIXct(sprintf("%d-01-01 00:00:00", year + 1), tz = "UTC")
  )
  starts + as.numeric(difftime(nexts, starts, units = "secs")) / 2
}

#' Load a deer's HR binary raster and align it to a template
#'
#' Reads HRbin_<id>_<season>_<year>.tif from hr_folder, resamples it onto the
#' template's grid (they share a grid since both were derived from env_raster),
#' and fills any cells outside the HR raster's extent with 0. Returns a
#' one-layer SpatRaster named "HR_bin" that can be added directly to an
#' env raster.
#'
#' @param id Deer ID
#' @param season Season string (e.g. "fa", "nb")
#' @param year Year (integer)
#' @param template SpatRaster defining the target extent and grid
#' @param hr_folder Folder containing HRbin_<id>_<season>_<year>.tif files
#' @return SpatRaster with a single layer named "HR_bin"
load_hr_raster <- function(id, season, year, template, hr_folder = "data/HR") {
  fname <- sprintf("HRbin_%s_%s_%d.tif", id, season, year)
  hr <- terra::rast(file.path(hr_folder, fname))
  hr_aligned <- terra::resample(hr, template, method = "near")
  hr_aligned <- terra::subst(hr_aligned, NA, 0)
  names(hr_aligned) <- "HR_bin"
  hr_aligned
}

#' Load per-deer signed distance-to-HR-edge raster
#'
#' Positive inside HR, negative outside. Resampled with bilinear interpolation
#' (the field is continuous). Cells falling outside the source HR distance
#' raster's extent get the minimum observed value — a conservative "far
#' outside" fill to avoid NA propagation in the linear predictor.
#'
#' @param id Deer ID
#' @param season Season string (e.g. "fa", "nb")
#' @param year Year (integer)
#' @param template SpatRaster defining the target extent and grid
#' @param hr_folder Folder containing HRedge_<id>_<season>_<year>.tif files
#' @return SpatRaster with a single layer named "HR_edge"
load_hr_edge_raster <- function(
  id,
  season,
  year,
  template,
  hr_folder = "data/HR"
) {
  fname <- sprintf("HRedge_%s_%s_%d.tif", id, season, year)
  hr_edge <- terra::rast(file.path(hr_folder, fname))
  hr_edge_aligned <- terra::resample(hr_edge, template, method = "bilinear")
  min_val <- terra::global(hr_edge_aligned, "min", na.rm = TRUE)[1, 1]
  hr_edge_aligned <- terra::subst(hr_edge_aligned, NA, min_val)
  names(hr_edge_aligned) <- "HR_edge"
  hr_edge_aligned
}

#' Load per-deer distance-to-HR-center raster
#'
#' Distance (in CRS 6610 units, i.e. metres) from each cell to the ctmm μ of
#' the fitted home-range model — this is the "where the deer lives" centroid,
#' not a geometric polygon centroid. Always non-negative. Resampled with
#' bilinear interpolation; cells outside the source raster's extent are filled
#' with the max observed value (a conservative "far from center" stand-in).
#'
#' @param id Deer ID
#' @param season Season string (e.g. "fa", "nb")
#' @param year Year (integer)
#' @param template SpatRaster defining the target extent and grid
#' @param hr_folder Folder containing HRcenter_<id>_<season>_<year>.tif files
#' @return SpatRaster with a single layer named "HR_center"
load_hr_center_raster <- function(
  id,
  season,
  year,
  template,
  hr_folder = "data/HR"
) {
  fname <- sprintf("HRcenter_%s_%s_%d.tif", id, season, year)
  hr_c <- terra::rast(file.path(hr_folder, fname))
  hr_c_aligned <- terra::resample(hr_c, template, method = "bilinear")
  max_val <- terra::global(hr_c_aligned, "max", na.rm = TRUE)[1, 1]
  hr_c_aligned <- terra::subst(hr_c_aligned, NA, max_val)
  names(hr_c_aligned) <- "HR_center"
  hr_c_aligned
}

