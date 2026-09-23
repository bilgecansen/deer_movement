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

#' The annual raster year for a deer's season
#'
#' The annual landcover and LANDFIRE rasters are built per deer-year, which
#' starts at the rut: file `<y>` holds breeding (Oct 15 y to Dec 14 y) and
#' non_breeding (Dec 15 y to Apr 30 y+1), then fawning (May 1 y+1 to Jun 30
#' y+1) and post_fawning (Jul 1 y+1 to Oct 14 y+1). Its crop classes come from
#' CDL y for the first two bands and CDL y+1 for the last two.
#'
#' A deer's `year` field is the calendar year of its steps, winter deer being
#' labelled by the December they start in. Autumn and winter deer therefore sit
#' in deer-year `year`, while spring and summer deer sit in `year - 1`: a deer
#' feeding in May 2017 walked through deer-year 2016's fawning season.
#'
#' @param year Deer's year field (calendar)
#' @param season Deer season code: one of "br", "nb", "fa", "pf"
#' @return Integer deer-year, i.e. the annual raster to read
deer_year <- function(year, season) {
  if (!season %in% c("br", "nb", "fa", "pf")) {
    stop(sprintf("Unknown season '%s' (expected one of br/nb/fa/pf)", season))
  }
  as.integer(year) - as.integer(season %in% c("fa", "pf"))
}

#' Load a season-specific annual landcover stack
#'
#' Reads the landcover file for the deer-year matching `year` and `season`
#' (see deer_year), selects the band for that `season`, and returns a
#' SpatRaster carrying:
#'   * `wiscland` — the categorical landcover band (factor; drives the amt
#'     `wiscland_*` design columns at fit time)
#'   * one binary indicator layer per non-reference predictor class, each named
#'     exactly like the class so the redistribution kernel can read
#'     `<class>_end` covariates at simulation / scoring time.
#'   * `forest_edge` — signed distance (m) to the nearest forest edge, from
#'     this same band (see forest_edge_distance)
#'
#' `open_water` is the water/exclusion class and `forest` is the reference
#' level; neither gets an indicator layer.
#'
#' @param year   Deer's year field (calendar); the file read is the deer-year
#' @param season Deer season code: one of "br", "nb", "fa", "pf"
#' @param folder Folder of annual landcover rasters
#' @param ref    Reference (intercept) class — emitted with no indicator layer
#' @return SpatRaster: `wiscland` + one binary layer per non-reference class +
#'   `forest_edge`
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
  fname <- sprintf("landcover_%d.tif", deer_year(year, season))
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

  do.call(c, c(list(lc), bins, list(forest_edge_distance(lc))))
}

#' Signed distance to the nearest forest edge
#'
#' Forest is FOREST_CLASSES (forest + wetland_forested); every other class,
#' open water included, is non-forest. Each cell gets the distance (m) from its
#' centre to the nearest cell centre of the other kind, positive inside forest
#' and negative outside. Half a cell is then taken off, so the edge itself,
#' which lies between two cells, sits at 0 and the cells either side of it read
#' +15 and -15.
#'
#' The distances are exact Euclidean, from a nearest-neighbour search (FNN)
#' over cell centres. terra::distance() is not used: it overstates ~0.1% of
#' them, by up to ~8 m, when the nearest cell touches its own kind only
#' diagonally. The search only needs edge cells as candidates, because the
#' nearest cell of the other kind always has a neighbour, diagonal included,
#' of the searching cell's kind: otherwise its neighbour towards that cell
#' would be closer and of the same kind.
#'
#' Computed over the whole band, never a crop, because a crop would hide any
#' nearest edge beyond its border. Cells within a few km of the map boundary
#' can still have their nearest edge off the map, which overstates their
#' distance; every deer endpoint lies at least ~4 km inside the boundary.
#'
#' @param lc Categorical landcover band (the `wiscland` layer)
#' @return SpatRaster with a single layer named "forest_edge"
forest_edge_distance <- function(lc) {
  forest <- Reduce(`|`, lapply(FOREST_CLASSES, function(cl) lc == cl))
  f <- terra::values(forest, mat = FALSE)
  is_f <- f %in% c(TRUE, 1)
  is_nf <- f %in% c(FALSE, 0)

  # Edge cells: any 8-neighbour of the other kind
  nb_max <- terra::values(
    terra::focal(forest, w = 3, fun = "max", na.rm = TRUE),
    mat = FALSE
  )
  nb_min <- terra::values(
    terra::focal(forest, w = 3, fun = "min", na.rm = TRUE),
    mat = FALSE
  )
  edge_nf <- is_nf & nb_max %in% 1
  edge_f <- is_f & nb_min %in% 0

  xy <- terra::xyFromCell(lc, seq_len(terra::ncell(lc)))
  nearest <- function(pool, from) {
    nn <- FNN::get.knnx(
      xy[pool, , drop = FALSE],
      xy[from, , drop = FALSE],
      k = 1
    )
    nn$nn.dist[, 1]
  }
  half <- terra::res(lc)[1] / 2
  d <- rep(NA_real_, terra::ncell(lc))
  if (any(is_f) && any(edge_nf)) {
    d[is_f] <- nearest(edge_nf, is_f) - half
  }
  if (any(is_nf) && any(edge_f)) {
    d[is_nf] <- half - nearest(edge_f, is_nf)
  }

  out <- terra::setValues(terra::rast(lc), d)
  names(out) <- "forest_edge"
  out
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

#' Representative timestamp for each monthly NDVI composite: the midpoint of
#' the month.
#'
#' The two sides of the pipeline select an NDVI layer by different rules, and
#' they have to agree:
#'   * fit time — extract_step_variables() uses
#'     amt::extract_covariates_var_time with when = "any", i.e. the layer
#'     nearest in time to the step;
#'   * kernel time — simulate_movement() / onestep_logscore*() index the stack
#'     by lubridate::month(t1_), i.e. the step's calendar month.
#'
#' With each layer stamped at its month's midpoint, the nearest layer is the
#' step's own month. The match cannot be exact at every hour, because months
#' differ in length and the stamps are UTC while the tracks are local time, so
#' a few steps on the first or last day of a month can still pick the
#' neighbouring layer. check K5b and check_wrangle.R W3.7 assert that any
#' disagreement stays on those days.
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

#' Load a year's LANDFIRE-derived covariates as model-ready layers
#'
#' Reads the `year` band of library/landfire/oak_mast.tif and oak_dist.tif and
#' the five scores in famd_<year>.tif, all on the landcover grid, and returns:
#'   * `oak_mast`      — 1 on oak mast habitat, 0 everywhere else
#'   * `oak_dist`      — distance (m) to the nearest oak mast cell; 0 inside it
#'   * `famd1`..`famd5` — FAMD scores of vegetation structure, NA where the
#'                        source has none
#'
#' The FAMD scores exist only on forested cells and stay NA elsewhere. They are
#' ordination coordinates, like PCA scores, so no fill value means "not
#' forest": 0 is the centre of the score space, an average forest cell. A fit
#' that names a FAMD term therefore drops every row whose endpoint is off the
#' mask (check_wrangle.R W4.3 and W6 count them). The FAMD mask is close to,
#' but not the same as, forest + wetland_forested in the landcover bands: it
#' differs in the fawning bands, and it leaves out ~1,600 forest cells whose
#' LANDFIRE vegetation type has no label (2016-2019).
#'
#' Bands are keyed by deer-year, like the landcover files (see deer_year).
#'
#' @param year   Deer's year field (calendar); the band read is the deer-year
#' @param season Deer season code: one of "br", "nb", "fa", "pf"
#' @param folder Folder of LANDFIRE rasters
#' @return SpatRaster: oak_mast, oak_dist, famd1..famd5
load_landfire <- function(year, season, folder = "library/landfire") {
  dyear <- deer_year(year, season)
  band <- as.character(dyear)
  famd_file <- file.path(folder, sprintf("famd_%d.tif", dyear))
  if (!file.exists(famd_file)) {
    stop(sprintf("Missing LANDFIRE file for %d: %s", dyear, famd_file))
  }
  mast <- terra::rast(file.path(folder, "oak_mast.tif"))
  dist <- terra::rast(file.path(folder, "oak_dist.tif"))
  if (!band %in% names(mast) || !band %in% names(dist)) {
    stop(sprintf("No %s band in the LANDFIRE oak rasters in %s", band, folder))
  }

  # The source mask is 1 on oak mast cells and NA everywhere else
  oak_mast <- terra::subst(mast[[band]], NA, 0)
  names(oak_mast) <- "oak_mast"
  oak_dist <- dist[[band]]
  names(oak_dist) <- "oak_dist"

  famd <- terra::rast(famd_file)
  # Source bands are dim1..dim5; a bare dim1_end column would say nothing
  names(famd) <- sub("^dim", "famd", names(famd))

  c(oak_mast, oak_dist, famd)
}

#' Load LANDFIRE elevation, northness and eastness
#'
#' Reads library/landfire/elevation.tif (metres) and aspect.tif (whole degrees
#' clockwise from true north, 0-359; -1 = flat), both already on the landcover
#' grid (see scripts/prep_landfire_topo.R), and turns aspect into
#' northness = cos(aspect) and eastness = sin(aspect). Flat cells face no
#' direction, so both are 0 there. None of this changes by year or season, so
#' one call serves every deer.
#'
#' @param folder Folder of LANDFIRE rasters
#' @return SpatRaster with layers `elevation`, `northness` and `eastness`
load_topo <- function(folder = "library/landfire") {
  elevation <- terra::rast(file.path(folder, "elevation.tif"))
  names(elevation) <- "elevation"

  aspect <- terra::rast(file.path(folder, "aspect.tif"))
  rad <- aspect * pi / 180
  northness <- terra::ifel(aspect < 0, 0, cos(rad))
  names(northness) <- "northness"
  eastness <- terra::ifel(aspect < 0, 0, sin(rad))
  names(eastness) <- "eastness"

  c(elevation, northness, eastness)
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

