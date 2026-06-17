# Helper functions -------------------------------------------------------------

# Landcover factor levels, reference (intercept) level first. Sourced from the
# annual `library/landcover/landcover_<year>.tif` rasters, whose seasonal bands
# carry 11 categories. `open_water` is intentionally omitted: it is the water /
# step-exclusion mask (see make_water_mask), not a predictor. `forest` is the
# reference level, so it gets no indicator layer and no coefficient.
LANDCOVER_LEVELS <- c(
  "forest",
  "corn",
  "soybeans",
  "alfalfa_hay",
  "small_grains",
  "other_ag",
  "wetland_forested",
  "wetland_open",
  "grassland",
  "developed"
)

# Landcover classes on which NDVI is a biologically meaningful covariate: the
# open-canopy ag / grassland classes whose greenness tracks forage state. On
# forest, wetlands and developed land NDVI is noise (closed canopy / standing
# water / impervious), so the breeding NDVI smooths (fit_GAM models 4 & 5) are
# switched off there. prepare_gam_data() turns this set into the is_veg /
# ndvi_veg / veg_class columns the NDVI block consumes.
NDVI_VEG_CLASSES <- c(
  "corn",
  "soybeans",
  "alfalfa_hay",
  "small_grains",
  "other_ag",
  "grassland"
)

# Crop buffer (metres, CRS 6610) added around each deer's track when cropping
# rasters for random-step generation, covariate extraction, HR raster creation,
# simulation, and scoring. Sized to clear the longest single observed/simulated
# step (~2.3 km max) plus room for simulated-path drift before the HR_center /
# HR_edge terms pull the path back. Keep create_hr_rasters' buffer >= this so HR
# rasters fully cover the runtime crop window.
CROP_BUFFER_M <- 3000

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
  terra::time(r) <- as.Date(sprintf("%d-%02d-01", year, 1:12))
  r
}

#' Generate random steps with environmental covariates (single-deer)
#'
#' Operates on a single-row tibble. Generates `n_pts` water-free random
#' control steps per observed step and writes them to `data[[output_col]]`
#' as a list-column.
#'
#' @param data Single-row dataframe with movement steps
#' @param n_pts Number of random points per step
#' @param water Binary raster for water bodies
#' @param model Sampling design for the random points:
#'   * "issf" — step lengths / turning angles drawn from gamma and von Mises
#'     distributions fitted to the data (amt defaults). Movement covariates are
#'     then corrections to these tentative distributions (the standard iSSF).
#'     Also usable for GAMs with a *parametric* movement kernel.
#'   * "nonp" — points drawn uniformly over a disc of radius R = max observed
#'     step length (turning angle ~ Unif(-pi, pi), step length =
#'     sqrt(Unif(0, R^2))), the availability required for *non-parametric*
#'     movement smooths s(sl_) / s(ta_) (Klappstein et al. 2024, Section 3.1).
#'     (Also valid for parametric kernels, just less efficient.)
make_random_pt_extraction <- function(
  data,
  n_pts,
  water,
  model = c("issf", "nonp"),
  stp_col = "stp",
  output_col = "random.stp"
) {
  model <- match.arg(model)
  stopifnot(nrow(data) == 1)

  # Crop water raster to this deer's buffer
  data_step <- data[[stp_col]][[1]]
  crop_extent <- sf::st_buffer(
    sf::st_as_sf(data_step, coords = c('x1_', 'y1_'), crs = 6610),
    CROP_BUFFER_M
  )
  water_local <- terra::crop(water, crop_extent)

  # Start with buffer for water removal
  n_random <- ceiling(n_pts * 10)

  # Generate random steps. "issf" samples from fitted gamma / von Mises
  # distributions (amt defaults); "nonp" samples uniformly over a disc of
  # radius R = max observed step length (amt draws from the rand_* pools with
  # replacement, so a fixed pool gives good coverage at any step count).
  rs <- if (model == "issf") {
    amt::random_steps(data_step, n_control = n_random)
  } else {
    R <- max(data_step$sl_, na.rm = TRUE)
    amt::random_steps(
      data_step,
      n_control = n_random,
      rand_sl = sqrt(stats::runif(1e5, 0, R^2)),
      rand_ta = stats::runif(1e5, -pi, pi)
    )
  }

  random_pts <- rs |>
    amt::extract_covariates(water_local, where = "end")

  # Filter and select final random points
  res <- random_pts |>
    # Select random steps not on water
    dplyr::filter(case_ == FALSE, Water == 0) |>
    # sample n_pts for each step
    dplyr::slice_sample(n = n_pts, by = step_id_) |>
    # Combine with observed steps
    dplyr::bind_rows(random_pts |> dplyr::filter(case_ == TRUE)) |>
    dplyr::ungroup() |>
    dplyr::select(-Water)

  # Warning if there are steps with less than random n_pts (+ 1 is the obs step)
  final_counts <- table(res$step_id_)
  if (any(final_counts < (n_pts + 1))) {
    message(paste(
      "\nNote: Individual",
      data$id,
      "has steps with missing random points."
    ))
    warning(
      paste("Data for individual", data$id, "is incomplete."),
      immediate. = TRUE
    )
  }

  data[[output_col]] <- list(res)
  data
}

#' Extract environmental variables for each step (single-deer)
#'
#' Operates on a single-row tibble. Builds a cropped env stack
#' (env + NDVI + the three HR rasters keyed by id/season/year) and extracts
#' covariates at every step start and end of the deer's random-steps tibble
#' in `data[[random_col]]`.
#'
#' Produces design-matrix columns including:
#'   * env layers (wiscland, ele, east, east2, north, dist) at start and end
#'   * HR_bin_start/end, HR_edge_start/end (m), HR_center_start/end (km)
#'   * ndvi_start/end (time-matched)
#'
#' @param data Single-row dataframe; must contain id, season, year and the
#'   random-steps column named by `random_col`. The HR raster files
#'   `HRbin_<id>_<season>_<year>.tif`, `HRedge_<id>_<season>_<year>.tif`, and
#'   `HRcenter_<id>_<season>_<year>.tif` must exist in `hr_folder`.
#' @param env Landcover stack for the deer's year+season (from load_landcover)
#' @param ndvi 12-layer monthly NDVI stack for the deer's year (from load_ndvi)
#' @param hr_folder Folder containing per-deer HR rasters
extract_step_variables <- function(
  data,
  env,
  ndvi,
  hr_folder = "data/HR",
  random_col = "random.stp",
  output_col = "stp.var"
) {
  stopifnot(nrow(data) == 1)

  data_step <- data[[random_col]][[1]]

  crop_extent <- sf::st_buffer(
    sf::st_as_sf(data_step, coords = c('x1_', 'y1_'), crs = 6610),
    CROP_BUFFER_M
  )

  env_cropped <- terra::crop(env, crop_extent)
  ndvi_local <- terra::crop(ndvi, crop_extent)

  # Per-deer HR rasters, aligned to env_cropped via the load_hr_*_raster
  # helpers (which handle resample + NA fill)
  hr_bin <- load_hr_raster(
    data$id,
    data$season,
    data$year,
    env_cropped,
    hr_folder
  )
  hr_edge <- load_hr_edge_raster(
    data$id,
    data$season,
    data$year,
    env_cropped,
    hr_folder
  )
  hr_center <- load_hr_center_raster(
    data$id,
    data$season,
    data$year,
    env_cropped,
    hr_folder
  )

  # Convert to log scale
  hr_center_log <- log1p(hr_center)
  names(hr_center_log) <- "HR_center_log"

  lc_levels <- LANDCOVER_LEVELS

  data_ssf <- data_step |>
    amt::extract_covariates(env_cropped, where = 'both') |>
    amt::extract_covariates(hr_bin, where = 'both') |>
    amt::extract_covariates(hr_edge, where = 'both') |>
    amt::extract_covariates(hr_center, where = 'both') |>
    amt::extract_covariates(hr_center_log, where = 'both') |>
    amt::extract_covariates_var_time(
      ndvi_local,
      max_time = lubridate::days(31),
      when = "any",
      where = "both",
      name_covar = "ndvi"
    ) |>
    # Setting the intercept to the 'central.hardwoods' land type
    dplyr::mutate(
      wiscland_start = factor(wiscland_start, levels = lc_levels),
      wiscland_end = factor(wiscland_end, levels = lc_levels)
    ) |>
    amt::time_of_day(include.crepuscule = FALSE, where = "both") |>
    dplyr::mutate(
      tod_start_day = as.integer(tod_start_ == "day"),
      tod_start_night = as.integer(tod_start_ == "night"),
      # Continuous time of day at the step start (decimal hour, 0-24) for the
      # GAM cyclic spline s(tod_, bs = "cc"); complements the day/night factor.
      tod_ = lubridate::hour(t1_) + lubridate::minute(t1_) / 60,
      days = lubridate::yday(t2_) - min(lubridate::yday(t2_)) + 1
    )

  data[[output_col]] <- list(data_ssf)
  data
}

#' Load a deer's HR binary raster and align it to a template
#'
#' Reads HRbin_<id>_<season>_<year>.tif from hr_folder, resamples it onto the
#' template's grid (they share a grid since both were derived from env_raster),
#' and fills any cells outside the HR raster's extent with 0. Returns a
#' one-layer SpatRaster named "HR_bin" that can be added directly to an env
#' raster.
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

#' Warp a raster to match a template's CRS / extent / resolution
#'
#' Wraps gdalwarp via sf::gdal_utils. Same job as terra::project()/resample(),
#' but multi-threaded, streaming-capable, and with full control over output
#' datatype, nodata, compression, and tiling.
#'
#' Default datatype is INT2U — covers up to 65535, fits Wiscland-style
#' 4-digit codes and most WI categorical rasters in the smallest size.
#' Override with FLT4S for continuous data (NDVI, elevation, etc.).
#'
#' Before warping, the function checks the input's data range against the
#' chosen datatype's bounds; out-of-range values raise an error (rather than
#' being silently clipped by gdalwarp). It also warns when float-storage
#' input is being converted to an integer datatype, since fractional values
#' would be truncated.
#'
#' @param r           SpatRaster or path to input raster
#' @param template    SpatRaster or path defining target CRS / extent / resolution
#' @param outfile     Output file path (string)
#' @param method      GDAL resampling method ("near", "bilinear", "cubic", ...).
#'                    Default "near".
#' @param datatype    Terra-style datatype string. Options:
#'                    INT1U, INT2U, INT2S, INT4U, INT4S, FLT4S, FLT8S.
#'                    Default "INT2U".
#' @param nodata      Value treated as NA in input AND output. NA disables the
#'                    -srcnodata/-dstnodata flags. Default NA.
#' @param compress    GDAL compression scheme. Default "LZW".
#' @param threads     Worker threads for warping. Default "ALL_CPUS".
#' @param verify      If TRUE, run terra::compareGeom() against template after
#'                    warp and warn on misalignment. Default TRUE.
#' @param check_range If TRUE, scan input min/max and verify they fit the
#'                    target datatype's bounds before warping. Adds a one-time
#'                    minmax compute cost; turn off for very large rasters
#'                    you've already validated. Default TRUE.
#' @return SpatRaster wrapping the warped output file (lazy reference).
warp_to_template <- function(
  r,
  template,
  outfile,
  method = "near",
  datatype = "INT2U",
  nodata = NA,
  compress = "LZW",
  threads = "ALL_CPUS",
  verify = TRUE,
  check_range = TRUE
) {
  # Terra datatype -> GDAL -ot. The two must match or downstream stacks
  # will silently disagree on dtype.
  gdal_dt <- c(
    INT1U = "Byte",
    INT2U = "UInt16",
    INT2S = "Int16",
    INT4U = "UInt32",
    INT4S = "Int32",
    FLT4S = "Float32",
    FLT8S = "Float64"
  )
  if (!datatype %in% names(gdal_dt)) {
    stop(
      "datatype must be one of: ",
      paste(names(gdal_dt), collapse = ", ")
    )
  }
  ot <- gdal_dt[[datatype]]

  # Resolve input to a SpatRaster handle (for the range / type check below).
  if (inherits(r, "SpatRaster")) {
    r_in <- r
  } else if (is.character(r)) {
    r_in <- terra::rast(r)
  } else {
    stop("r must be a SpatRaster or a file path string")
  }

  # Datatype compatibility check — catches silent clipping that would
  # otherwise happen during gdalwarp / writeRaster.
  if (check_range) {
    dt_bounds <- list(
      INT1U = c(0, 255),
      INT2U = c(0, 65535),
      INT2S = c(-32768, 32767),
      INT4U = c(0, 4294967295),
      INT4S = c(-2147483648, 2147483647),
      FLT4S = c(-3.4e38, 3.4e38),
      FLT8S = c(-Inf, Inf)
    )
    bounds <- dt_bounds[[datatype]]

    mm <- terra::minmax(r_in, compute = TRUE)
    data_min <- suppressWarnings(min(mm[1, ], na.rm = TRUE))
    data_max <- suppressWarnings(max(mm[2, ], na.rm = TRUE))

    if (is.finite(data_min) && data_min < bounds[1]) {
      stop(sprintf(
        "Input min %g is below datatype %s lower bound %g — 
        would silently clip. Pick a wider datatype.",
        data_min,
        datatype,
        bounds[1]
      ))
    }
    if (is.finite(data_max) && data_max > bounds[2]) {
      stop(sprintf(
        "Input max %g is above datatype %s upper bound %g — 
        would silently clip. Pick a wider datatype.",
        data_max,
        datatype,
        bounds[2]
      ))
    }

    # Float-storage -> integer datatype: would truncate fractional values.
    # Storage check first (cheap); if storage is float, also confirm the
    # data actually has fractional values via the min/max we already have.
    if (datatype %in% c("INT1U", "INT2U", "INT2S", "INT4U", "INT4S")) {
      in_dtype <- terra::datatype(r_in)
      if (any(grepl("FLT", in_dtype))) {
        has_fractional <- any(mm != round(mm), na.rm = TRUE)
        if (has_fractional) {
          warning(sprintf(
            "Input has float storage (%s) with fractional values; 
            target datatype %s is integer — values will be truncated.",
            paste(unique(in_dtype), collapse = ","),
            datatype
          ))
        }
      }
    }
  }

  # Resolve template (path or SpatRaster) and pull its grid spec.
  if (is.character(template)) {
    template <- terra::rast(template)
  }
  tmpl_ext <- as.vector(terra::ext(template))
  tmpl_res <- terra::res(template)

  # Use EPSG:CODE form when available (compact, unambiguous); fall back to
  # full WKT otherwise. -t_srs accepts both.
  crs_info <- terra::crs(template, describe = TRUE)
  tmpl_crs_str <- if (!is.na(crs_info$code)) {
    paste0(crs_info$authority, ":", crs_info$code)
  } else {
    terra::crs(template)
  }

  # gdalwarp reads from disk, so a SpatRaster argument gets written to a temp
  # file first. Path inputs go through directly.
  cleanup_dir <- NULL
  if (inherits(r, "SpatRaster")) {
    cleanup_dir <- tempfile()
    dir.create(cleanup_dir)
    infile <- file.path(cleanup_dir, "in.tif")
    wopt <- list(
      datatype = datatype,
      gdal = sprintf("COMPRESS=%s", compress)
    )
    if (!is.na(nodata)) {
      wopt$NAflag <- nodata
    }
    terra::writeRaster(r, infile, overwrite = TRUE, wopt = wopt)
  } else {
    infile <- r
  }

  # -te + -tr force cell-for-cell alignment with the template (more
  # deterministic than letting gdalwarp infer the output grid).
  opts <- c(
    "-overwrite",
    "-t_srs",
    tmpl_crs_str,
    "-r",
    method,
    "-ot",
    ot,
    "-tr",
    as.character(tmpl_res[1]),
    as.character(tmpl_res[2]),
    "-te",
    as.character(tmpl_ext["xmin"]),
    as.character(tmpl_ext["ymin"]),
    as.character(tmpl_ext["xmax"]),
    as.character(tmpl_ext["ymax"]),
    "-wo",
    sprintf("NUM_THREADS=%s", threads),
    "-co",
    "TILED=YES",
    "-co",
    sprintf("COMPRESS=%s", compress),
    "-co",
    "BIGTIFF=IF_SAFER"
  )
  if (!is.na(nodata)) {
    opts <- c(
      opts,
      "-srcnodata",
      as.character(nodata),
      "-dstnodata",
      as.character(nodata)
    )
  }

  sf::gdal_utils(
    util = "warp",
    source = infile,
    destination = outfile,
    options = opts
  )

  # Copy RAT / band-metadata sidecar if present. gdalwarp doesn't move
  # .aux.xml, but pixel values are unchanged by reprojection so the original
  # RAT (and any layer-name / time metadata terra writes) stays valid.
  src_aux <- paste0(infile, ".aux.xml")
  out_aux <- paste0(outfile, ".aux.xml")
  if (file.exists(src_aux)) {
    file.copy(src_aux, out_aux, overwrite = TRUE)
  }

  # Clean up temp input
  if (!is.null(cleanup_dir)) {
    unlink(cleanup_dir, recursive = TRUE)
  }

  # Verify alignment
  warped <- terra::rast(outfile)
  if (verify && !terra::compareGeom(template, warped, stopOnError = FALSE)) {
    warning(
      "Output does NOT match template grid; downstream stacking will fail."
    )
  }

  warped
}

#' Fit a single iSSF model for one individual
#' @param ssf_data Step selection data for one deer (from stp.var.train)
#' @param formula Model formula string
#' @return List with iss (fitted model), coeff (tidy coefficients), aic (AIC value)
fit_mod <- function(ssf_data, formula) {
  iss <- tryCatch(
    ssf_data |> amt::fit_issf(as.formula(formula), model = TRUE),
    error = function(err) "Error"
  )

  coeff <- if (is.character(iss)) NA else broom::tidy(iss$model)
  aic <- if (is.character(iss)) NA_real_ else AIC(iss$model)

  list(iss = iss, coeff = coeff, aic = aic)
}

#' Reshape iSSF step data into mgcv Cox-PH (GAM-SSF) form
#'
#' Following Klappstein et al. (2024), an (integrated) SSF is fit as a
#' stratified Cox proportional-hazards model in mgcv. This adds the three
#' columns that the cox.ph fit needs, leaving every existing covariate column
#' (sl_, ta_, wiscland_end, ndvi_end, HR_*, tod_*, ...) untouched:
#'   * times   — constant "event time" (all strata share it; the within-stratum
#'               tie of one event + N censored points reproduces conditional
#'               logistic regression).
#'   * stratum — integer stratification index (one per observed step + its
#'               random points), from step_id_.
#'   * obs     — used/available indicator (1 = observed step, 0 = random),
#'               passed to gam() as `weights` (the cox.ph censoring indicator).
#'
#' @param ssf_data Step-selection data for one deer (a deer's stp.var tibble)
#' @return The same data with `times`, `stratum`, `obs` columns added
prepare_gam_data <- function(ssf_data) {
  ssf_data$times <- 1
  ssf_data$stratum <- as.integer(factor(ssf_data$step_id_))
  ssf_data$obs <- as.integer(ssf_data$case_)

  # NDVI block support for breeding models 4 & 5 (see make_formulas in
  # fit_GAM.R). NDVI is meaningful only on NDVI_VEG_CLASSES:
  #   * is_veg    — 1/0 indicator passed as `by=` to the NDVI smooths so they are
  #                 exactly zero on non-veg steps (and on steps lacking NDVI).
  #   * ndvi_veg  — NDVI kept as-is wherever it exists (so the smooth basis is not
  #                 skewed by a pile of placeholders); only NA is filled with an
  #                 in-range stand-in, which is_veg zeroes out anyway. This also
  #                 stops na.omit from dropping NDVI-less steps the movement /
  #                 landcover terms still need.
  #   * veg_class — factor with only the veg levels; non-veg / NDVI-less steps are
  #                 parked on the first level but made inert by is_veg = 0, so they
  #                 never enter that class's per-class curve.
  has_ndvi <- ssf_data$wiscland_end %in% NDVI_VEG_CLASSES &
    !is.na(ssf_data$ndvi_end)
  ssf_data$is_veg <- as.numeric(has_ndvi)
  ssf_data$ndvi_veg <- ifelse(
    is.na(ssf_data$ndvi_end),
    stats::median(ssf_data$ndvi_end, na.rm = TRUE),
    ssf_data$ndvi_end
  )
  ssf_data$veg_class <- factor(
    ifelse(has_ndvi, as.character(ssf_data$wiscland_end), NDVI_VEG_CLASSES[1]),
    levels = NDVI_VEG_CLASSES
  )

  ssf_data
}

#' Fit a single GAM-SSF model for one individual (mgcv cox.ph)
#'
#' GAM analogue of fit_mod(). `formula` is the RHS only (e.g.
#' "s(sl_) + s(ta_, bs = 'cc') + s(ndvi_end)"); the
#' "cbind(times, stratum) ~" Cox-PH response is prepended here, and the
#' used/available indicator `obs` is passed as the cox.ph censoring weight.
#' Fitting is via REML, as recommended in Klappstein et al. (2024).
#'
#' The cyclic time-of-day spline s(tod_, bs = "cc") is given explicit knots at
#' 0 and 24 so the curve wraps at midnight; the knot is ignored for formulas
#' without a tod_ smooth.
#'
#' `select = TRUE` adds the Marra & Wood (2011) double penalty, which also
#' penalises each smooth's null space — so an unsupported smooth can shrink past
#' linear all the way to zero (i.e. be dropped). gam_smooth_diag() records when
#' that happens (status "near-linear" / "removed") and whether any smooth is
#' pressed against its basis ceiling (status "k-bound" -> raise k).
#'
#' `drop.unused.levels = FALSE` retains every declared landcover level even when
#' a deer never visited it, so the saved model can be predicted on any class at
#' simulation / scoring time (unvisited classes fall back to the population
#' average; see the inline note on the gam() call).
#'
#' @param gam_data Output of prepare_gam_data() for one deer
#' @param formula  RHS-only model formula string
#' @param select   Use the double-penalty shrinkage (default TRUE)
#' @return List with gam (fitted model or "Error"), coeff (tidy parametric
#'   terms), aic (AIC value), and smooth_diag (per-smooth edf / k' audit, or
#'   NULL if the fit failed)
fit_gam_mod <- function(gam_data, formula, select = TRUE) {
  full_formula <- stats::as.formula(
    paste("cbind(times, stratum) ~", formula)
  )

  gfit <- tryCatch(
    mgcv::gam(
      full_formula,
      data = gam_data,
      family = mgcv::cox.ph(),
      weights = obs,
      method = "REML",
      select = select,
      knots = list(tod_ = c(0, 24)),
      # Keep ALL declared factor levels (LANDCOVER_LEVELS / NDVI_VEG_CLASSES),
      # not just the ones this deer visited. Classes with no training steps get a
      # coefficient pinned at ~0 in the re / fs smooths -- which is exactly the
      # population-average response -- so prediction (simulation / scoring) on a
      # landcover class the deer never visited works with a plain predict():
      # the re and per-class fs fall back to 0 and only the shared smooths (e.g.
      # s(ndvi_veg, by = is_veg)) carry the effect. Verified to leave the fits
      # for visited classes (coefficients and variance components) unchanged.
      drop.unused.levels = FALSE
    ),
    error = function(err) "Error"
  )

  coeff <- if (is.character(gfit)) {
    NA
  } else {
    tryCatch(broom::tidy(gfit, parametric = TRUE), error = function(e) NA)
  }
  aic <- if (is.character(gfit)) NA_real_ else AIC(gfit)
  smooth_diag <- if (is.character(gfit)) NULL else gam_smooth_diag(gfit)

  list(gam = gfit, coeff = coeff, aic = aic, smooth_diag = smooth_diag)
}

#' Per-smooth basis / shrinkage diagnostic for a fitted GAM-SSF
#'
#' Reads each smooth's effective degrees of freedom (edf) and basis ceiling (k')
#' from mgcv::k.check() and labels its status. We deliberately ignore k.check()'s
#' k-index and p-value columns: those are residual-based and unreliable for a
#' cox.ph SSF (the same reason Klappstein et al. 2024 advise against PH
#' residuals). The trustworthy signal is edf vs k':
#'   * "removed"     edf <= removed_edf  — shrunk to ~zero (only with select)
#'   * "near-linear" edf <= linear_edf   — collapsed to the penalty null space
#'                                         (a straight line, or a constant for
#'                                         cyclic terms = no time-of-day effect)
#'   * "k-bound"     edf >= k_bound_ratio * k' — pressed against the ceiling;
#'                                         consider a higher k
#'   * "ok"          otherwise
#'
#' @param gfit Fitted mgcv gam
#' @param k_bound_ratio edf/k' at/above which a smooth is flagged "k-bound"
#' @param linear_edf edf at/below which a smooth is flagged "near-linear"
#' @param removed_edf edf at/below which a smooth is flagged "removed"
#' @return data.frame(smooth, k_prime, edf, edf_ratio, status), or NULL
gam_smooth_diag <- function(
  gfit,
  k_bound_ratio = 0.8,
  linear_edf = 1.5,
  removed_edf = 0.5
) {
  kc <- tryCatch(suppressWarnings(mgcv::k.check(gfit)), error = function(e) {
    NULL
  })
  if (is.null(kc) || nrow(kc) == 0) {
    return(NULL)
  }
  df <- data.frame(
    smooth = rownames(kc),
    k_prime = kc[, "k'"],
    edf = round(kc[, "edf"], 3),
    row.names = NULL
  )
  df$edf_ratio <- round(df$edf / df$k_prime, 3)
  df$status <- ifelse(
    df$edf <= removed_edf,
    "removed",
    ifelse(
      df$edf <= linear_edf,
      "near-linear",
      ifelse(df$edf_ratio >= k_bound_ratio, "k-bound", "ok")
    )
  )
  df
}

# GAM redistribution kernel ----------------------------------------------------
# amt::redistribution_kernel() only accepts fit_clogit (iSSF) objects: it reads
# fixed coefficients via coef() and a fixed model.matrix. Our movement models are
# mgcv cox.ph GAMs whose effects are penalised smooths evaluated through
# predict.gam(). The functions below reproduce amt's discrete-landscape kernel
# (kernel_setup + ssf_weights), swapping only the `model.matrix %*% coefs` step
# for `predict(gam, type = "link")`. The disc geometry, the tentative-kernel
# compensation (phi - log(sl_)), the exp() stabilisation, and the normalise-to-
# sum-1 step all match amt, so the result is interchangeable with the iSSF path
# (Klappstein et al. 2024; derivation in docs/kernel_exponent.html).
#
# Parametric design only: the GAMs are fit on the gamma / von Mises steps
# (stp.var), so the tentative distributions h ride along on that tibble as
# attr(stp.var, "sl_") / attr(stp.var, "ta_") and are passed in as sl_distr /
# ta_distr. Compensation re-adds h and the polar->planar Jacobian, giving the
# cell weight  h(sl_, ta_) * exp(eta) / sl_  before normalisation.

#' Log tentative movement kernel at each candidate cell (GAM)
#'
#' Reimplements amt's internal movement_kernel1() so we don't depend on an
#' unexported function. Returns the log density of the tentative step-length and
#' turning-angle distributions (the cloud the random points were sampled from),
#' evaluated per cell. Only gamma / exponential step lengths and a von Mises
#' turning angle are supported (the families amt::random_steps produces).
#'
#' @param xy       Candidate-step tibble carrying sl_ and ta_ columns
#' @param sl_distr Tentative step-length distribution (amt sl_distr, e.g.
#'                 attr(stp.var, "sl_"))
#' @param ta_distr Tentative turning-angle distribution (amt ta_distr) or NULL
#' @return numeric vector of log tentative-kernel values, one per row of xy
gam_movement_kernel <- function(xy, sl_distr, ta_distr = NULL) {
  phi <- switch(
    sl_distr$name,
    gamma = -(1 / sl_distr$params$scale) *
      xy$sl_ +
      log(xy$sl_) * (sl_distr$params$shape - 1),
    exp = -xy$sl_ * sl_distr$params$rate,
    stop("Unsupported step-length distribution: ", sl_distr$name)
  )
  if (!is.null(ta_distr) && ta_distr$name == "vonmises") {
    phi <- phi + cos(xy$ta_) * ta_distr$params$kappa
  }
  phi
}

#' Default covariate-extraction callback for a GAM redistribution kernel
#'
#' GAM analogue of the cov_fun used by simulate_movement(). Extracts every map
#' layer at the step start and end and adds the continuous time-of-day covariate
#' tod_ (decimal hour at the step start) that the cyclic spline s(tod_, bs = "cc")
#' needs. wiscland_start/end are re-levelled to LANDCOVER_LEVELS so predict.gam
#' sees the same factor levels as the fit (cells whose class is outside those
#' levels, e.g. open_water, become NA and so get weight 0).
#'
#' @param xy  Candidate-step tibble (from the kernel grid) with t1_ set
#' @param map SpatRaster carrying all covariate layers the GAM references
#' @return xy with covariate columns (and tod_) added
gam_cov_fun <- function(xy, map) {
  xy |>
    amt::extract_covariates(map, where = "both") |>
    dplyr::mutate(
      tod_ = lubridate::hour(t1_) + lubridate::minute(t1_) / 60,
      wiscland_start = factor(wiscland_start, levels = LANDCOVER_LEVELS),
      wiscland_end = factor(wiscland_end, levels = LANDCOVER_LEVELS)
    )
}

#' Redistribution-kernel weights for a fitted GAM (analogue of ssf_weights)
#'
#' Computes the unnormalised kernel weight at every candidate cell:
#'   w = exp( eta + phi - log(sl_) )
#' where eta = predict(gam, type = "link") is the cox.ph linear predictor (the
#' fitted movement corrections + habitat), phi is the log tentative movement
#' kernel (gam_movement_kernel), and -log(sl_) is the polar->planar Jacobian.
#' With compensate.movement = FALSE the phi - log(sl_) term is dropped (uniform /
#' non-parametric availability), leaving w = exp(eta). Exponentiation uses the
#' same mean-centring stabilisation as amt; non-finite weights become 0.
#'
#' @param xy        Candidate-step tibble after covariate extraction
#' @param gam_fit   Fitted mgcv gam (family = cox.ph)
#' @param sl_distr  Tentative step-length distribution (required if compensating)
#' @param ta_distr  Tentative turning-angle distribution or NULL
#' @param compensate.movement Add the tentative-kernel + Jacobian terms (TRUE for
#'   the parametric gamma / von Mises design)
#' @return numeric vector of non-negative weights, one per row of xy
gam_kernel_weights <- function(
  xy,
  gam_fit,
  sl_distr = NULL,
  ta_distr = NULL,
  compensate.movement = TRUE
) {
  w <- as.numeric(stats::predict(
    gam_fit,
    newdata = xy,
    type = "link",
    na.action = stats::na.pass
  ))

  if (compensate.movement) {
    if (is.null(sl_distr)) {
      stop(
        "compensate.movement = TRUE needs sl_distr (the tentative step-length distribution)."
      )
    }
    phi <- gam_movement_kernel(xy, sl_distr, ta_distr)
    w <- w + phi - log(xy$sl_)
  }

  w <- exp(w - mean(w[is.finite(w)], na.rm = TRUE))
  w[!is.finite(w)] <- 0
  w
}

#' Build a redistribution kernel from a fitted GAM
#'
#' GAM analogue of amt::redistribution_kernel() for a discrete landscape. Lays a
#' candidate point at every map cell within `max.dist` of the start, computes
#' each cell's step length / turning angle relative to the start, runs `fun` to
#' attach covariates, scores the cells with gam_kernel_weights(), and returns
#' either the normalised kernel raster (`as.rast = TRUE`; each cell is the
#' probability the next step lands there) or a sampled endpoint (`as.rast =
#' FALSE`; used by path simulation). The geometry mirrors amt's kernel_setup.
#'
#' @param x         Fitted mgcv gam (family = cox.ph)
#' @param map       SpatRaster carrying every covariate layer the GAM references
#' @param start     amt sim_start (from amt::make_start): x_, y_, t_, ta_, dt
#' @param fun       Covariate-extraction callback fun(xy, map); default gam_cov_fun
#' @param sl_distr  Tentative step-length distribution (attr(stp.var, "sl_"))
#' @param ta_distr  Tentative turning-angle distribution (attr(stp.var, "ta_"))
#' @param max.dist  Disc radius (m); default = 0.99 quantile of sl_distr (matches
#'   amt::get_max_dist)
#' @param compensate.movement Re-add the tentative kernel + Jacobian (TRUE for the
#'   parametric design)
#' @param normalize Scale the raster to sum to 1 (as.rast only)
#' @param as.rast   TRUE -> normalised SpatRaster; FALSE -> sampled endpoint(s)
#' @param n.sample  Endpoints to draw when as.rast = FALSE
#' @param tolerance.outside Max fraction of candidate cells allowed beyond the map
#'   extent before bailing out (returns NULL)
#' @return list(args, redistribution.kernel) of class
#'   "redistribution_kernel_gam", or NULL if the kernel ran off the map / has no
#'   mass
redistribution_kernel_gam <- function(
  x,
  map,
  start,
  fun = gam_cov_fun,
  sl_distr = NULL,
  ta_distr = NULL,
  max.dist = NULL,
  compensate.movement = TRUE,
  normalize = TRUE,
  as.rast = TRUE,
  n.sample = 1,
  tolerance.outside = 0.01
) {
  arguments <- as.list(environment())

  if (compensate.movement && is.null(sl_distr)) {
    stop(
      "compensate.movement = TRUE needs sl_distr (the tentative step-length distribution)."
    )
  }

  # Default disc radius: 0.99 quantile of the tentative step length (matches
  # amt::get_max_dist.fit_clogit).
  if (is.null(max.dist)) {
    if (is.null(sl_distr)) {
      stop(
        "Provide max.dist, or sl_distr so it can default to the 0.99 step-length quantile."
      )
    }
    max.dist <- ceiling(do.call(
      paste0("q", sl_distr$name),
      c(list(p = 0.99), sl_distr$params)
    ))
  }

  # --- Candidate grid: disc of radius max.dist around the start (kernel_setup) -
  pt <- sf::st_sf(
    geom = sf::st_sfc(sf::st_point(
      as.numeric(start[1, c("x_", "y_")])
    ))
  )
  disc <- sf::st_buffer(pt, dist = max.dist)
  r1 <- terra::rasterize(terra::vect(disc), terra::crop(map, disc))
  xy_crds <- terra::crds(r1)

  k <- tibble::tibble(x2_ = xy_crds[, 1], y2_ = xy_crds[, 2])
  ta <- atan2(k$y2_ - start$y_[1], k$x2_ - start$x_[1])
  ta <- ifelse(ta < 0, 2 * pi + ta, ta)
  ta <- (ta - start$ta_[1]) %% (2 * pi)
  ta <- ifelse(ta > pi, (2 * pi - ta) * -1, ta)
  k$ta_ <- ta
  k$sl_ <- sqrt((k$x2_ - start$x_[1])^2 + (k$y2_ - start$y_[1])^2)
  k$x1_ <- start$x_[1]
  k$y1_ <- start$y_[1]
  class(k) <- c("steps_xyt", "steps_xy", class(k))
  attr(k, "crs") <- attr(start, "crs")

  # Bail out if too much of the kernel sits beyond the map (mirror amt)
  bb <- as.vector(terra::ext(map))
  frac_out <- mean(
    k$x2_ < bb["xmin"] |
      k$x2_ > bb["xmax"] |
      k$y2_ < bb["ymin"] |
      k$y2_ > bb["ymax"]
  )
  if (frac_out > tolerance.outside) {
    warning(sprintf(
      "%.2f%% of candidate steps end outside the map (tolerance %.2f%%); returning NULL.",
      100 * frac_out,
      100 * tolerance.outside
    ))
    return(NULL)
  }

  # Times + covariates
  k$t1_ <- start$t_
  k$t2_ <- start$t_ + start$dt
  k <- fun(k, map)

  # Weights
  w <- gam_kernel_weights(k, x, sl_distr, ta_distr, compensate.movement)

  if (!any(w > 0)) {
    warning("Redistribution kernel has no positive mass; returning NULL.")
    return(NULL)
  }

  rk <- if (!as.rast) {
    idx <- sample.int(nrow(k), size = n.sample, prob = w)
    dplyr::select(k[idx, ], x_ = x2_, y_ = y2_, t2_)
  } else {
    r <- terra::rast(data.frame(k[, c("x2_", "y2_")], w = w))
    if (normalize) {
      r <- r / terra::global(r, "sum", na.rm = TRUE)[1, 1]
    }
    names(r) <- "kernel"
    r
  }

  res <- list(args = arguments, redistribution.kernel = rk)
  class(res) <- c("redistribution_kernel_gam", "list")
  res
}

#' Simulate a movement path from a GAM redistribution kernel
#'
#' Standalone GAM analogue of amt::simulate_path(). Given a kernel built by
#' redistribution_kernel_gam(..., as.rast = FALSE), it walks `n.steps` forward,
#' rebuilding the kernel at each new position + bearing and drawing one
#' endpoint. The per-step rebuild is essential: the kernel depends on the
#' current location (covariates change) and the previous turning angle (the von
#' Mises correction is relative to the last heading), so a path is *not* a
#' repeated draw from one static kernel.
#'
#' I/O mirrors amt::simulate_path exactly, so the two are interchangeable inside
#' simulate_movement(): it returns an (n.steps + 1)-row tibble (x_, y_, t_, dt)
#' whose first row is the start, with the time grid spaced by the start's dt. If
#' a rebuilt kernel runs off the map or has no mass (redistribution_kernel_gam
#' returns NULL), the walk stops early and the path so far is returned (trailing
#' rows stay NA) — again matching amt.
#'
#' @param kernel  A redistribution_kernel_gam object (built with as.rast =
#'   FALSE). Its `args` carry everything needed to rebuild the kernel each step
#'   (x, map, fun, sl_distr, ta_distr, compensate.movement, ...).
#' @param n.steps Number of steps to simulate
#' @param start   sim_start to begin from (default: the kernel's own start)
#' @param verbose Warn when the walk stops early (default FALSE)
#' @return tibble(x_, y_, t_, dt) with n.steps + 1 rows
simulate_path_gam <- function(
  kernel,
  n.steps = 100,
  start = kernel$args$start,
  verbose = FALSE
) {
  mod <- kernel$args

  xy <- tibble::tibble(
    x_ = rep(NA_real_, n.steps + 1),
    y_ = NA_real_,
    t_ = start$t_ + start$dt * (0:n.steps),
    dt = start$dt
  )
  xy$x_[1] <- start$x_
  xy$y_[1] <- start$y_

  for (i in 1:n.steps) {
    # Rebuild the kernel at the current (position, bearing) and draw one point.
    rk <- redistribution_kernel_gam(
      x = mod$x,
      map = mod$map,
      start = start,
      fun = mod$fun,
      sl_distr = mod$sl_distr,
      ta_distr = mod$ta_distr,
      max.dist = mod$max.dist,
      compensate.movement = mod$compensate.movement,
      normalize = TRUE,
      as.rast = FALSE,
      n.sample = 1,
      tolerance.outside = mod$tolerance.outside
    )

    if (is.null(rk)) {
      if (verbose) {
        warning(sprintf(
          "Simulation stopped after %d steps: kernel ran off the map or had no mass.",
          i - 1
        ))
      }
      return(xy)
    }

    rk <- rk$redistribution.kernel
    new.ta <- atan2(rk$y_[1] - start$y_[1], rk$x_[1] - start$x_[1])
    xy$x_[i + 1] <- rk$x_[1]
    xy$y_[i + 1] <- rk$y_[1]

    # Next start: the point just drawn, with bearing = the realised turn angle.
    # dt / crs mirror amt::simulate_path (dt falls back to make_start's default;
    # the output time grid above already uses the original start's dt).
    start <- amt::make_start(
      as.numeric(xy[i + 1, c("x_", "y_")]),
      new.ta,
      time = xy$t_[i],
      crs = attr(kernel$args$start, "crs")
    )
  }

  xy
}

#' Simulate a single movement path (iSSF or GAM)
#'
#' Builds a redistribution kernel from `model` and `env_test` and simulates one
#' path that follows the observed step timing in `stp_data`. The caller is
#' responsible for ensuring `env_test` already carries every covariate layer the
#' model references (HR_edge, HR_center_log, the per-class landcover indicators,
#' etc.).
#'
#' `method` picks the route; the burst / month orchestration below is shared,
#' only the kernel builder and path simulator differ:
#'   * "issf" — amt::redistribution_kernel() + amt::simulate_path(); `model` is a
#'     fitted iSSF (amt make_issf_model / fit_clogit).
#'   * "gam"  — redistribution_kernel_gam() + simulate_path_gam(); `model` is an
#'     mgcv cox.ph GAM, and the tentative movement distributions sl_distr /
#'     ta_distr (e.g. attr(stp.var, "sl_") / attr(stp.var, "ta_")) are needed
#'     when compensate.movement = TRUE (the parametric design).
#'
#' Two code paths inside, chosen automatically based on whether the model has
#' any NDVI term:
#'   * NDVI path — for each burst, walks month-by-month, swapping `env_test$ndvi`
#'     to that month's layer before building a fresh kernel. Required because
#'     NDVI is the only time-varying covariate.
#'   * Simple path — for each burst, builds one kernel and simulates every step
#'     in the burst in a single call. Avoids redundant kernel rebuilds when the
#'     model has no NDVI terms.
#'
#' @param stp_data Step data dataframe (provides timestamps + burst structure)
#' @param env_test Environmental rasters for the simulation extent
#' @param ndvi_test NDVI rasters indexed by month (ignored if the model has no
#'   NDVI terms)
#' @param model Fitted movement model — an iSSF (method = "issf") or an mgcv
#'   cox.ph GAM (method = "gam")
#' @param method Which kernel / simulator to use: "issf" (default) or "gam"
#' @param sl_distr Tentative step-length distribution for the GAM kernel
#'   (attr(stp.var, "sl_")); required for method = "gam" with
#'   compensate.movement = TRUE
#' @param ta_distr Tentative turning-angle distribution for the GAM kernel
#'   (attr(stp.var, "ta_")) or NULL
#' @param compensate.movement GAM only: re-add the tentative kernel + Jacobian
#'   (TRUE for the parametric gamma / von Mises design)
simulate_movement <- function(
  stp_data,
  env_test,
  ndvi_test,
  model,
  method = c("issf", "gam"),
  sl_distr = NULL,
  ta_distr = NULL,
  compensate.movement = TRUE
) {
  method <- match.arg(method)

  # Route-specific pieces. `model_formula` feeds the NDVI check below;
  # `build_kernel` makes a sampled-endpoint kernel at a start point; `run_path`
  # simulates a path from it. Everything after this block (burst / month
  # orchestration) is shared, so the two routes never duplicate that logic.
  if (method == "issf") {
    model_formula <- model$model$formula

    cov_fun <- function(xy, map) {
      xy |>
        amt::extract_covariates(map, where = "both") |>
        amt::time_of_day(include.crepuscule = FALSE, where = "both") |>
        dplyr::mutate(
          tod_start_day = as.integer(tod_start_ == "day"),
          tod_start_night = as.integer(tod_start_ == "night"),
          days = lubridate::yday(t2_) - min(lubridate::yday(t2_)) + 1
        )
    }

    build_kernel <- function(map, start) {
      amt::redistribution_kernel(
        x = model,
        map = map,
        fun = cov_fun,
        start = start,
        landscape = "discrete",
        as.rast = FALSE
      )
    }

    run_path <- function(kernel, n_steps) {
      amt::simulate_path(kernel, n = n_steps)
    }
  } else {
    if (compensate.movement && is.null(sl_distr)) {
      stop(
        "method = 'gam' with compensate.movement = TRUE needs sl_distr (the tentative step-length distribution, e.g. attr(stp.var, 'sl_'))."
      )
    }

    model_formula <- model$formula

    build_kernel <- function(map, start) {
      redistribution_kernel_gam(
        x = model,
        map = map,
        start = start,
        fun = gam_cov_fun,
        sl_distr = sl_distr,
        ta_distr = ta_distr,
        compensate.movement = compensate.movement,
        as.rast = FALSE
      )
    }

    run_path <- function(kernel, n_steps) {
      simulate_path_gam(kernel, n.steps = n_steps)
    }
  }

  # Detect whether this model needs monthly NDVI swaps. If no term mentions
  # "ndvi", the month sub-loop is wasted work and we use the simple
  # one-kernel-per-burst path.
  needs_ndvi <- any(grepl(
    "ndvi",
    all.vars(model_formula),
    ignore.case = TRUE
  ))

  data_step <- stp_data
  bursts <- unique(data_step$burst_)

  # Simulate each burst separately, then combine
  sim_all_bursts <- foreach(b = bursts, .combine = "rbind") %do%
    {
      burst_data <- data_step |> dplyr::filter(burst_ == b)

      if (needs_ndvi) {
        # ---- NDVI path: month-by-month within the burst -----------------
        burst_data$month_group <- lubridate::month(burst_data$t1_)
        months <- unique(burst_data$month_group)

        # we fill sim_burst with simulated paths
        sim_burst <- NULL

        for (mo in months) {
          mo_data <- burst_data |> dplyr::filter(month_group == mo)
          n_steps <- nrow(mo_data)

          # Set NDVI for this month
          env_test$ndvi <- terra::resample(
            ndvi_test[[mo]],
            env_test,
            method = 'near'
          )

          # Start from previous chunk's last point, or burst start
          if (is.null(sim_burst)) {
            start_row <- mo_data[1, c('x1_', 'y1_', 't1_')]
            start_pt <- amt::make_track(
              start_row,
              .x = x1_,
              .y = y1_,
              .t = t1_,
              crs = terra::crs(env_test)
            )
          } else {
            start_row <- sim_burst[nrow(sim_burst), ]
            start_pt <- amt::make_track(
              start_row,
              .x = x_,
              .y = y_,
              .t = t_,
              crs = terra::crs(env_test)
            )
          }

          start_pt <- start_pt |>
            amt::make_start() |>
            amt::mutate(dt = lubridate::hours(4))

          kernel <- build_kernel(env_test, start_pt)

          sim_result <- tryCatch(
            run_path(kernel, n_steps),
            error = function(err) NA
          )

          if (any(is.na(sim_result))) {
            return(NULL)
          }

          sim_burst <- dplyr::bind_rows(
            sim_burst,
            sim_result |> dplyr::select(x_, y_, t_)
          )
        }

        sim_burst |> dplyr::mutate(burst_ = b)
      } else {
        # ---- Simple path: one kernel for the whole burst ----------------
        n_steps <- nrow(burst_data)

        start_row <- burst_data[1, c('x1_', 'y1_', 't1_')]
        start_pt <- amt::make_track(
          start_row,
          .x = x1_,
          .y = y1_,
          .t = t1_,
          crs = terra::crs(env_test)
        ) |>
          amt::make_start() |>
          amt::mutate(dt = lubridate::hours(4))

        kernel <- build_kernel(env_test, start_pt)

        sim_result <- tryCatch(
          run_path(kernel, n_steps),
          error = function(err) NA
        )

        if (any(is.na(sim_result))) {
          return(NULL)
        }

        sim_result |>
          dplyr::select(x_, y_, t_) |>
          dplyr::mutate(burst_ = b)
      }
    }

  sim_all_bursts
}

#' Rename categorical landcover coefficients to match binary raster layer names
#' @param coef_names Character vector of coefficient names from a fitted model
#' @param prefix The categorical variable prefix to look for (default: "wiscland_end")
#' @param levels Non-reference levels of the categorical variable
#' @return Character vector with renamed coefficients
rename_landcover_coefs <- function(
  coef_names,
  prefix = "wiscland_end",
  levels = setdiff(LANDCOVER_LEVELS, "forest")
) {
  # Sort levels longest first to avoid partial matches
  levels <- levels[order(nchar(levels), decreasing = TRUE)]

  for (level in levels) {
    coef_names <- gsub(
      paste0(prefix, level),
      paste0(level, "_end"),
      coef_names,
      fixed = TRUE
    )
  }

  coef_names
}

#' Calculate mean Energy Score between observed and simulated paths
#' @param obs Observed path dataframe with x1_, y1_, t1_ columns
#' @param sim Simulated paths dataframe with x_, y_, t_, nsim columns
#' @return Mean energy score across all matched time steps
calc_energy_score <- function(obs, sim) {
  # Round times to nearest hour for matching
  obs_times <- lubridate::round_date(obs$t1_, unit = "hour")
  sim_times <- lubridate::round_date(sim$t_, unit = "hour")

  # Find shared time steps
  shared_times <- intersect(obs_times, sim_times)

  if (length(shared_times) == 0) {
    warning("No matching time steps between observed and simulated paths")
    return(NA_real_)
  }

  # Filter to shared times
  obs_matched <- obs[obs_times %in% shared_times, ]
  sim_matched <- sim[sim_times %in% shared_times, ]

  es_per_step <- purrr::map_dbl(1:nrow(obs_matched), function(t) {
    obs_xy <- c(obs_matched$x1_[t], obs_matched$y1_[t])
    obs_time <- lubridate::round_date(obs_matched$t1_[t], unit = "hour")

    # Get all sim locations at this time step
    sim_at_t <- sim_matched[
      lubridate::round_date(sim_matched$t_, unit = "hour") == obs_time,
    ]

    if (nrow(sim_at_t) == 0) {
      return(NA_real_)
    }

    # es_sample expects: obs = d-vector, dat = d x n matrix
    sim_matrix <- t(cbind(sim_at_t$x_, sim_at_t$y_))
    scoringRules::es_sample(obs_xy, dat = sim_matrix)
  })

  mean(es_per_step, na.rm = TRUE)
}

#' One-step-ahead log score of observed path under a given model
#'
#' For each observed step, builds the redistribution kernel anchored at the
#' observed start, evaluates it as a raster, and reads off the kernel value at
#' the observed endpoint to give a per-step log probability. The caller is
#' responsible for ensuring `env_test` already carries every covariate layer
#' the model formula references (HR_edge, HR_center_log, etc.).
#'
#' @param stp_data Observed step data for one deer (x1_, y1_, t1_, x2_, y2_, t2_, burst_)
#' @param env_test Cropped environmental rasters (with all model covariates as layers)
#' @param ndvi_test Cropped NDVI rasters (indexed by month)
#' @param issf_train Precomputed iSSF model (from amt::make_issf_model)
#' @return data.frame with burst_, step_index, t1_, logp
onestep_logscore <- function(
  stp_data,
  env_test,
  ndvi_test,
  issf_train
) {
  bursts <- unique(stp_data$burst_)

  results <- foreach(b = bursts, .combine = "rbind") %do%
    {
      burst_data <- stp_data |> dplyr::filter(burst_ == b)
      current_month <- NULL
      step_results <- vector("list", nrow(burst_data))

      for (i in seq_len(nrow(burst_data))) {
        # Update NDVI only when month changes
        mo <- lubridate::month(burst_data$t1_[i])
        if (is.null(current_month) || mo != current_month) {
          env_test$ndvi <- terra::resample(
            ndvi_test[[mo]],
            env_test,
            method = "near"
          )
          current_month <- mo
        }

        # Build start point from observed location
        start_pt <- burst_data[i, c("x1_", "y1_", "t1_")] |>
          amt::make_track(
            .x = x1_,
            .y = y1_,
            .t = t1_,
            crs = terra::crs(env_test)
          ) |>
          amt::make_start() |>
          amt::mutate(dt = lubridate::hours(4))

        # Build redistribution kernel as raster
        kernel_rast <- tryCatch(
          amt::redistribution_kernel(
            x = issf_train,
            map = env_test,
            fun = function(xy, map) {
              xy |>
                amt::extract_covariates(map, where = "both") |>
                amt::time_of_day(
                  include.crepuscule = FALSE,
                  where = "both"
                ) |>
                dplyr::mutate(
                  tod_start_day = as.integer(tod_start_ == "day"),
                  tod_start_night = as.integer(tod_start_ == "night"),
                  days = lubridate::yday(t2_) - min(lubridate::yday(t2_)) + 1
                )
            },
            start = start_pt,
            landscape = "discrete",
            as.rast = TRUE
          ),
          error = function(e) NULL
        )

        if (is.null(kernel_rast)) {
          step_results[[i]] <- data.frame(
            burst_ = b,
            step_index = i,
            t1_ = burst_data$t1_[i],
            logp = NA_real_
          )
          next
        }

        # Extract probability at observed endpoint
        obs_pt <- cbind(burst_data$x2_[i], burst_data$y2_[i])
        p <- terra::extract(kernel_rast$redistribution.kernel, obs_pt)[1, 1]

        lp <- if (!is.na(p) && p > 0) log(p) else NA_real_

        step_results[[i]] <- data.frame(
          burst_ = b,
          step_index = i,
          t1_ = burst_data$t1_[i],
          logp = lp
        )
      }

      dplyr::bind_rows(step_results)
    }

  results
}

#' One-step-ahead log score of observed path under a fitted GAM
#'
#' GAM analogue of onestep_logscore(). For each observed step it builds the GAM
#' redistribution kernel (redistribution_kernel_gam) anchored at the observed
#' start, evaluates it as a normalised raster, and reads off the kernel value at
#' the observed endpoint to give a per-step log probability. Structure, NDVI
#' month-swapping, and the NA / out-of-disc handling all mirror the iSSF version
#' exactly; only the kernel builder differs (redistribution_kernel_gam +
#' gam_cov_fun, with the tentative gamma / von Mises distributions for the
#' parametric movement compensation). The caller is responsible for ensuring
#' `env_test` already carries every covariate layer the model references
#' (HR_edge, HR_center_log, the landcover band, etc.).
#'
#' @param stp_data Observed step data for one deer (x1_, y1_, t1_, x2_, y2_, t2_, burst_)
#' @param env_test Cropped environmental rasters (with all model covariates as layers)
#' @param ndvi_test Cropped NDVI rasters (indexed by month)
#' @param gam_train Fitted mgcv cox.ph GAM (results_gam[[m]]$gam)
#' @param sl_distr Tentative step-length distribution (attr(stp.var, "sl_"));
#'   required when compensate.movement = TRUE
#' @param ta_distr Tentative turning-angle distribution (attr(stp.var, "ta_")) or NULL
#' @param compensate.movement Re-add the tentative kernel + Jacobian (TRUE for the
#'   parametric gamma / von Mises design; FALSE for a uniform-disc / nonp design)
#' @return data.frame with burst_, step_index, t1_, logp
onestep_logscore_gam <- function(
  stp_data,
  env_test,
  ndvi_test,
  gam_train,
  sl_distr = NULL,
  ta_distr = NULL,
  compensate.movement = TRUE
) {
  if (compensate.movement && is.null(sl_distr)) {
    stop(
      "compensate.movement = TRUE needs sl_distr (the tentative step-length distribution, e.g. attr(stp.var, 'sl_'))."
    )
  }

  bursts <- unique(stp_data$burst_)

  results <- foreach(b = bursts, .combine = "rbind") %do%
    {
      burst_data <- stp_data |> dplyr::filter(burst_ == b)
      current_month <- NULL
      step_results <- vector("list", nrow(burst_data))

      for (i in seq_len(nrow(burst_data))) {
        # Update NDVI only when month changes
        mo <- lubridate::month(burst_data$t1_[i])
        if (is.null(current_month) || mo != current_month) {
          env_test$ndvi <- terra::resample(
            ndvi_test[[mo]],
            env_test,
            method = "near"
          )
          current_month <- mo
        }

        # Build start point from observed location
        start_pt <- burst_data[i, c("x1_", "y1_", "t1_")] |>
          amt::make_track(
            .x = x1_,
            .y = y1_,
            .t = t1_,
            crs = terra::crs(env_test)
          ) |>
          amt::make_start() |>
          amt::mutate(dt = lubridate::hours(4))

        # Build GAM redistribution kernel as raster
        kernel_rast <- tryCatch(
          redistribution_kernel_gam(
            x = gam_train,
            map = env_test,
            start = start_pt,
            fun = gam_cov_fun,
            sl_distr = sl_distr,
            ta_distr = ta_distr,
            compensate.movement = compensate.movement,
            as.rast = TRUE
          ),
          error = function(e) NULL
        )

        if (is.null(kernel_rast)) {
          step_results[[i]] <- data.frame(
            burst_ = b,
            step_index = i,
            t1_ = burst_data$t1_[i],
            logp = NA_real_
          )
          next
        }

        # Extract probability at observed endpoint
        obs_pt <- cbind(burst_data$x2_[i], burst_data$y2_[i])
        p <- terra::extract(kernel_rast$redistribution.kernel, obs_pt)[1, 1]

        lp <- if (!is.na(p) && p > 0) log(p) else NA_real_

        step_results[[i]] <- data.frame(
          burst_ = b,
          step_index = i,
          t1_ = burst_data$t1_[i],
          logp = lp
        )
      }

      dplyr::bind_rows(step_results)
    }

  results
}

#' Score agreement between two ctmm models via integrated absolute difference
#' of their model-implied semi-variance functions.
#'
#' Score = 1 - ∫|γ_a - γ_b| / max(∫γ_a, ∫γ_b), integrated in log-Δt so all
#' timescales weigh equally. The denominator uses the larger of the two
#' curve areas, which makes the score symmetric in (a, b) and bounded in
#' [0, 1] for typical SVF shapes (both monotonically rising, non-negative).
#'
#' 1 = curves identical at every lag.
#' 0 = curves disjoint (one is zero where the other is nonzero).
#'
#' Time-lag grid is taken from variogram(ref) so the comparison spans the
#' timescales actually probed by the data.
#'
#' @param ms_data   Fitted ctmm model (typically on the observed track)
#' @param ms_model  Fitted ctmm model (typically on the candidate)
#' @param ref       Telemetry object used to determine the time-lag grid
#'                  (typically the observed track)
#' @return Numeric scalar in [0, 1]; NA if both curves have zero area or the
#'   time-lag grid is degenerate.
svf_score <- function(ms_data, ms_model, ref) {
  # Time lags to evaluate at — bin midpoints from the empirical variogram
  # of the reference telemetry (drop zero lags).
  v <- ctmm::variogram(ref)
  dt_grid <- v$lag[v$lag > 0]

  if (length(dt_grid) < 2) {
    return(NA_real_)
  }

  # Model-implied SVF functions (ctmm internal accessor)
  svf_data <- ctmm:::svf.func(ms_data, moment = TRUE)$svf
  svf_model <- ctmm:::svf.func(ms_model, moment = TRUE)$svf

  curve_data <- vapply(dt_grid, svf_data, numeric(1))
  curve_model <- vapply(dt_grid, svf_model, numeric(1))

  # Integrate in log-Δt so short and long timescales weigh equally
  log_dt <- log(dt_grid)
  diff_curve <- abs(curve_model - curve_data)

  # Trapezoidal integration on (possibly uneven) x-grid
  trap <- function(x, y) sum(diff(x) * (y[-1] + y[-length(y)]) / 2)

  iad <- trap(log_dt, diff_curve)
  data_area <- trap(log_dt, curve_data)
  model_area <- trap(log_dt, curve_model)

  # Normalize by the larger of the two curve areas. Makes the score
  # symmetric and avoids the over/undershoot asymmetry of the data-only
  # normalization.
  norm <- max(data_area, model_area)

  if (norm == 0 || !is.finite(norm)) {
    return(NA_real_)
  }

  1 - iad / norm
}

#' Estimate overlap of utilization distributions
#' @param data Observed paths
#' @param sim Simulated paths
#' @param n_sim number of simulated paths
#' Estimate UD overlap and SVF agreement between observed and simulated paths.
#' @param data Observed paths
#' @param sim Simulated paths
#' @param n_sim number of simulated paths
overlap_ud <- function(data, sim, n_sim) {
  z1_starts <- data |>
    dplyr::select(x = x1_, y = y1_, timestamp = t1_)

  z1_ends <- data |>
    dplyr::group_by(burst_) |>
    dplyr::slice_tail(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(x = x2_, y = y2_, timestamp = t2_)

  z1 <- dplyr::bind_rows(z1_starts, z1_ends) |>
    dplyr::arrange(timestamp) |>
    dplyr::distinct(timestamp, .keep_all = TRUE) |>
    sf::st_as_sf(coords = c("x", "y"), crs = 6610) |>
    sf::st_transform(4326) |>
    dplyr::mutate(
      longitude = sf::st_coordinates(geometry)[, 1],
      latitude = sf::st_coordinates(geometry)[, 2],
      timestamp = timestamp
    ) |>
    sf::st_drop_geometry() |>
    as.data.frame() |>
    ctmm::as.telemetry()

  z2 <- purrr::map(1:n_sim, function(k) {
    tel <- sim |>
      dplyr::filter(nsim == k) |>
      dplyr::select("x_", "y_", "t_") |>
      sf::st_as_sf(coords = c("x_", "y_"), crs = 6610) |>
      sf::st_transform(4326) |>
      dplyr::mutate(
        longitude = sf::st_coordinates(geometry)[, 1],
        latitude = sf::st_coordinates(geometry)[, 2],
        timestamp = t_
      ) |>
      sf::st_drop_geometry() |>
      as.data.frame() |>
      ctmm::as.telemetry()

    ctmm::projection(tel) <- ctmm::projection(z1)
    tel
  })

  # Fit ctmm to observed track (still needed for AKDE bandwidth + SVF curve)
  invisible(capture.output(
    ms1 <- ctmm::ctmm.select(
      z1,
      ctmm::ctmm.guess(z1, interactive = FALSE),
      verbose = FALSE,
      cores = 1
    )
  ))

  # Fit each simulated track using ms1's structure as a guess
  ms2 <- purrr::map(z2, function(z) {
    tryCatch(
      {
        invisible(capture.output(fit <- ctmm::ctmm.fit(z, ms1)))
        fit
      },
      error = function(e) NULL
    )
  })

  keep <- !purrr::map_lgl(ms2, is.null)
  ms2 <- ms2[keep]
  z2 <- z2[keep]
  invisible(capture.output(ms2_avg <- mean(ms2)))

  # ---- UD overlap (unchanged) ----------------------------------------------
  z1_uds <- ctmm::akde(z1, ms1)
  invisible(capture.output(
    z2_uds <- ctmm::akde(z2, ms2, grid = list(r = z1_uds$r, dr = z1_uds$dr)) |>
      mean()
  ))
  bat_uds <- ctmm::overlap(list(z1_uds, z2_uds))$CI[1, 2, 2]

  # ---- SVF score (replaces bat_ctmm) ---------------------------------------
  svf <- svf_score(ms1, ms2_avg, ref = z1)

  list(
    bat_uds = bat_uds,
    svf_score = svf
  )
}

#' @description
#' Plot simulated paths vs. actual (test) path for a single deer.
#' Usage: source this file, then call plot_deer_paths(row_no).

plot_deer_paths <- function(
  row_no,
  data_path = "data_deer_1_119.rds",
  sim_dir = "results",
  model = NULL, # NULL = facet across all non-NA models; or integer index
  sim_alpha = 0.25
) {
  # Load deer data and pick the row
  deer_mvt <- readRDS(data_path)
  if (row_no < 1 || row_no > nrow(deer_mvt)) {
    stop(sprintf("row_no %d out of range (1:%d)", row_no, nrow(deer_mvt)))
  }
  deer_row <- deer_mvt[row_no, ]

  # Actual (test) path
  stp_test <- deer_row$stp_test[[1]]
  obs <- stp_test %>%
    dplyr::select(x_ = x1_, y_ = y1_, t_ = t1_, burst_) %>%
    dplyr::mutate(type = "observed")

  # Simulated paths
  sim_file <- file.path(sim_dir, sprintf("results_sim_%d.rds", row_no))
  if (!file.exists(sim_file)) {
    stop(sprintf("Simulation file not found: %s", sim_file))
  }
  results_sim <- readRDS(sim_file)

  # Build a long data frame of simulated paths across models
  sim_df <- purrr::imap_dfr(results_sim, function(sim_m, m) {
    if (is.null(sim_m) || (length(sim_m) == 1 && is.na(sim_m))) {
      return(NULL)
    }
    sim_m %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(model = as.integer(m))
  })

  if (nrow(sim_df) == 0) {
    stop("No non-NA simulations found for this deer.")
  }

  if (!is.null(model)) {
    sim_df <- sim_df %>% dplyr::filter(model == !!model)
    if (nrow(sim_df) == 0) {
      stop(sprintf("No simulations for model %d", model))
    }
  }

  # Group id so each (model, nsim) draws as its own line
  sim_df <- sim_df %>%
    dplyr::mutate(path_id = paste(model, nsim, sep = "_"))

  # Plot
  p <- ggplot() +
    geom_path(
      data = obs,
      aes(x = x_, y = y_, group = burst_),
      colour = "orange",
      linewidth = 0.7
    ) +
    geom_path(
      data = sim_df,
      aes(x = x_, y = y_, group = path_id),
      colour = "black",
      alpha = sim_alpha,
      linewidth = 0.3
    ) +
    coord_equal() +
    labs(
      title = sprintf(
        "Deer row %d (id: %s) - observed vs simulated (test)",
        row_no,
        if ("id" %in% names(deer_row)) as.character(deer_row$id) else ""
      ),
    ) +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.title = element_blank())

  if (is.null(model)) {
    p <- p + facet_wrap(~model)
  }

  p
}


#' Normalize an stp_test-shaped path to (x, y, burst_, path_id, frame).
#' Assigns sequential frame indices so that each burst starts `burst_gap`
#' frames after the previous burst ended, creating a visual pause during
#' which the meteorite tail can fade out before the next burst head appears.
.normalize_stp <- function(path, id, burst_gap = 5) {
  df <- tibble::tibble(
    x = path$x1_,
    y = path$y1_,
    burst_ = path$burst_
  )

  # Preserve original order while grouping rows within their burst
  df$.order <- seq_len(nrow(df))
  burst_ids <- unique(df$burst_)

  out <- purrr::imap_dfr(burst_ids, function(b, i) {
    sub <- df[df$burst_ == b, ]
    sub$burst_idx <- i
    sub
  })
  out <- out[order(out$.order), ]
  out$.order <- NULL

  # Assign frames: continuous within a burst, `burst_gap` pause between bursts
  rle_b <- rle(out$burst_idx)
  starts <- c(
    1,
    cumsum(rle_b$lengths[-length(rle_b$lengths)]) +
      1 +
      seq_len(length(rle_b$lengths) - 1) * burst_gap
  )
  frame_vec <- integer(nrow(out))
  pos <- 1
  for (k in seq_along(rle_b$lengths)) {
    n_k <- rle_b$lengths[k]
    frame_vec[pos:(pos + n_k - 1)] <- starts[k]:(starts[k] + n_k - 1)
    pos <- pos + n_k
  }
  out$frame <- frame_vec

  out$path_id <- id
  out
}


#' Animate one or two deer paths (stp_test-shaped) with a meteorite-tail
#' moving head. Single panel, colored by path_id.
#'
#' @param path   stp_test-like tibble (x1_, y1_, burst_, ...)
#' @param path2  Optional second path of the same shape
#' @param labels Legend labels for path / path2
#' @param colors Colors for path / path2
#' @param file   Output file (e.g. "plots/deer_42.mp4"). ".mp4" -> av,
#'               ".gif" -> gifski. If NULL, returns animation without saving.
#' @param fps    Frames per second (default 20)
#' @param duration Optional total duration (s); overrides fps-per-step
#' @param wake_length Tail length as fraction of total frames (default 0.05)
#' @param point_size Size of moving head
#' @param width,height,res Output dimensions (pixels) and resolution
#' @param burst_gap Empty frames inserted between bursts for tail fade-out
#' @param step_duration Seconds of animation per data step. Determines clip
#'                      length as n_steps * step_duration. Ignored if
#'                      `duration` is set.
animate_deer_path <- function(
  path,
  path2 = NULL,
  path3 = NULL,
  labels = c("observed", "simulated", "simulated2"),
  colors = c("orange", "blue", "green"),
  file = NULL,
  fps = 30,
  duration = NULL,
  step_duration = 0.15,
  wake_length = 0.05,
  point_size = 3,
  width = 800,
  height = 800,
  res = 150,
  burst_gap = 1
) {
  if (!requireNamespace("gganimate", quietly = TRUE)) {
    stop("Please install gganimate.")
  }

  # Workaround: amt::simulate_path returns n+1 rows (includes start point),
  # so simulated bursts are one step longer than observed bursts. Drop the
  # last row of each burst in sim paths if that pattern is detected.
  trim_sim <- function(sim_path, obs_path) {
    obs_counts <- table(obs_path$burst_)
    sim_counts <- table(sim_path$burst_)
    shared <- intersect(names(obs_counts), names(sim_counts))
    if (
      length(shared) > 0 &&
        all(sim_counts[shared] == obs_counts[shared] + 1)
    ) {
      sim_path <- sim_path %>%
        dplyr::group_by(burst_) %>%
        dplyr::slice(-dplyr::n()) %>%
        dplyr::ungroup()
    }
    sim_path
  }

  # Normalize inputs ---------------------------------------------------------
  df1 <- .normalize_stp(path, labels[1], burst_gap = burst_gap)
  dfs <- list(df1)
  path_levels <- labels[1]

  if (!is.null(path2)) {
    path2 <- trim_sim(path2, path)
    dfs[[length(dfs) + 1]] <- .normalize_stp(
      path2,
      labels[2],
      burst_gap = burst_gap
    )
    path_levels <- c(path_levels, labels[2])
  }
  if (!is.null(path3)) {
    path3 <- trim_sim(path3, path)
    dfs[[length(dfs) + 1]] <- .normalize_stp(
      path3,
      labels[3],
      burst_gap = burst_gap
    )
    path_levels <- c(path_levels, labels[3])
  }

  # Truncate all to the shared number of frames so they end together
  shared_n <- min(vapply(dfs, function(d) max(d$frame), numeric(1)))
  dfs <- lapply(dfs, function(d) d[d$frame <= shared_n, ])
  df <- dplyr::bind_rows(dfs)

  df$path_id <- factor(df$path_id, levels = path_levels)
  # Unique group per (path, burst) so lines/reveal break between bursts
  df$group_id <- interaction(df$path_id, df$burst_, drop = TRUE)

  n_steps <- max(df$frame)

  # Build animation ----------------------------------------------------------
  # Attach color directly (avoid scale_color_manual which can conflict with
  # shadow_wake's color interpolation on some gganimate versions).
  color_map <- stats::setNames(colors[seq_along(path_levels)], path_levels)
  df$color <- color_map[as.character(df$path_id)]

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = x, y = y, group = group_id, colour = color)
  ) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_colour_identity() +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none") +
    gganimate::transition_time(frame) +
    gganimate::shadow_wake(wake_length = wake_length)

  # Render -------------------------------------------------------------------
  anim_args <- list(
    plot = p,
    fps = fps,
    width = width,
    height = height,
    res = res,
    device = "png"
  )
  if (is.null(duration)) {
    duration <- n_steps * step_duration
  }
  anim_args$nframes <- max(round(duration * fps), 30)

  if (!is.null(file)) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "mp4") {
      if (!requireNamespace("av", quietly = TRUE)) {
        stop("Install 'av' for mp4 output.")
      }
      anim_args$renderer <- gganimate::av_renderer()
    } else if (ext == "gif") {
      if (!requireNamespace("gifski", quietly = TRUE)) {
        stop("Install 'gifski' for gif output.")
      }
      anim_args$renderer <- gganimate::gifski_renderer()
    } else {
      stop("file must end in .mp4 or .gif")
    }
  }

  anim <- do.call(gganimate::animate, anim_args)

  if (!is.null(file)) {
    dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
    gganimate::anim_save(file, animation = anim)
  }

  invisible(anim)
}
