# Track / step-data wrangling --------------------------------------------------
#
# Turning observed relocations into step-selection data: random-point
# generation and the covariate extraction that attaches landcover, NDVI and
# home-range values to each observed and random step. Shared by the iSSF and
# GAM paths, which both fit on the output.
#
# Part of the helper library split out of scripts/helper_functions.R, which
# now sources every file in this folder. Scripts keep sourcing that one
# aggregator, so nothing here needs to be sourced directly.

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

