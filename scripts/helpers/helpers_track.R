# Track / step-data wrangling --------------------------------------------------
#
# Turning observed relocations into step-selection data: random-point
# generation and the covariate extraction that attaches landcover, NDVI,
# LANDFIRE and home-range values to each observed and random step. Shared by
# the amt and GAM paths, which both fit on the output.
#
# Part of the helper library split out of scripts/helper_functions.R, which
# now sources every file in this folder. Scripts keep sourcing that one
# aggregator, so nothing here needs to be sourced directly.

#' Crop buffer for one deer-year
#'
#' CROP_BUFFER_M, widened to the deer-year's longest step plus CROP_MARGIN_M
#' when that is larger, so no random endpoint (none is further from its start
#' than the longest step) can fall outside the window.
#'
#' @param sl Step lengths (m) of the deer-year's observed steps
#' @return Buffer in metres
deer_crop_buffer <- function(sl) {
  max(CROP_BUFFER_M, max(sl, na.rm = TRUE) + CROP_MARGIN_M)
}

#' Raster crop window around a deer's steps
#'
#' `buffer_m` around the start of every step in `steps`. Random-point
#' generation and covariate extraction both crop to it, from the same steps
#' and the same deer_crop_buffer(), so the window that decides which random
#' endpoints are kept is the window their covariates are read from.
#' terra::crop() uses the extent.
#'
#' @param steps Steps with x1_ / y1_ columns (CRS 6610)
#' @param buffer_m Buffer in metres, from deer_crop_buffer()
#' @return sf polygons, one buffer per distinct start point
step_crop_window <- function(steps, buffer_m) {
  starts <- unique(data.frame(x = steps$x1_, y = steps$y1_))
  sf::st_buffer(
    sf::st_as_sf(starts, coords = c("x", "y"), crs = 6610),
    buffer_m
  )
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
#'   * "amt" — step lengths / turning angles drawn from gamma and von Mises
#'     distributions fitted to the data (amt defaults). Movement covariates are
#'     then corrections to these tentative distributions (the standard amt).
#'     Also usable for GAMs with a *parametric* movement kernel.
#'   * "nonp" — the design the xgb path reads; see nonp_random_points().
#' @param stp_col Column holding the observed steps
#' @param output_col Column to write the random steps to
make_random_pt_extraction <- function(
  data,
  n_pts,
  water,
  model = c("amt", "nonp"),
  stp_col = "stp",
  output_col = "random.stp"
) {
  model <- match.arg(model)
  stopifnot(nrow(data) == 1)

  data_step <- data[[stp_col]][[1]]
  buffer_m <- deer_crop_buffer(data_step$sl_)

  res <- if (model == "amt") {
    amt_random_points(data_step, n_pts, water, buffer_m)
  } else {
    nonp_random_points(data_step, n_pts, water, buffer_m)
  }

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

#' Water at the end of every candidate step
#'
#' The water raster is cropped to the window extract_step_variables() reads,
#' built from the same steps: the ones random_steps() kept, which are the ones
#' that get random points. A random endpoint outside the window would have no
#' Water value and be dropped, so every kept endpoint has covariates.
#'
#' @param rs Output of amt::random_steps()
#' @param water Binary raster for water bodies
#' @param buffer_m Crop buffer, from deer_crop_buffer()
#' @return rs with a Water column
water_at_end <- function(rs, water, buffer_m) {
  window <- step_crop_window(dplyr::filter(rs, case_), buffer_m)
  amt::extract_covariates(rs, terra::crop(water, window), where = "end")
}

#' Random points from fitted gamma / von Mises distributions ("amt" design)
#'
#' Ten times `n_pts` candidates per step, the ones on water dropped, then
#' `n_pts` drawn from the rest.
#'
#' @param data_step Observed steps
#' @param n_pts Random points per step
#' @param water Binary raster for water bodies
#' @param buffer_m Crop buffer, from deer_crop_buffer()
#' @return Observed and random steps
amt_random_points <- function(data_step, n_pts, water, buffer_m) {
  pts <- water_at_end(
    amt::random_steps(data_step, n_control = n_pts * 10),
    water,
    buffer_m
  )
  pts |>
    dplyr::filter(case_ == FALSE, Water == 0) |>
    dplyr::slice_sample(n = n_pts, by = step_id_) |>
    dplyr::bind_rows(pts |> dplyr::filter(case_ == TRUE)) |>
    dplyr::ungroup() |>
    dplyr::select(-Water)
}

#' Random points uniform in distance, in two segments ("nonp" design)
#'
#' Per deer-year, take its NONP_INNER_QUANTILE step length q and its longest
#' step R. Each step gets NONP_INNER_SHARE of its points with distance
#' uniform on (0, q) and the rest uniform on (q, R), turning angle uniform.
#' Uniform in distance puts as many points in each 100 m near the start as
#' far out, where a disc uniform in area puts most of them far away, and the
#' split keeps a long tail of steps from thinning the inner points.
#'
#' Each point then stands for a different amount of ground, recorded as
#'
#'   log_w = log(distance) + log(segment width / points in the segment)
#'           + log(land_share)
#'
#' the log of the land the point stands for, up to a constant shared by the
#' step. Within a segment a point's ground grows with its distance (a ring
#' further out is wider); the segment term spreads the segment over the
#' points it was given; land_share is the share of that step's candidates in
#' the segment that fell on land, since points on water are dropped and the
#' segment's points are then spread over its land alone. A model fitted to
#' these points adds log_w to every endpoint's score, a fixed offset, which
#' makes it the model it would be with points spread evenly over the land.
#' The observed step takes the log_w of the segment its distance falls in.
#'
#' @param data_step Observed steps
#' @param n_pts Random points per step
#' @param water Binary raster for water bodies
#' @param buffer_m Crop buffer, from deer_crop_buffer()
#' @return Observed and random steps, with sample_seg (1 inner, 2 outer),
#'   land_share and log_w
nonp_random_points <- function(data_step, n_pts, water, buffer_m) {
  sl <- data_step$sl_
  edges <- c(
    0,
    stats::quantile(sl, NONP_INNER_QUANTILE, names = FALSE, na.rm = TRUE),
    max(sl, na.rm = TRUE)
  )
  n_seg <- round(n_pts * NONP_INNER_SHARE)
  n_seg <- c(n_seg, n_pts - n_seg)
  stopifnot(all(n_seg >= 1), all(diff(edges) > 0))

  segs <- lapply(1:2, function(k) {
    n_cand <- n_seg[k] * 10
    pts <- water_at_end(
      amt::random_steps(
        data_step,
        n_control = n_cand,
        rand_sl = stats::runif(1e5, edges[k], edges[k + 1]),
        rand_ta = stats::runif(1e5, -pi, pi)
      ),
      water,
      buffer_m
    )
    cand <- dplyr::filter(pts, case_ == FALSE)
    share <- tapply(cand$Water %in% 0, cand$step_id_, mean)
    # A segment with no land candidate for a step keeps no points there; the
    # floor only matters if that step's observed endpoint falls in it.
    share <- pmax(share, 0.5 / n_cand)
    list(
      obs = dplyr::filter(pts, case_ == TRUE),
      kept = cand |>
        dplyr::filter(Water == 0) |>
        dplyr::slice_sample(n = n_seg[k], by = step_id_) |>
        dplyr::mutate(sample_seg = k),
      share = share
    )
  })

  obs <- segs[[1]]$obs |>
    dplyr::mutate(sample_seg = ifelse(sl_ <= edges[2], 1L, 2L))
  res <- dplyr::bind_rows(segs[[1]]$kept, segs[[2]]$kept, obs) |>
    dplyr::ungroup() |>
    dplyr::select(-Water)

  id <- as.character(res$step_id_)
  res$land_share <- ifelse(
    res$sample_seg == 1L,
    segs[[1]]$share[id],
    segs[[2]]$share[id]
  )
  stopifnot(all(res$sl_ > 0), !anyNA(res$land_share))
  res$log_w <- log(res$sl_) +
    log(diff(edges)[res$sample_seg] / n_seg[res$sample_seg]) +
    log(res$land_share)
  res
}

#' Extract environmental variables for each step (single-deer)
#'
#' Operates on a single-row tibble. Builds a cropped env stack
#' (env + NDVI + LANDFIRE + the three HR rasters keyed by id/season/year) and
#' extracts covariates at every step start and end of the deer's random-steps
#' tibble in `data[[random_col]]`.
#'
#' Produces design-matrix columns including:
#'   * wiscland and its per-class indicators at start and end
#'   * forest_edge_start/end (m; + inside forest, - outside)
#'   * HR_bin_start/end, HR_edge_start/end (m), HR_center_start/end (km)
#'   * ndvi_start/end (time-matched)
#'   * oak_mast, oak_dist (m) and famd1..famd5 at start and end
#'   * elevation (m), northness and eastness at start and end
#'
#' @param data Single-row dataframe; must contain id, season, year and the
#'   random-steps column named by `random_col`. The HR raster files
#'   `HRbin_<id>_<season>_<year>.tif`, `HRedge_<id>_<season>_<year>.tif`, and
#'   `HRcenter_<id>_<season>_<year>.tif` must exist in `hr_folder`.
#' @param env Landcover stack for the deer's year+season (from load_landcover)
#' @param ndvi 12-layer monthly NDVI stack for the deer's year (from load_ndvi)
#' @param landfire LANDFIRE layers for the deer's year (from load_landfire)
#' @param topo Elevation, northness and eastness (from load_topo)
#' @param hr_folder Folder containing per-deer HR rasters
#' @param stp_col Column holding the observed steps, which set the crop
#'   buffer exactly as they did in make_random_pt_extraction()
extract_step_variables <- function(
  data,
  env,
  ndvi,
  landfire,
  topo,
  hr_folder = "data/HR",
  random_col = "random.stp",
  output_col = "stp.var",
  stp_col = "stp"
) {
  stopifnot(nrow(data) == 1)

  data_step <- data[[random_col]][[1]]

  crop_extent <- step_crop_window(
    data_step,
    deer_crop_buffer(data[[stp_col]][[1]]$sl_)
  )

  env_cropped <- terra::crop(env, crop_extent)
  ndvi_local <- terra::crop(ndvi, crop_extent)
  landfire_local <- terra::crop(landfire, crop_extent)
  topo_local <- terra::crop(topo, crop_extent)

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
      # Days since the deer-year's first step (day 1), counted on calendar
      # dates in the track's time zone, so a winter track counts on past
      # New Year
      days = as.numeric(lubridate::date(t2_) - min(lubridate::date(t2_))) + 1
    ) |>
    amt::extract_covariates(landfire_local, where = "both") |>
    amt::extract_covariates(topo_local, where = "both")

  data[[output_col]] <- list(data_ssf)
  data
}

#' Age in months at the start of each track
#'
#' The source table records age once, at collaring, as ">2yrs" or "20mo".
#' A deer collared at 20 months is aged forward from the first fix of that
#' animal anywhere in the table. Collaring cannot come after that fix, so
#' the estimate can only understate the true age. Deer collared at over two
#' years are adults in every track and get NA.
#'
#' @param tracks The source track table (library/SW_filtered_deer.RData)
#' @return tibble(key, age_at_collar, age_months), one row per track, keyed
#'   <id>_<season>_<year> as in data/tracks/
deer_track_age <- function(tracks) {
  first_fix <- vapply(tracks$stp, function(s) as.numeric(min(s$t1_)),
                      numeric(1))
  out <- tibble::tibble(
    id = tracks$id,
    key = sprintf("%s_%s_%d", tracks$id, tracks$season,
                  as.integer(tracks$year)),
    age_at_collar = tracks$age.at.col1,
    t0 = first_fix
  )
  stopifnot(all(out$age_at_collar %in% c(">2yrs", "20mo")))
  out <- dplyr::mutate(dplyr::group_by(out, id), t_first = min(t0))
  out$age_months <- ifelse(
    out$age_at_collar == "20mo",
    20 + (out$t0 - out$t_first) / (86400 * 30.44),
    NA_real_
  )
  dplyr::ungroup(out)[c("key", "age_at_collar", "age_months")]
}

#' The canonical deer list for analysis
#'
#' Keeps the deer-season-years every analysis uses:
#'   * females only
#'   * deer collared at 20 months only once they are MIN_AGE_MONTHS old at
#'     the start of the track (deer_track_age); deer collared at over two
#'     years are always kept
#'   * TRAIN_YEARS as training and calibration, TEST_YEAR as test
#'   * only deer-seasons with both a training and a test track
#'   * none of EXCLUDED_DEER_SEASONS
#' keep, excursion, unstable_hr_center, dump and weird are not used.
#'
#' @param tracks The source track table (library/SW_filtered_deer.RData),
#'   whole: ages are counted from each deer's first fix in it
#' @return The rows of `tracks` on the list, with a `set` column, "train" or
#'   "test"
canonical_deer <- function(tracks) {
  age <- deer_track_age(tracks)
  year <- as.integer(tracks$year)
  deer_season <- paste(tracks$id, tracks$season, sep = "_")
  on_list <- tracks$sex == "Female" &
    (tracks$age.at.col1 == ">2yrs" | age$age_months >= MIN_AGE_MONTHS) &
    year %in% c(TRAIN_YEARS, TEST_YEAR) &
    !deer_season %in% EXCLUDED_DEER_SEASONS
  out <- tracks[on_list, ]
  out$set <- ifelse(as.integer(out$year) %in% TRAIN_YEARS, "train", "test")
  out_ds <- paste(out$id, out$season, sep = "_")
  both <- intersect(out_ds[out$set == "train"], out_ds[out$set == "test"])
  out[out_ds %in% both, ]
}

