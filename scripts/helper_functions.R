# Helper functions -------------------------------------------------------------

#' Generate random steps with environmental covariates (single-deer)
#'
#' Operates on a single-row tibble. Generates `n_pts` water-free random
#' control steps per observed step and writes them to `data[[output_col]]`
#' as a list-column.
#'
#' @param data Single-row dataframe with movement steps
#' @param n_pts Number of random points per step
#' @param water Binary raster for water bodies
make_random_pt_extraction <- function(
  data,
  n_pts,
  water = water_binary,
  stp_col = "stp",
  output_col = "random.stp"
) {
  stopifnot(nrow(data) == 1)

  # Crop water raster to this deer's buffer
  data_step <- data[[stp_col]][[1]]
  crop_extent <- sf::st_buffer(
    sf::st_as_sf(data_step, coords = c('x1_', 'y1_'), crs = 6610),
    5000
  )
  water_local <- terra::crop(water, crop_extent)

  # Start with buffer for water removal
  n_random <- ceiling(n_pts * 10)

  # Generate random steps
  random_pts <- data_step |>
    amt::random_steps(n_control = n_random) |>
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
#' @param env Environmental rasters
#' @param ndvi_list NDVI rasters
#' @param hr_folder Folder containing per-deer HR rasters
extract_step_variables <- function(
  data,
  env = env_raster,
  ndvi_list = ndvi_rasters,
  hr_folder = "data/HR",
  random_col = "random.stp",
  output_col = "stp.var"
) {
  stopifnot(nrow(data) == 1)

  data_step <- data[[random_col]][[1]]

  crop_extent <- sf::st_buffer(
    sf::st_as_sf(data_step, coords = c('x1_', 'y1_'), crs = 6610),
    5000
  )

  env_cropped <- terra::crop(env, crop_extent)
  ndvi_year <- ndvi_list[[paste0('ndvi_', data$year)]]
  ndvi_local <- terra::crop(ndvi_year, crop_extent)

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

  lc_levels <- c(
    "central.hardwoods",
    "oak",
    "agriculture",
    "grassland",
    "other"
  )

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

#' Simulate a single movement path
#'
#' Builds a redistribution kernel from `issf_train` and `env_test` and
#' simulates one path that follows the observed step timing in `stp_data`.
#' The caller is responsible for ensuring `env_test` already carries every
#' covariate layer the model formula references (HR_edge, HR_center_log, etc.).
#'
#' Two code paths inside, chosen automatically based on whether the model has
#' any NDVI coefficients:
#'   * NDVI path — for each burst, walks month-by-month, swapping
#'     `env_test$ndvi` to that month's layer before building a fresh kernel.
#'     Required because NDVI is the only time-varying covariate.
#'   * Simple path — for each burst, builds one kernel and simulates every step
#'     in the burst in a single call. Avoids redundant kernel rebuilds when
#'     the model has no NDVI terms (model-formula coefficients drive the
#'     detection).
#'
#' @param stp_data Step data dataframe (provides timestamps + burst structure)
#' @param env_test Environmental rasters for the simulation extent
#' @param ndvi_test NDVI rasters indexed by month (ignored if the model has no
#'   NDVI terms)
#' @param issf_train Fitted iSSF model
simulate_movement <- function(
  stp_data,
  env_test,
  ndvi_test,
  issf_train
) {
  # Detect whether this model needs monthly NDVI swaps. If no coefficient name
  # mentions "ndvi", the month sub-loop is wasted work and we use the simple
  # one-kernel-per-burst path.
  needs_ndvi <- any(grepl(
    "ndvi",
    all.vars(issf_train$model$formula),
    ignore.case = TRUE
  ))

  data_step <- stp_data
  bursts <- unique(data_step$burst_)

  # Covariate-extraction callback used by both code paths
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

          kernel <- amt::redistribution_kernel(
            x = issf_train,
            map = env_test,
            fun = cov_fun,
            start = start_pt,
            landscape = "discrete",
            as.rast = FALSE
          )

          sim_result <- tryCatch(
            amt::simulate_path(kernel, n = n_steps),
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

        kernel <- amt::redistribution_kernel(
          x = issf_train,
          map = env_test,
          fun = cov_fun,
          start = start_pt,
          landscape = "discrete",
          as.rast = FALSE
        )

        sim_result <- tryCatch(
          amt::simulate_path(kernel, n = n_steps),
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
  levels = c(
    "central.hardwoods",
    "oak",
    "agriculture",
    "grassland",
    "other"
  )
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

#' Score agreement between two ctmm models via integrated absolute difference
#' of their model-implied semi-variance functions.
#'
#' Returns a value in roughly [0, 1] where 1 = curves identical at every lag
#' and 0 = curves maximally different. Time-lag grid is taken from
#' variogram(ref) so the comparison spans the timescales actually probed by
#' the data; integration is done in log-Δt so all timescales weigh equally.
#'
#' @param ms_a   First fitted ctmm model
#' @param ms_b   Second fitted ctmm model
#' @param ref    Telemetry object used to determine the time-lag grid
#'               (typically the observed track)
#' @return Numeric scalar; NA if the curves are identically zero
svf_score <- function(ms_a, ms_b, ref) {
  # Time lags to evaluate at — bin midpoints from the empirical variogram
  # of the reference telemetry (drop zero lags).
  v <- ctmm::variogram(ref)
  dt_grid <- v$lag[v$lag > 0]

  if (length(dt_grid) < 2) {
    return(NA_real_)
  }

  # Model-implied SVF functions (ctmm internal accessor)
  svf_a <- ctmm:::svf.func(ms_a, moment = TRUE)$svf
  svf_b <- ctmm:::svf.func(ms_b, moment = TRUE)$svf

  curve_a <- vapply(dt_grid, svf_a, numeric(1))
  curve_b <- vapply(dt_grid, svf_b, numeric(1))

  # Integrate in log-Δt so short and long timescales weigh equally
  log_dt <- log(dt_grid)
  diff_curve <- abs(curve_a - curve_b)
  mean_curve <- (curve_a + curve_b) / 2

  # Trapezoidal integration on (possibly uneven) x-grid
  trap <- function(x, y) sum(diff(x) * (y[-1] + y[-length(y)]) / 2)

  iad <- trap(log_dt, diff_curve)
  norm <- trap(log_dt, mean_curve)

  if (norm == 0 || !is.finite(norm)) {
    return(NA_real_)
  }

  1 - iad / (2 * norm)
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
  z1 <- data |>
    dplyr::select("x1_", "y1_", "t1_") |>
    sf::st_as_sf(coords = c("x1_", "y1_"), crs = 6610) |>
    sf::st_transform(4326) |>
    dplyr::mutate(
      longitude = sf::st_coordinates(geometry)[, 1],
      latitude = sf::st_coordinates(geometry)[, 2],
      timestamp = t1_
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
