# Scoring metrics (both paths) -------------------------------------------------
#
# Model-agnostic scores computed from observed vs. simulated paths: energy
# score, semi-variance function agreement, and utilisation-distribution
# overlap. Consumed by the compute_es_* and run_udoverlap_* runners on both the
# amt and GAM sides. Model-specific scores (log score) live with their path.
#
# Part of the helper library split out of scripts/helper_functions.R, which
# now sources every file in this folder. Scripts keep sourcing that one
# aggregator, so nothing here needs to be sourced directly.

#' Calculate mean Energy Score between observed and simulated paths
#' @param obs Observed path dataframe with x1_, y1_, t1_ columns
#' @param sim Simulated paths dataframe with x_, y_, t_, nsim columns
#' @return Mean energy score across all matched time steps
calc_energy_score <- function(obs, sim) {
  es <- energy_score_steps(obs, sim)
  if (nrow(es) == 0) {
    warning("No matching time steps between observed and simulated paths")
    return(NA_real_)
  }
  mean(es$es, na.rm = TRUE)
}

#' Energy score of the observed position at each matched time step
#'
#' The per-step values that calc_energy_score() averages. An observed step
#' is matched to the simulated positions whose time rounds to the same hour
#' as its start.
#'
#' @param obs Observed path dataframe with x1_, y1_, t1_ columns
#' @param sim Simulated paths dataframe with x_, y_, t_, nsim columns
#' @return tibble(t, es), one row per matched observed step; no rows when no
#'   time matches
energy_score_steps <- function(obs, sim) {
  # Round times to nearest hour for matching
  obs_times <- lubridate::round_date(obs$t1_, unit = "hour")
  sim_times <- lubridate::round_date(sim$t_, unit = "hour")

  # Find shared time steps
  shared_times <- intersect(obs_times, sim_times)

  if (length(shared_times) == 0) {
    return(tibble::tibble(t = obs$t1_[0], es = numeric(0)))
  }

  # Filter to shared times
  obs_matched <- obs[obs_times %in% shared_times, ]
  sim_matched <- sim[sim_times %in% shared_times, ]
  sim_hours <- sim_times[sim_times %in% shared_times]

  es_per_step <- purrr::map_dbl(1:nrow(obs_matched), function(t) {
    obs_xy <- c(obs_matched$x1_[t], obs_matched$y1_[t])
    obs_time <- lubridate::round_date(obs_matched$t1_[t], unit = "hour")

    # Get all sim locations at this time step
    sim_at_t <- sim_matched[sim_hours == obs_time, ]

    if (nrow(sim_at_t) == 0) {
      return(NA_real_)
    }

    # es_sample expects: obs = d-vector, dat = d x n matrix
    sim_matrix <- t(cbind(sim_at_t$x_, sim_at_t$y_))
    scoringRules::es_sample(obs_xy, dat = sim_matrix)
  })

  tibble::tibble(t = obs_matched$t1_, es = es_per_step)
}

#' Share of observed steps longer than the energy score at their start
#'
#' Each observed step is paired with the energy score of its start
#' position against the simulated positions at that time
#' (energy_score_steps()). The share of steps whose length exceeds that
#' score is the share of steps in which the simulated positions are, on
#' this score, closer to the real deer than the distance it then moved.
#' Each step faces its own score, so a few fixes where the simulations are
#' far off count only against their own steps.
#'
#' @param obs Observed steps with x1_, y1_, t1_ and sl_
#' @param sim Simulated paths dataframe with x_, y_, t_, nsim columns
#' @return Share in [0, 1]; NA when no step has a score
calc_p_step <- function(obs, sim) {
  es <- energy_score_steps(obs, sim)
  sl <- obs$sl_[match(es$t, obs$t1_)]
  ok <- !is.na(sl) & !is.na(es$es)
  if (!any(ok)) {
    return(NA_real_)
  }
  mean(sl[ok] > es$es[ok])
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
  # symmetric and avoids the over/undershoot asymmetry of the
  # data-only normalization.
  norm <- max(data_area, model_area)

  if (norm == 0 || !is.finite(norm)) {
    return(NA_real_)
  }

  1 - iad / norm
}

#' ctmm telemetry from projected coordinates
#'
#' Built the way overlap_ud() builds its tracks: EPSG 6610 coordinates go to
#' longitude/latitude and back into ctmm's own projection.
#'
#' @param x,y Coordinates in EPSG 6610 (metres)
#' @param t Times (POSIXct)
#' @param projection Optional ctmm projection to put the track in, so that
#'   several tracks share one coordinate frame
#' @return A ctmm telemetry object
telemetry_6610 <- function(x, y, t, projection = NULL) {
  tel <- tibble::tibble(x = x, y = y, timestamp = t) |>
    sf::st_as_sf(coords = c("x", "y"), crs = 6610) |>
    sf::st_transform(4326) |>
    dplyr::mutate(
      longitude = sf::st_coordinates(geometry)[, 1],
      latitude = sf::st_coordinates(geometry)[, 2]
    ) |>
    sf::st_drop_geometry() |>
    as.data.frame() |>
    ctmm::as.telemetry()
  if (!is.null(projection)) {
    ctmm::projection(tel) <- projection
  }
  tel
}

#' Agreement between the observed and the simulated paths' variograms
#'
#' Compares the variograms themselves rather than ctmm models fitted to the
#' paths. The observed path's empirical variogram is set against the mean of
#' the simulated paths' empirical variograms (ctmm's mean() for
#' variograms, which averages them lag by lag, weighted by DOF). Both use
#' ctmm's exact algorithm (fast = FALSE): the default fast one can sit
#' well below the plain average over pairs of fixes.
#'
#' The simulated paths must already be cut to the observed fix times, so
#' every variogram is built from the same times and shares one lag grid.
#' ctmm leaves out lag bins that are nearly empty, and which ones can differ
#' between tracks whose times differ by a minute, so bins are matched by
#' lag. Bins with no pairs of fixes (DOF < 1) in either curve are dropped.
#'
#' Score = 1 - integral |g_obs - g_sim| / max(integral g_obs,
#' integral g_sim), integrated over log lag, as in svf_score(), and only up
#' to half the observed path's span: a single path's variogram rests on few
#' independent pairs beyond that.
#'
#' @param obs Observed steps with x1_, y1_, t1_, x2_, y2_, t2_, burst_
#' @param sim Simulated paths (nsim, x_, y_, t_), cut to the observed fix
#'   times
#' @param dt Lag bin width in seconds
#' @return Numeric scalar in [0, 1]; 1 = identical curves
variogram_agreement <- function(obs, sim, dt = 4 * 3600) {
  last <- !duplicated(obs$burst_, fromLast = TRUE)
  fixes <- dplyr::bind_rows(
    tibble::tibble(x = obs$x1_, y = obs$y1_, t = obs$t1_),
    tibble::tibble(x = obs$x2_[last], y = obs$y2_[last], t = obs$t2_[last])
  ) |>
    dplyr::arrange(t) |>
    dplyr::distinct(t, .keep_all = TRUE)
  z1 <- telemetry_6610(fixes$x, fixes$y, fixes$t)
  v_obs <- ctmm::variogram(z1, dt = dt, fast = FALSE, trace = FALSE)
  v_sims <- lapply(sort(unique(sim$nsim)), function(k) {
    p <- sim[sim$nsim == k, ]
    z <- telemetry_6610(p$x_, p$y_, p$t_, ctmm::projection(z1))
    ctmm::variogram(z, dt = dt, fast = FALSE, trace = FALSE)
  })
  v_sim <- mean(v_sims)

  obs_bins <- tibble::tibble(bin = round(v_obs$lag / dt), lag = v_obs$lag,
                             g_obs = v_obs$SVF, dof_obs = v_obs$DOF)
  sim_bins <- tibble::tibble(bin = round(v_sim$lag / dt), g_sim = v_sim$SVF,
                             dof_sim = v_sim$DOF)
  half_span <- diff(range(z1$t)) / 2
  v <- dplyr::inner_join(obs_bins, sim_bins, by = "bin") |>
    dplyr::filter(lag > 0, lag <= half_span, dof_obs >= 1, dof_sim >= 1)
  if (nrow(v) < 2) {
    return(NA_real_)
  }

  # Trapezoidal integration on the (uneven) log-lag grid
  trap <- function(x, y) sum(diff(x) * (y[-1] + y[-length(y)]) / 2)
  x <- log(v$lag)
  norm <- max(trap(x, v$g_obs), trap(x, v$g_sim))
  if (norm == 0 || !is.finite(norm)) {
    return(NA_real_)
  }
  1 - trap(x, abs(v$g_obs - v$g_sim)) / norm
}

#' Estimate UD overlap and SVF agreement between observed and simulated paths.
#' @param data Observed paths
#' @param sim Simulated paths
#' @param n_sim number of simulated paths
#' @return list(bat_uds, svf_agree): Bhattacharyya UD overlap and SVF agreement,
#'   both in [0, 1] with 1 = perfect agreement.
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

  # ---- SVF agreement (replaces bat_ctmm) ------------------------------------
  # Named svf_agree, not svf_score, so the value is never confused with the
  # svf_score() function that produces it. Files written before this rename
  # carry the element as `svf_score`; the filter scripts read either name.
  svf <- svf_score(ms1, ms2_avg, ref = z1)

  list(
    bat_uds = bat_uds,
    svf_agree = svf
  )
}


#' Wasserstein-1 distance between two samples of one variable
#'
#' The area between the two empirical cumulative distribution functions,
#' the same as the mean absolute difference between their quantile
#' functions: how far, in the variable's own units, the values of one sample
#' have to move on average to become the other. The samples may differ in
#' size.
#'
#' @param a,b Numeric samples; NA values are dropped
#' @return The distance, in the units of the samples
wasserstein1 <- function(a, b) {
  a <- sort(a[!is.na(a)])
  b <- sort(b[!is.na(b)])
  stopifnot(length(a) > 0, length(b) > 0)
  x <- sort(unique(c(a, b)))
  if (length(x) < 2) {
    return(0)
  }
  # Both CDFs are flat between consecutive values of the pooled sample
  fa <- findInterval(x[-length(x)], a) / length(a)
  fb <- findInterval(x[-length(x)], b) / length(b)
  sum(abs(fa - fb) * diff(x))
}
