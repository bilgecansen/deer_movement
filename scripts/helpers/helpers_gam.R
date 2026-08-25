# GAM-specific helpers ---------------------------------------------------------
#
# Everything for the mgcv cox.ph (Klappstein et al. 2024) path: reshaping step
# data into stratified Cox-PH form, fitting and diagnosing the smooths, and the
# custom redistribution kernel used for both simulation and log scoring. The
# iSSF equivalents live in helpers_issf.R.
#
# Part of the helper library split out of scripts/helper_functions.R, which
# now sources every file in this folder. Scripts keep sourcing that one
# aggregator, so nothing here needs to be sourced directly.

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
      # Keep ALL declared landcover levels (LANDCOVER_LEVELS), not just the ones
      # this deer visited. A class with no training steps gets its fs / re
      # coefficients pinned at ~0 -- i.e. the population-average response -- so
      # prediction (simulation / scoring) on a landcover class the deer never
      # visited works with a plain predict(): the per-class fs deviation falls
      # back to 0 and the global s(ndvi_end) carries the response (this is the
      # Model GS "predict unobserved levels" property). Verified to leave the
      # fits for visited classes (coefficients and variance components) unchanged.
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
#'   * "removed"     edf <= removed_edf  — shrunk to ~zero. With select = TRUE
#'                                         any smooth can reach this; with
#'                                         select = FALSE only smooths with no
#'                                         unpenalised null space (cc cyclic, re,
#'                                         fs) can, when REML drives their
#'                                         variance to ~0 (no time-of-day
#'                                         modulation / no between-class variance)
#'   * "near-linear" edf <= linear_edf   — collapsed to the penalty null space
#'                                         (a straight line = no nonlinearity);
#'                                         the floor for ordinary tp/cr smooths
#'                                         (HR, global NDVI) when select = FALSE
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

      # Incoming heading: the ABSOLUTE bearing of the preceding step in this
      # burst. The kernel turns a candidate endpoint into a turning angle with
      # ta_ = bearing - start$ta_, so start$ta_ must carry that heading -- it is
      # a reference direction, not a turning angle, despite the name.
      #
      # make_start() on a single point cannot know the heading and returns 0,
      # i.e. due east. Scoring every step from a start built that way evaluates
      # each one as though the deer had just been travelling east, so cos(ta_) is
      # measured off the wrong baseline. Measured on 7193_fa_2020: per-step log p
      # off by up to 0.81 (sd 0.44), and because each model fits its own cos(ta_)
      # coefficient the error does not cancel in delta_logp -- it shifted by
      # ~4.4 for a median-length deer, against a gate-3 threshold of 3.
      #
      # The first step of a burst has no preceding step, so no heading exists for
      # it. Rather than invent one, it is skipped (logp = NA) -- inventing a
      # heading is what produced the bug.
      burst_data$prev_head <- dplyr::lag(
        atan2(
          burst_data$y2_ - burst_data$y1_,
          burst_data$x2_ - burst_data$x1_
        )
      )

      current_month <- NULL
      step_results <- vector("list", nrow(burst_data))

      for (i in seq_len(nrow(burst_data))) {
        # No heading recoverable (first step of the burst): not scoreable.
        if (is.na(burst_data$prev_head[i])) {
          step_results[[i]] <- data.frame(
            burst_ = b,
            step_index = i,
            t1_ = burst_data$t1_[i],
            logp = NA_real_,
            status = "skipped_no_heading"
          )
          next
        }
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

        # Build the start from the observed location AND the observed incoming
        # heading. make_start() is given ta_ explicitly here; called without it
        # (as on a bare one-row track) it silently defaults to 0.
        start_pt <- amt::make_start(
          c(burst_data$x1_[i], burst_data$y1_[i]),
          ta_ = burst_data$prev_head[i],
          time = burst_data$t1_[i],
          dt = lubridate::hours(4),
          crs = terra::crs(env_test)
        )

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
            logp = NA_real_,
            status = "failed_kernel"
          )
          next
        }

        # Extract probability at observed endpoint.
        #
        # An NA here is not a numerical accident: the kernel raster spans only
        # the candidate disc of radius max.dist (the 0.99 quantile of the
        # tentative gamma), so an observed step LONGER than that lands outside it
        # and cannot be scored. Verified across three deer -- the count of
        # unscoreable steps equalled the count of steps exceeding max.dist
        # exactly, and the longest scored step was always just under it. This is
        # a property of the disc, not of the cropped raster: the failing steps
        # sat no closer to the map edge than the scored ones, and the
        # tolerance.outside warning never fired.
        obs_pt <- cbind(burst_data$x2_[i], burst_data$y2_[i])
        p <- terra::extract(kernel_rast$redistribution.kernel, obs_pt)[1, 1]

        st <- if (is.na(p)) {
          "failed_outside_disc"
        } else if (p <= 0) {
          "failed_zero_density"
        } else {
          "ok"
        }
        lp <- if (identical(st, "ok")) log(p) else NA_real_

        step_results[[i]] <- data.frame(
          burst_ = b,
          step_index = i,
          t1_ = burst_data$t1_[i],
          logp = lp,
          status = st
        )
      }

      dplyr::bind_rows(step_results)
    }

  results
}

