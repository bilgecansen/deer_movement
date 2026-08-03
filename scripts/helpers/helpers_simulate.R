# Simulation orchestration (both paths) ----------------------------------------
#
# simulate_movement(): the seam between the iSSF and GAM paths. Its `method`
# argument picks the redistribution kernel and path simulator -- amt's for
# iSSF, ours for GAM -- while the burst/month orchestration, NDVI swapping and
# start-point handling around it are shared. This is the only helper that
# reaches into both paths, which is why it sits in its own file rather than in
# either one.
#
# Part of the helper library split out of scripts/helper_functions.R, which
# now sources every file in this folder. Scripts keep sourcing that one
# aggregator, so nothing here needs to be sourced directly.

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

