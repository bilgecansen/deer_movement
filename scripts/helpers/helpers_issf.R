# iSSF-specific helpers --------------------------------------------------------
#
# Fitting, coefficient tidying and one-step-ahead log scoring for the
# conditional-logistic (amt::fit_issf) path. The GAM equivalents live in
# helpers_gam.R; the two are kept deliberately separate rather than merged
# behind a common interface.
#
# Part of the helper library split out of scripts/helper_functions.R, which
# now sources every file in this folder. Scripts keep sourcing that one
# aggregator, so nothing here needs to be sourced directly.

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

