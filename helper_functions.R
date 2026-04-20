# Helper functions -------------------------------------------------------------

#' Generate random steps with environmental covariates
#' @param data Dataframe with movement steps
#' @param n_pts Number of random points per step
#' @param water binary raster for water bodies
make_random_pt_extraction <- function(
  data,
  n_pts,
  water = water_binary,
  stp_col = "stp",
  output_col = "random.stp"
) {
  # Pre-crop and wrap water raster per deer
  deer_water <- purrr::map(1:nrow(data), function(i) {
    step_data <- data[i, ][[stp_col]][[1]]
    crop_extent <- sf::st_buffer(
      sf::st_as_sf(step_data, coords = c('x1_', 'y1_'), crs = 6610),
      5000
    )
    terra::wrap(terra::crop(water, crop_extent))
  })

  random.stp <- furrr::future_map(
    1:nrow(data),
    function(i) {
      print(i)

      water_local <- terra::unwrap(deer_water[[i]])

      z <- data[i, ]

      data_step <- z[[stp_col]][[1]]

      # Start with buffer for water removal
      n_random <- ceiling(n_pts * 10)

      # Generate random steps
      random_pts <- data_step |>
        amt::random_steps(n_control = n_random) |>
        amt::extract_covariates(water_local, where = "end")

      # Filter and select final random points
      res <- random_pts |>
        # Select points not on water
        dplyr::filter(case_ == FALSE, Water == 0) |>
        # sample n_pts for each step
        dplyr::slice_sample(n = n_pts, by = step_id_) |>
        dplyr::bind_rows(random_pts |> dplyr::filter(case_ == TRUE)) |>
        dplyr::ungroup() |>
        dplyr::select(-Water)

      # Warning if there are steps with less than random n_pts
      # + 1 is the original step
      final_counts <- table(res$step_id_)
      if (any(final_counts < (n_pts + 1))) {
        # Immediate feedback in console
        message(paste(
          "\nNote: Individual",
          z$id,
          "has steps with missing random points."
        ))
        # Formal warning
        warning(
          paste("Data for individual", z$id, "is incomplete."),
          immediate. = TRUE
        )
      }

      res
    },
    .options = furrr_options(
      packages = c("amt", "terra", "sf", "tidyverse"),
      stdout = FALSE,
      seed = T
    )
  )

  data[[output_col]] <- random.stp
  data
}

#' Extract environmental variables for each step
#' @param data Dataframe with random steps
#' @param env Environmental rasters
#' @param ndvi_list NDVI rasters
extract_step_variables <- function(
  data,
  env = env_raster,
  ndvi_list = ndvi_rasters,
  random_col = "random.stp",
  output_col = "stp.var"
) {
  # Pre-crop and wrap rasters per deer, bundled with deer data
  deer_inputs <- purrr::map(1:nrow(data), function(i) {
    data_row <- data[i, ]
    step_data <- data_row[[random_col]][[1]]

    crop_extent <- sf::st_buffer(
      sf::st_as_sf(step_data, coords = c('x1_', 'y1_'), crs = 6610),
      5000
    )

    ndvi_year <- ndvi_list[[paste0('ndvi_', data_row$year)]]

    list(
      env = terra::wrap(terra::crop(env, crop_extent)),
      ndvi = terra::wrap(terra::crop(ndvi_year, crop_extent)),
      data_row = data_row,
      step_data = step_data
    )
  })

  stp.var <- furrr::future_map(
    deer_inputs,
    function(input) {
      env_local <- terra::unwrap(input$env)
      ndvi_local <- terra::unwrap(input$ndvi)
      data_row <- input$data_row
      step_data <- input$step_data

      # Add HR distance to environment
      median_pt <- terra::vect(
        cbind(data_row$x_median, data_row$y_median),
        crs = terra::crs(env_local)
      )

      ## Any layer in env works for distance calculation
      env_local$HR <- NA
      env_local$HR <- terra::distance(env_local$HR, median_pt) / 1000

      lc_levels <- c(
        "central.hardwoods",
        "oak",
        "agriculture",
        "grassland",
        "other"
      )

      # Process step data
      data_ssf <- step_data |>
        amt::extract_covariates(env_local, where = 'both') |>
        amt::extract_covariates_var_time(
          ndvi_local,
          max_time = lubridate::days(31),
          when = "any",
          where = "both",
          name_covar = "ndvi"
        ) |>
        # Setting the intercept to the 'forest' land type
        dplyr::mutate(
          wiscland_start = factor(wiscland_start, levels = lc_levels),
          wiscland_end = factor(wiscland_end, levels = lc_levels)
        ) |>
        amt::time_of_day(include.crepuscule = F, where = "both") |>
        dplyr::mutate(
          tod_start_day = as.integer(tod_start_ == "day"),
          tod_start_night = as.integer(tod_start_ == "night"),
          days = lubridate::yday(t2_) - min(lubridate::yday(t2_)) + 1
        )
      data_ssf
    },
    .options = furrr_options(
      packages = c("amt", "terra", "sf", "tidyverse"),
      stdout = FALSE,
      seed = T
    )
  )

  data[[output_col]] <- stp.var
  data
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
#' @param stp_data Step data dataframe
#' @param x_median Median x coordinate of home range center
#' @param y_median Median y coordinate of home range center
#' @param env_test Environmental rasters for test data
#' @param ndvi_test ndvi rasters for test data
#' @param issf_train Fitted iSSF model
simulate_movement <- function(
  stp_data,
  x_median,
  y_median,
  env_test,
  ndvi_test,
  issf_train
) {
  # Calculate distance to home range center
  median_pt <- terra::vect(
    cbind(x_median, y_median),
    crs = terra::crs(env_test[[1]])
  )

  env_test$HR <- NA
  env_test$HR <- terra::distance(env_test$HR, median_pt) / 1000

  step_data <- stp_data
  bursts <- unique(step_data$burst_)

  # Simulate each burst separately, then combine
  sim_all_bursts <- foreach(b = bursts, .combine = "rbind") %do%
    {
      burst_data <- step_data |> dplyr::filter(burst_ == b)

      # Split burst into monthly groups
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
          fun = function(xy, map) {
            xy |>
              amt::extract_covariates(map, where = "both") |>
              amt::time_of_day(include.crepuscule = FALSE, where = "both") |>
              dplyr::mutate(
                tod_start_day = as.integer(tod_start_ == "day"),
                tod_start_night = as.integer(tod_start_ == "night"),
                days = lubridate::yday(t2_) - min(lubridate::yday(t2_)) + 1
              )
          },
          start = start_pt,
          landscape = "discrete",
          as.rast = F
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

#' Run simulate_movement across all deer and replications for a given model
#' @param deer_inputs_base List of per-deer data (cropped rasters, stp_test, coordinates)
#' @param model_sims List of pre-built simulation models per deer (NULL for skipped deer)
#' @param n_sim Number of simulations per deer
run_simulations <- function(
  deer_inputs_base,
  model_sims,
  n_sim = 10
) {
  # Combine base inputs with model-specific info
  deer_inputs <- purrr::map(1:length(deer_inputs_base), function(i) {
    input <- deer_inputs_base[[i]]
    input$skip <- is.null(model_sims[[i]])
    input$model_sim <- model_sims[[i]]
    input
  })

  furrr::future_map(
    deer_inputs,
    function(input) {
      if (input$skip) {
        return(NA)
      }

      env_local <- terra::unwrap(input$crop_env)
      ndvi_local <- terra::unwrap(input$crop_ndvi)

      foreach(h = 1:n_sim, .combine = "rbind") %do%
        {
          res <- simulate_movement(
            stp_data = input$stp_test,
            x_median = input$x_median,
            y_median = input$y_median,
            env_test = env_local,
            ndvi_test = ndvi_local,
            issf_train = input$model_sim
          )
          res$nsim <- h
          res
        }
    },
    .options = furrr_options(
      packages = c("amt", "terra", "sf", "dplyr", "lubridate", "foreach"),
      stdout = FALSE,
      seed = TRUE
    )
  )
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

#' One-step-ahead Energy Score between observed and simulated paths
#' @param stp_data Observed step data for one deer (x1_, y1_, t1_, x2_, y2_, t2_, burst_)
#' @param x_median Home range center x coordinate
#' @param y_median Home range center y coordinate
#' @param env_test Cropped environmental rasters
#' @param ndvi_test Cropped NDVI rasters (indexed by month)
#' @param issf_train Precomputed iSSF model (from amt::make_issf_model)
#' @param n_sim Number of single-step simulations per observed step
#' @return data.frame with burst_, step_index, t1_, energy_score, spread, n_sim_success
onestep_energy_score <- function(
  stp_data,
  x_median,
  y_median,
  env_test,
  ndvi_test,
  issf_train,
  n_sim = 10
) {
  # Setup: compute HR distance layer once
  median_pt <- terra::vect(
    cbind(x_median, y_median),
    crs = terra::crs(env_test[[1]])
  )
  env_test$HR <- NA
  env_test$HR <- terra::distance(env_test$HR, median_pt) / 1000

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

        # Build redistribution kernel
        kernel <- tryCatch(
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
            as.rast = FALSE
          ),
          error = function(e) NULL
        )

        if (is.null(kernel)) {
          step_results[[i]] <- data.frame(
            burst_ = b,
            step_index = i,
            t1_ = burst_data$t1_[i],
            energy_score = NA_real_,
            spread = NA_real_,
            n_sim_success = 0L
          )
          next
        }

        # Simulate n_sim single steps
        sim_endpoints <- matrix(NA_real_, nrow = 2, ncol = n_sim)
        n_success <- 0L

        for (s in seq_len(n_sim)) {
          sim_result <- tryCatch(
            amt::simulate_path(kernel, n = 1),
            error = function(e) NULL
          )
          if (!is.null(sim_result) && nrow(sim_result) >= 2) {
            n_success <- n_success + 1L
            sim_endpoints[1, n_success] <- sim_result$x_[2]
            sim_endpoints[2, n_success] <- sim_result$y_[2]
          }
        }

        # Compute energy score and spread
        es <- NA_real_
        spread <- NA_real_
        if (n_success >= 2) {
          sim_valid <- sim_endpoints[, seq_len(n_success), drop = FALSE]
          obs_xy <- c(burst_data$x2_[i], burst_data$y2_[i])
          es <- scoringRules::es_sample(obs_xy, dat = sim_valid)
          spread <- mean(dist(t(sim_valid)))
        }

        step_results[[i]] <- data.frame(
          burst_ = b,
          step_index = i,
          t1_ = burst_data$t1_[i],
          energy_score = es,
          spread = spread,
          n_sim_success = n_success
        )
      }

      dplyr::bind_rows(step_results)
    }

  results
}

#' One-step-ahead log-likelihood of observed path under a given model
#' @param stp_data Observed step data for one deer (x1_, y1_, t1_, x2_, y2_, t2_, burst_)
#' @param x_median Home range center x coordinate
#' @param y_median Home range center y coordinate
#' @param env_test Cropped environmental rasters
#' @param ndvi_test Cropped NDVI rasters (indexed by month)
#' @param issf_train Precomputed iSSF model (from amt::make_issf_model)
#' @return data.frame with burst_, step_index, t1_, logp
onestep_loglik <- function(
  stp_data,
  x_median,
  y_median,
  env_test,
  ndvi_test,
  issf_train
) {
  # Setup: compute HR distance layer once
  median_pt <- terra::vect(
    cbind(x_median, y_median),
    crs = terra::crs(env_test[[1]])
  )
  env_test$HR <- NA
  env_test$HR <- terra::distance(env_test$HR, median_pt) / 1000

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

#' Estimate overlap of utilization distributions
#' @param data Observed paths
#' @param sim Simulated paths
#' @param n_sim number of simulated paths
overlap_ud <- function(data, sim, n_sim) {
  z1 <- data |>
    dplyr::select("x1_", "y1_", "t1_") |>
    sf::st_as_sf(coords = c("x1_", "y1_"), crs = 6610) |>
    sf::st_transform(4326) |> # Transform to WGS84 lat/long
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
      sf::st_transform(4326) |> # Transform to WGS84 lat/long
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

  # Fit single guestimated ctmm model (silenced)
  invisible(capture.output(
    ms1 <- ctmm::ctmm.select(
      z1,
      ctmm::ctmm.guess(z1, interactive = FALSE),
      verbose = F,
      cores = 1
    )
  ))

  ms2 <- purrr::map(z2, function(z) {
    tryCatch(
      {
        invisible(capture.output(
          fit <- ctmm::ctmm.fit(z, ms1)
        ))
        fit
      },
      error = function(e) NULL
    )
  })

  # Remove failed fits and their corresponding telemetry data
  keep <- !purrr::map_lgl(ms2, is.null)
  ms2 <- ms2[keep]
  z2 <- z2[keep]
  invisible(capture.output(ms2_avg <- mean(ms2)))

  z1_uds <- ctmm::akde(z1, ms1)
  invisible(capture.output(
    z2_uds <- ctmm::akde(z2, ms2, grid = list(r = z1_uds$r, dr = z1_uds$dr)) |>
      mean()
  ))

  # Estimate the overlap of UDs
  bat_uds <- ctmm::overlap(list(z1_uds, z2_uds))$CI[1, 2, 2]
  bat_ctmm <- ctmm::overlap(list(ms1, ms2_avg))$CI[1, 2, 2]

  list(
    bat_uds = bat_uds,
    bat_ctmm = bat_ctmm
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
