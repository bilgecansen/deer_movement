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
  random.stp <- foreach(i = 1:nrow(data)) %do%
    {
      print(i)

      z <- data[i, ]

      data_step <- z[[stp_col]][[1]]

      # Start with buffer for water removal
      n_random <- ceiling(n_pts * 10)

      # Generate random steps
      random_pts <- data_step %>%
        amt::random_steps(n_control = n_random) %>%
        amt::extract_covariates(water, where = "end")

      # Filter and select final random points
      res <- random_pts %>%
        # Select points not on water
        filter(case_ == FALSE, Water == 0) %>%
        # sample n_pts for each step
        slice_sample(n = n_pts, by = step_id_) %>%
        bind_rows(random_pts %>% filter(case_ == TRUE)) %>%
        ungroup() %>%
        select(-Water)

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
          paste("Data for indiviudal", z$id, "is incomplete."),
          immediate. = TRUE
        )
      }

      res
    }

  data[[output_col]] <- random.stp
  data
}

#' Extract environmental variables for each step
#' @param data Dataframe with random steps
#' @param env_list Environmental rasters
#' @param ndvi_list NDVI rasters
extract_step_variables <- function(
  data,
  env = env_raster$wiscland,
  ndvi_list = ndvi_rasters,
  random_col = "random.stp",
  output_col = "stp.var"
) {
  stp.var <- foreach(i = 1:nrow(data)) %do%
    {
      data_row <- data[i, ]

      # Add HR distance to environment
      median_pt <- vect(
        cbind(data_row$x_median, data_row$y_median),
        crs = crs(env)
      )

      ## Any layer in env works for distance calculation
      env$HR <- NA
      env$HR <- values(terra::distance(env$wiscland, median_pt)) / 1000

      ndvi <- ndvi_list[[paste0('ndvi_', data_row$year)]]

      lc_levels <- c(
        "forest",
        "agriculture",
        "grassland",
        "hay",
        "oak",
        "central.hardwoods",
        "other.forest",
        "other"
      )

      # Process step data
      data_ssf <- data_ssf <- data_row[[random_col]][[1]] %>%
        amt::extract_covariates(env, where = 'both') %>%
        amt::extract_covariates_var_time(
          ndvi,
          max_time = days(31),
          when = "any",
          where = "both",
          name_covar = "ndvi"
        ) %>%
        # Setting the intercept to the 'forest' land type
        mutate(
          wiscland_start = factor(wiscland_start, levels = lc_levels),
          wiscland_end = factor(wiscland_end, levels = lc_levels)
        ) %>%
        amt::time_of_day(include.crepuscule = F, where = "both") %>%
        mutate(
          tod_start_day = as.integer(tod_start_ == "day"),
          tod_start_night = as.integer(tod_start_ == "night"),
          days = lubridate::yday(t2_) - min(lubridate::yday(t2_)) + 1
        )
      data_ssf
    }

  data[[output_col]] <- stp.var
  data
}

#' Fit iSSF models to all individuals
#' @param data Dataframe with step variables in stp.var column
#' @param formula Model formula string
#' @return List with iss (fitted models), coeff (tidy coefficients), aic (AIC values)
fit_mods <- function(data, formula) {
  iss <- map(data$stp.var.train, function(ssf_data) {
    tryCatch(
      ssf_data %>% amt::fit_issf(as.formula(formula), model = TRUE),
      error = function(err) "Error"
    )
  })

  coeff <- map(iss, function(mod) {
    if (is.character(mod)) {
      return(NA)
    }
    broom::tidy(mod$model)
  })

  aic <- map_dbl(iss, function(mod) {
    if (is.character(mod)) {
      return(NA_real_)
    }
    AIC(mod$model)
  })

  list(iss = iss, coeff = coeff, aic = aic)
}

#' Simulate a single movement path
#' @param deer_data Deer movement data (e.g., deer_mvt)
#' @param env_test Environmental rasters for test data
#' @param ndvi_test ndvi rasters for test data
#' @param issf_train Fitted iSSF model
simulate_movement <- function(
  deer_data,
  env_test,
  ndvi_test,
  issf_train,
  stp_col = "stp_test"
) {
  # Calculate distance to home range center
  median_pt <- terra::vect(
    cbind(deer_data$x_median, deer_data$y_median),
    crs = crs(env_test[[1]])
  )

  env_test$HR <- NA
  env_test$HR <- terra::distance(env_test$HR, median_pt) / 1000

  step_data <- deer_data[[stp_col]][[1]]
  bursts <- unique(step_data$burst_)

  # Simulate each burst separately, then combine
  sim_all_bursts <- foreach(b = bursts, .combine = "rbind") %do%
    {
      burst_data <- step_data %>% filter(burst_ == b)

      # Split burst into monthly groups
      burst_data$month_group <- format(burst_data$t1_, "%Y-%m")
      months <- unique(burst_data$month_group)

      sim_burst <- NULL

      for (mo in months) {
        mo_data <- burst_data %>% filter(month_group == mo)
        n_steps <- nrow(mo_data)

        # Set NDVI for this month
        mo_mid <- median(mo_data$t1_)
        env_test$ndvi <- resample(
          ndvi_test[[which.min(abs(
            as.Date(mo_mid, tz = 'America/Chicago') - terra::time(ndvi_test)
          ))]],
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
            crs = crs(env_test)
          )
        } else {
          start_row <- sim_burst[nrow(sim_burst), ]
          start_pt <- amt::make_track(
            start_row,
            .x = x_,
            .y = y_,
            .t = t_,
            crs = crs(env_test)
          )
        }

        start_pt <- start_pt |>
          amt::make_start() |>
          amt::mutate(dt = hours(4))

        kernel <- amt::redistribution_kernel(
          x = issf_train,
          map = env_test,
          fun = function(xy, map) {
            xy %>%
              amt::extract_covariates(map, where = "both") %>%
              amt::time_of_day(include.crepuscule = FALSE, where = "both") %>%
              mutate(
                tod_start_day = as.integer(tod_start_ == "day"),
                tod_start_night = as.integer(tod_start_ == "night"),
                days = lubridate::yday(t2_) - min(lubridate::yday(t2_)) + 1
              )
          },
          start = start_pt,
          landscape = "discrete",
          as.rast = TRUE
        )

        sim_result <- tryCatch(
          amt::simulate_path(kernel, n = n_steps),
          error = function(err) NA
        )

        if (any(is.na(sim_result))) {
          return(NULL)
        }

        sim_burst <- bind_rows(
          sim_burst,
          sim_result %>% select(x_, y_, t_)
        )
      }

      sim_burst %>% mutate(burst_ = b)
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
    "agriculture",
    "grassland",
    "hay",
    "oak",
    "central.hardwoods",
    "other.forest",
    "other"
  )
) {
  # Sort levels longest first to avoid partial matches
  # (e.g., "other" matching inside "other.forest")
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
#' @param deer_data Deer movement data (e.g., deer_mvt)
#' @param model_results List of model results from fitting step
#' @param env Environment raster to pass to simulate_movement
#' @param ndvi_list Named list of NDVI rasters
#' @param n_deer Number of deers to simulate
#' @param n_sim Number of simulations per deer
run_simulations <- function(
  deer_data,
  model_results,
  env,
  ndvi_list,
  n_deer = nrow(deer_data),
  n_sim = 10
) {
  # Pre-crop rasters per deer in main process (small memory footprint)
  deer_crops <- lapply(1:n_deer, function(i) {
    crop_extent <- st_buffer(
      st_as_sf(deer_data$stp_test[[i]], coords = c('x1_', 'y1_'), crs = 6610),
      5000
    )

    ndvi_year <- ndvi_list[[paste("ndvi", deer_data[i, ]$year, sep = "_")]]

    list(
      env = terra::wrap(crop(env, crop_extent)),
      ndvi = terra::wrap(crop(ndvi_year, crop_extent))
    )
  })

  # Pre-extract per-deer data so each worker only receives its own slice
  deer_rows <- lapply(1:n_deer, function(i) deer_data[i, ])
  deer_coeff <- model_results$coeff
  deer_iss <- model_results$iss
  deer_train <- deer_data$stp.var.train

  foreach(
    i = 1:n_deer,
    crop_i = deer_crops,
    deer_i = deer_rows,
    coeff_i = deer_coeff,
    iss_i = deer_iss,
    train_i = deer_train,
    .packages = c("amt", "terra", "sf", "dplyr", "foreach"),
    .export = c("simulate_movement", "rename_landcover_coefs")
  ) %dopar%
    {
      if (length(coeff_i) == 1) {
        return(NA)
      }

      # Each worker only unwraps its own small cropped raster
      env_local <- terra::unwrap(crop_i$env)
      ndvi_local <- terra::unwrap(crop_i$ndvi)

      coefs <- iss_i$model$coefficients
      names(coefs) <- rename_landcover_coefs(names(coefs))
      coefs <- coefs[!is.na(coefs)]

      # Fix interaction term ordering: model.matrix() may order terms in
      # interactions (e.g., cos(ta_):tod_end_day) differently than the
      # coefficient names from the fitted model (e.g., tod_end_day:cos(ta_)).
      # We build a dummy model, check what model.matrix() actually produces
      # using the training data, and swap any mismatched interaction terms.
      dummy_sim <- make_issf_model(coefs = coefs)
      mm_names <- colnames(model.matrix(
        amt:::ssf_formula(dummy_sim$model$formula),
        data = train_i
      ))

      # For each coef name,
      # if it doesn't match mm_names but swapping does, swap it
      for (idx in seq_along(coefs)) {
        if (
          grepl(":", names(coefs)[idx]) && !(names(coefs)[idx] %in% mm_names)
        ) {
          parts <- strsplit(names(coefs)[idx], ":")[[1]]
          perms <- combinat::permn(parts)
          for (p in perms) {
            candidate <- paste(p, collapse = ":")
            if (candidate %in% mm_names) {
              names(coefs)[idx] <- candidate
              break
            }
          }
        }
      }

      model_sim <- make_issf_model(
        coefs = coefs,
        sl = iss_i$sl_,
        ta = iss_i$ta_
      )

      foreach(h = 1:n_sim, .combine = "rbind") %do%
        {
          res <- simulate_movement(
            deer_i,
            env_local,
            ndvi_local,
            model_sim,
            stp_col = "stp_test"
          )
          res$nsim <- h
          res
        }
    }
}

#' Calculate mean Energy Score between observed and simulated paths
#' @param obs Observed path dataframe with x1_, y1_, t1_ columns
#' @param sim Simulated paths dataframe with x_, y_, t_, nsim columns
#' @return Mean energy score across all matched time steps
calc_energy_score <- function(obs, sim) {
  # Round times to nearest hour for matching
  obs_times <- round_date(obs$t1_, unit = "hour")
  sim_times <- round_date(sim$t_, unit = "hour")

  # Find shared time steps
  shared_times <- intersect(obs_times, sim_times)

  if (length(shared_times) == 0) {
    warning("No matching time steps between observed and simulated paths")
    return(NA_real_)
  }

  # Filter to shared times
  obs_matched <- obs[obs_times %in% shared_times, ]
  sim_matched <- sim[sim_times %in% shared_times, ]

  es_per_step <- foreach(t = seq_len(nrow(obs_matched)), .combine = "c") %do%
    {
      obs_xy <- c(obs_matched$x1_[t], obs_matched$y1_[t])
      obs_time <- round_date(obs_matched$t1_[t], unit = "hour")

      # Get all sim locations at this time step
      sim_at_t <- sim_matched[
        round_date(sim_matched$t_, unit = "hour") == obs_time,
      ]

      if (nrow(sim_at_t) == 0) {
        return(NA_real_)
      }

      # es_sample expects: obs = d-vector, dat = d x n matrix
      sim_matrix <- t(cbind(sim_at_t$x_, sim_at_t$y_))
      scoringRules::es_sample(obs_xy, dat = sim_matrix)
    }

  mean(es_per_step, na.rm = TRUE)
}

#' Estimate overlap of utilization distributions
#' @param data Observed paths
#' @param sim Simulated paths
#' @param n_sim number of simulated paths
overlap_ud <- function(data, sim, n_sim) {
  z1 <- data |>
    select("x1_", "y1_", "t1_") |>
    st_as_sf(coords = c("x1_", "y1_"), crs = 6610) |>
    st_transform(4326) |> # Transform to WGS84 lat/long
    mutate(
      longitude = st_coordinates(geometry)[, 1],
      latitude = st_coordinates(geometry)[, 2],
      timestamp = t1_
    ) |>
    st_drop_geometry() |>
    as.data.frame() |>
    as.telemetry()

  z2 <- foreach(k = 1:n_sim) %do%
    {
      tel <- sim |>
        filter(nsim == k) |>
        select("x_", "y_", "t_") |>
        st_as_sf(coords = c("x_", "y_"), crs = 6610) |>
        st_transform(4326) |> # Transform to WGS84 lat/long
        mutate(
          longitude = st_coordinates(geometry)[, 1],
          latitude = st_coordinates(geometry)[, 2],
          timestamp = t_
        ) |>
        st_drop_geometry() |>
        as.data.frame() |>
        as.telemetry()

      projection(tel) <- projection(z1)

      tel
    }

  # Fit single guestimated ctmm model
  ms1 <- ctmm.select(
    z1,
    ctmm.guess(z1, interactive = FALSE),
    verbose = F,
    cores = 1
  )

  ms2 <- map(z2, function(z) {
    tryCatch(
      ctmm.fit(z, ms1),
      error = function(e) NULL
    )
  })

  # Remove failed fits and their corresponding telemetry data
  keep <- !map_lgl(ms2, is.null)
  ms2 <- ms2[keep]
  z2 <- z2[keep]
  ms2_avg <- mean(ms2)

  z1_uds <- akde(z1, ms1)
  z2_uds <- akde(z2, ms2, grid = list(r = z1_uds$r, dr = z1_uds$dr)) %>%
    mean()

  # Estimate the overlap of UDs
  bat_uds <- overlap(list(z1_uds, z2_uds))$CI[1, 2, 2]
  bat_ctmm <- overlap(list(ms1, ms2_avg))$CI[1, 2, 2]

  list(
    bat_uds = bat_uds,
    bat_ctmm = bat_ctmm
  )
}

#' Estimate proximity between observed and a simulated path
#' @param data Observed paths
#' @param sim Simulated paths
#' @param n_sim number of simulated paths
prox_path <- function(data, sim, n_sim) {
  z1 <- data |>
    select("x1_", "y1_", "t1_") |>
    st_as_sf(coords = c("x1_", "y1_"), crs = 6610) |>
    st_transform(4326) |> # Transform to WGS84 lat/long
    mutate(
      longitude = st_coordinates(geometry)[, 1],
      latitude = st_coordinates(geometry)[, 2],
      timestamp = t1_
    ) |>
    st_drop_geometry() |>
    as.data.frame() |>
    as.telemetry()

  z2 <- foreach(k = 1:n_sim) %do%
    {
      tel <- sim |>
        filter(nsim == k) |>
        select("x_", "y_", "t_") |>
        st_as_sf(coords = c("x_", "y_"), crs = 6610) |>
        st_transform(4326) |> # Transform to WGS84 lat/long
        mutate(
          longitude = st_coordinates(geometry)[, 1],
          latitude = st_coordinates(geometry)[, 2],
          timestamp = t_
        ) |>
        st_drop_geometry() |>
        as.data.frame() |>
        as.telemetry()

      projection(tel) <- projection(z1)

      tel
    }

  # Fit single guestiamted ctmm model
  ms1 <- ctmm.select(
    z1,
    ctmm.guess(z1, interactive = FALSE),
    verbose = TRUE,
    cores = 1
  )
  ms1 <- ms1[[1]]

  ms2 <- lapply(z2, function(z) {
    ctmm.select(
      z,
      ctmm.guess(z, interactive = FALSE),
      verbose = TRUE,
      cores = 1
    )
  })
  ms2 <- map(ms2, function(x) x[[1]])

  prox <- foreach(k = 1:n_sim, .combine = "rbind") %do%
    {
      proximity(
        list(z1, z2[[k]]),
        list(ms1, ms2[[k]])
      )
    }

  c(mean(prox[, 2]), length(which(prox[, 3] < 1)))
}
