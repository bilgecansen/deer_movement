# Helper functions -------------------------------------------------------------

#' Generate random steps with environmental covariates
#' @param data Dataframe with movement steps
#' @param n_pts Number of random points per step
#' @param water binary raster for water bodies
make_random_pt_extraction <- function(data, n_pts, water = water_binary) {
  random.stp <- foreach(i = 1:nrow(data)) %do%
    {
      print(i)

      z <- data[i, ]

      data_step <- z$stp[[1]]

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

  data$random.stp <- random.stp
  data
}

#' Extract environmental variables for each step
#' @param data Dataframe with random steps
#' @param env_list Environmental rasters
#' @param ndvi_list NDVI rasters
extract_step_variables <- function(
  data,
  env = env_raster$wiscland,
  ndvi_list = ndvi_rasters
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
      data_ssf <- data_row$random.stp[[1]] %>%
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
        amt::time_of_day(include.crepuscule = F) %>%
        mutate(
          tod_end_day = as.integer(tod_end_ == "day"),
          tod_end_night = as.integer(tod_end_ == "night"),
          days = lubridate::yday(t2_) - min(lubridate::yday(t2_)) + 1
        )
      data_ssf
    }

  data$stp.var <- stp.var

  data
}

#' Fit iSSF models to all individuals
#' @param data Dataframe with step variables in stp.var column
#' @param formula Model formula string
#' @return List with iss (fitted models), coeff (tidy coefficients), aic (AIC values)
fit_mods <- function(data, formula) {
  iss <- map(data$stp.var, function(ssf_data) {
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
  issf_train
) {
  # Calculate distance to home range center
  median_pt <- terra::vect(
    cbind(deer_data$x_median, deer_data$y_median),
    crs = crs(env_test[[1]])
  )

  env_test$HR <- NA
  env_test$HR <- terra::distance(env_test$HR, median_pt) / 1000

  sim_days <- unique(as.Date(deer_data$stp[[1]]$t1_, tz = 'CST6CDT'))

  sim_days_all <- seq.Date(
    first(sim_days),
    last(sim_days),
    by = 1
  )

  # Extract the first step to use as the initial step of the simulation
  sim_i <- deer_data$stp[[1]][1, c('x1_', 'y1_', 't1_')]
  names(sim_i) <- c('x_', 'y_', 't_')

  # Simulate each day one by one. To do this, we simulate 6 steps at a time
  # Because we are on a 4h time scale

  for (d in sim_days_all) {
    # Get the corresponding ndvi value for that day
    env_test$ndvi <- resample(
      ndvi_test[[which.min(abs(
        as.Date(d, origin = "1970-01-01") - terra::time(ndvi_test)
      ))]],
      env_test,
      method = 'near'
    )

    # Make the starting pts of the kernel and simulation
    start_pt_sim <- sim_i[nrow(sim_i), ] |>
      amt::make_track(.x = x_, .y = y_, .t = t_, crs = crs(env_test)) |>
      amt::make_start() |>
      amt::mutate(dt = hours(4))

    # Make the redistribution kernel
    kernel <- amt::redistribution_kernel(
      x = issf_train,
      map = env_test,
      fun = function(xy, map) {
        xy %>%
          amt::extract_covariates(map, where = "both") %>%
          amt::time_of_day(include.crepuscule = FALSE) %>%
          mutate(
            tod_end_day = as.integer(tod_end_ == "day"),
            tod_end_night = as.integer(tod_end_ == "night"),
            days = lubridate::yday(t2_) - min(lubridate::yday(t2_)) + 1
          )
      },
      start = start_pt_sim,
      landscape = "discrete",
      as.rast = TRUE
    )

    # Simulate the mvt using the redistribution kernel and the starting pt
    # if simulation fails, we have a NA

    # CHANGE TIMING TO ENSURE WE SIMULATE ON THE RIGHT DATE AND TIME !!!
    ## For each day, the 1st point should always be the begining of the day !
    ### Ensure the sim dates are all on the same day, remove it if not

    sim_result <- tryCatch(
      amt::simulate_path(kernel, n = 6),
      error = function(err) NA
    )

    if (!any(is.na(sim_result))) {
      sim_filtered <- sim_result %>%
        filter(
          as.Date(t_, tz = 'America/Chicago') ==
            as.Date(d, origin = "1970-01-01")
        )
      sim_i <- bind_rows(sim_i, sim_filtered[, -4])
    }
  }

  sim_i[-1, ]
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
      st_as_sf(deer_data$stp[[i]], coords = c('x1_', 'y1_'), crs = 6610),
      5000
    )

    ndvi_year <- ndvi_list[[paste("ndvi", deer_data[i, ]$year, sep = "_")]]

    list(
      env = terra::wrap(crop(env, crop_extent)),
      ndvi = terra::wrap(crop(ndvi_year, crop_extent))
    )
  })

  foreach(
    i = 1:n_deer,
    .packages = c("amt", "terra", "sf", "dplyr", "foreach"),
    .export = c("simulate_movement", "rename_landcover_coefs")
  ) %dopar%
    {
      if (length(model_results$coeff[[i]]) == 1) {
        return(NA)
      }

      # Each worker only unwraps a small cropped raster
      env_local <- terra::unwrap(deer_crops[[i]]$env)
      ndvi_local <- terra::unwrap(deer_crops[[i]]$ndvi)

      coefs <- model_results$iss[[i]]$model$coefficients
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
        data = deer_data$stp.var[[i]]
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
        sl = model_results$iss[[i]]$sl_,
        ta = model_results$iss[[i]]$ta_
      )

      foreach(h = 1:n_sim, .combine = "rbind") %do%
        {
          res <- simulate_movement(
            deer_data[i, ],
            env_local,
            ndvi_local,
            model_sim
          )
          res$nsim <- h
          res
        }
    }
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

  # Fit single guestiamted ctmm model
  ms1 <- ctmm.select(
    z1,
    ctmm.guess(z1, interactive = FALSE),
    verbose = TRUE,
    cores = 1
  )
  ms2 <- lapply(z2, function(z) {
    ctmm.select(
      z,
      ctmm.guess(z, interactive = FALSE),
      verbose = TRUE,
      cores = 1
    )
  })

  uds <- foreach(k = 1:length(z2)) %do%
    {
      akde(list(z1, z2[[k]]), list(ms1[[1]], ms2[[k]][[1]]))
    }

  # Estimate the overlap of UDs
  bat <- map_dbl(uds, function(x) overlap(x)$CI[1, 2, 2])

  list(
    data = z1,
    sim = z2,
    ctmm_data = ms1[[1]],
    ctmm_sim = ms2,
    uds = uds,
    bat_coeff = bat
  )
}

#' Estimate proximity between observed and a simulated path
#' @param res_ud output of overlap_ud function
#' @param n_sim Number of simulated paths
prox_path <- function(res_ud, n_sim) {
  z1 <- res_ud$data
  z2 <- res_ud$sim

  prox <- foreach(k = 1:n_sim, .combine = "rbind") %do%
    {
      proximity(
        list(z1, z2[[k]]),
        list(res_ud$ctmm_data, res_ud$ctmm_sim[[k]][[1]])
      )
    }

  c(mean(prox[, 2]), length(which(prox[, 3] < 1)))
}
