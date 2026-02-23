#' @description
#' Simplified deer movement analysis code
#' 1. View raw deer movement data
#' 2. Create random steps with environmental covariates
#' 3. Fit iSSF models
#' 4. Simulate movement from models

# Load packages ----------------------------------------------------------------
library(amt)
library(terra)
library(tidyverse)
library(foreach)
library(sf)


# Load data --------------------------------------------------------------------
env_rasters <- list(
  env_2017 = rast('Example_code/Env_2017.tif'),
  env_2018 = rast('Example_code/Env_2018.tif')
)

# Load NDVI data
ndvi_rasters <- list(
  ndvi_2017 = terra::unwrap(readRDS(
    "Example_code/NDVI_Daily_MODIS_250M_2017_fa.Rdata"
  )),

  ndvi_2018 = terra::unwrap(readRDS(
    "Example_code/NDVI_Daily_MODIS_250M_2018_fa.Rdata"
  ))
)

raw_data <- readRDS('Example_code/SW_filtered_deer.RData')

# Prepare sample data
deer_mvt <- raw_data %>%
  filter(year == 2017, sex == 'Female', season == 'fa') %>%
  slice(1:10)

# Helper functions -------------------------------------------------------------

#' Generate random steps with environmental covariates
#' @param data Dataframe with movement steps
#' @param n_pts Number of random points per step
#' @param Env_list List of environmental rasters by year
make_random_pt_extraction <- function(data, n_pts, env_list = env_rasters) {
  random.stp <- foreach(i = 1:nrow(data)) %do%
    {
      print(i)

      z <- data[i, ]

      # Extract water layer for this year
      water <- env_list[[paste0('env_', z$year)]]$Water
      data_step <- z$stp[[1]]

      # Start with buffer for water removal
      n_random <- ceiling(n_pts * 10)

      # Generate random steps
      random_pts <- data_step %>%
        random_steps(n_control = n_random) %>%
        extract_covariates(water, where = "end")

      # Filter and select final random points
      # Filter and select final random points
      res <- random_pts %>%
        # Select points not on water
        filter(case_ == FALSE, Water == 0) %>%
        # sample n_pts for each step
        slice_sample(n = 10, by = step_id_) %>%
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
  env_list = env_rasters,
  ndvi_list = ndvi_rasters
) {
  stp.var <- foreach(i = 1:nrow(data)) %do%
    {
      data_row <- data[i, ]
      env <- env_list[[paste0('env_', data_row$year)]]

      # Add HR distance to environment
      median_pt <- vect(
        cbind(data_row$x_median, data_row$y_median),
        crs = crs(env)
      )

      ## Any layer in env works for distance calculation
      env$HR <- values(distance(env$ele, median_pt)) / 1000

      daily_ndvi <- ndvi_list[[paste0('ndvi_', data_row$year)]]

      # Process step data
      data_ssf <- data_row$random.stp[[1]] %>%
        extract_covariates(env, where = 'both') %>%
        extract_covariates_var_time(
          daily_ndvi,
          max_time = days(1),
          when = "any",
          where = "both",
          name_covar = "ndvi"
        ) %>%
        # Setting the intercept to the 'forest' land type
        mutate(
          Coarse_scale_end = factor(
            Coarse_scale_end,
            levels = c("forest", 'agriculture', 'grassland', 'other')
          ),
          Coarse_scale_start = factor(
            Coarse_scale_start,
            levels = c("forest", 'agriculture', 'grassland', 'other')
          )
        )
      data_ssf
    }

  data$stp.var <- stp.var

  data
}

#' Fit iSSF model with error handling
#' @param ssf_data Data with cases and covariates
#' @param formula Model formula
fit_mods <- function(ssf_data, formula) {
  tryCatch(
    ssf_data %>% amt::fit_issf(as.formula(formula), model = TRUE),
    warning = function(w) 'Warning',
    error = function(err) 'Error'
  )
}

#' Fit iSSF models to all individuals
#' @param data Dataframe with step variables
#' @param var Column name containing step data
#' @param formula Model formula
fit_multi_mods <- function(data, var, formula) {
  data %>%
    mutate(
      iss = map(!!sym(var), ~ fit_mods(.x, formula)),
      tidy_res = map(iss, possibly(~ broom::tidy(.$model), otherwise = NA)),
      AIC = map(iss, possibly(~ AIC(.$model), otherwise = NA))
    )
}

#' Prepare environment for simulation day
prepare_env_for_day <- function(Env_crop, test_ndvi, date) {
  new_Env <- Env_crop
  new_Env$ndvi <- resample(
    test_ndvi[[which.min(abs(
      as.Date(date, origin = "1970-01-01") - time(test_ndvi)
    ))]],
    new_Env,
    method = 'near'
  )
  return(new_Env)
}

#' Simulate a single movement path
#' @param deer_test Step data for testing sims
#' @param env_test Environmental rasters for test data
#' @param issf_train Fitted iSSF model
#' @param formula_name Formula identifier
simulate_movement <- function(
  deer_test,
  env_test,
  issf_train,
  formula_name
) {
  # Calcualte distance to home range center
  median_pt <- terra::vect(
    cbind(deer_test$x_median, deer_test$y_median),
    crs = crs(env_test[[1]])
  )

  env_test$HR <- NA
  env_test$HR <- terra::distance(env_test$HR, median_pt) / 1000

  env_crop <- crop(
    env_test,
    st_buffer(
      st_as_sf(deer_test$stp[[1]], coords = c('x1_', 'y1_'), crs = 6610),
      5000
    )
  )

  sim_days <- unique(as.Date(deer_test$stp[[1]]$t1_, tz = 'CST6CDT'))

  sim_days_all <- seq.Date(
    first(sim_days),
    last(sim_days),
    by = 1
  )

  # Extract the first step to use as the initial step of the simulation
  sim_i <- deer_test$stp[[1]][1, c('x1_', 'y1_', 't1_')]
  names(sim_i) <- c('x_', 'y_', 't_')

  # Simulate each day one by one. To do this, we simulate 6 steps at a time
  # Because we are on a 4h time scale

  for (d in sim_days_all) {
    # Make the starting pts of the kernel and simulation
    start_pt_sim <- sim_i[nrow(sim_i), ] |>
      make_track(.x = x_, .y = y_, .t = t_) |>
      make_start() |>
      mutate(dt = hours(4))

    # Make the redistribution kernel
    kernel <- redistribution_kernel(
      x = issf_train,
      map = env_test$HR,
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
      list(simulate_path(kernel, n = 6), 'Work'),
      warning = function(w) list(simulate_path(kernel, n = 6), 'Warning'),
      error = function(err) list(NA, 'Error')
    )

    if (!any(is.na(sim_result[[1]]))) {
      sim_filtered <- sim_result[[1]] %>%
        filter(
          as.Date(t_, tz = 'America/Chicago') ==
            as.Date(d, origin = "1970-01-01")
        )
      sim_i <- bind_rows(sim_i, sim_filtered[, -4])
    }
  }

  sim_i[-1, ]
}

# Main workflow ----------------------------------------------------------------

# Set seed for reproducibility
set.seed(1234)

# Step 1: Generate random steps
cat("Generating random steps...\n")
deer_with_random <- make_random_pt_extraction(
  data = deer_mvt,
  n_pts = 10,
  env_list = env_rasters
) %>%
  select(id, season, year, age.at.col1, sex, indicator, random.stp)

## Join with original data
deer_mvt_random <- deer_mvt %>%
  left_join(
    deer_with_random,
    by = c('id', 'season', 'year', 'age.at.col1', 'sex', 'indicator')
  )

# Step 2: Extract environmental variables
cat("Extracting environmental variables...\n")
deer_mvt_var <-
  extract_step_variables(data = deer_mvt_random)

# Step 3: Fit models
cat("Fitting models...\n")
formula_df <- data.frame(
  formula = c(
    "case_ ~ log(sl_) + cos(ta_) + HR_end + strata(step_id_)",
    "case_ ~ log(sl_) + cos(ta_) + HR_end + I(HR_end^2) + strata(step_id_)"
  ),
  name = c("null_model_lin", "null_model_quad")
)

## Fit models for each formula
model_res <- foreach(i = 1:nrow(formula_df)) %do%
  {
    cat("Fitting formula", i, "\n")

    model <- fit_multi_mods(
      data = deer_mvt_var,
      var = 'stp.var',
      formula = formula_df[i, 1]
    )

    res <- list(
      iss = model$iss,
      coeff = model$tidy_res,
      aic = unlist(model$AIC)
    )

    res
  }

# Step 4: Simulate movement
cat("Simulating movement...\n")

sim_res_hr <- foreach(i = 1:10) %do%
  {
    coefs <- model_res[[1]]$iss[[i]]$model$coefficients
    coefs <- coefs[!is.na(coefs)]

    model_sim <- make_issf_model(
      coefs = coefs,
      sl = model_res[[1]]$iss[[i]]$sl_,
      ta = model_res[[1]]$iss[[i]]$ta_
    )

    foreach(h = 1:10, .combine = "rbind") %do%
      {
        res <- simulate_movement(
          deer_mvt[i, ],
          env_rasters[[1]],
          model_sim,
          formula_df[1, 1]
        )

        res$nsim <- h

        res
      }
  }

sim_res_hr_quad <- foreach(i = 1:10) %do%
  {
    coefs <- model_res[[2]]$iss[[i]]$model$coefficients
    coefs <- coefs[!is.na(coefs)]

    model_sim <- make_issf_model(
      coefs = coefs,
      sl = model_res[[2]]$iss[[i]]$sl_,
      ta = model_res[[2]]$iss[[i]]$ta_
    )

    foreach(h = 1:10, .combine = "rbind") %do%
      {
        res <- simulate_movement(
          deer_mvt[i, ],
          env_rasters[[1]],
          model_sim,
          formula_df[2, 1]
        )

        res$nsim <- h

        res
      }
  }

## Visualize results
cat("Creating visualization...\n")
library(patchwork)

plot_list <- foreach(i = 1:4) %do%
  {
    obs_data <- deer_mvt$stp[[i]]
    sim_hr <- sim_res_hr[[i]]
    sim_quad <- sim_res_hr_quad[[i]]

    x_range_lin <- range(c(obs_data$x1_, sim_hr$x_))
    y_range_lin <- range(c(obs_data$y1_, sim_hr$y_))

    x_range_quad <- range(c(obs_data$x1_, sim_quad$x_))
    y_range_quad <- range(c(obs_data$y1_, sim_quad$y_))

    # Linear HR model
    p_lin <- ggplot() +
      geom_spatraster(data = env_rasters$env_2017$Coarse_scale) +
      xlim(x_range_lin[1] - 1000, x_range_lin[2] + 1000) +
      ylim(y_range_lin[1] - 1000, y_range_lin[2] + 1000) +
      geom_path(
        data = sim_hr,
        aes(x = x_, y = y_, group = nsim),
        color = "black",
        linewidth = 0.3,
        alpha = 0.3
      ) +
      geom_path(
        data = obs_data,
        aes(x = x1_, y = y1_),
        color = 'white',
        linewidth = 0.8,
        alpha = 0.5
      ) +
      labs(title = paste0("Deer ", i, " — Linear HR")) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none"
      )

    # Quadratic HR model
    p_quad <- ggplot() +
      geom_spatraster(data = env_rasters$env_2017$Coarse_scale) +
      xlim(x_range_quad[1] - 1000, x_range_quad[2] + 1000) +
      ylim(y_range_quad[1] - 1000, y_range_quad[2] + 1000) +
      geom_path(
        data = sim_quad,
        aes(x = x_, y = y_, group = nsim),
        color = "black",
        linewidth = 0.3,
        alpha = 0.3
      ) +
      geom_path(
        data = obs_data,
        aes(x = x1_, y = y1_),
        color = 'white',
        linewidth = 0.8,
        alpha = 0.5
      ) +
      labs(title = paste0("Deer ", i, " — Quadratic HR")) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none"
      )

    list(p_lin, p_quad)
  }

combined_plot <- wrap_plots(
  unlist(plot_list, recursive = FALSE),
  ncol = 2,
  nrow = 4
)

ggsave(
  "Example_code/deer_movement_comparison.png",
  plot = combined_plot,
  width = 10,
  height = 16,
  dpi = 300
)

# Single deer plot
ggplot() +
  geom_spatraster(data = env_rasters$env_2017$Coarse_scale) +
  xlim(
    min(sim_res_hr_quad[[3]]$x_) - 1000,
    max(sim_res_hr_quad[[3]]$x_) + 1000
  ) +
  ylim(
    min(sim_res_hr_quad[[3]]$y_) - 1000,
    max(sim_res_hr_quad[[3]]$y_) + 1000
  ) +
  geom_path(
    data = deer_mvt$stp[[3]],
    aes(x = x1_, y = y1_),
    color = 'black',
    alpha = 0.5
  ) +
  geom_path(
    data = sim_res_hr_quad[[3]],
    aes(x = x_, y = y_, group = nsim),
    color = "orange",
    alpha = 0.1
  )

# Calculate overlap of utilization distributions

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

  z2 <- foreach(i = 1:n_sim) %do%
    {
      tel <- sim |>
        filter(nsim == i) |>
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

  ms1 <- ctmm.select(z1, ctmm.guess(z1), verbose = TRUE, cores = 4)

  ms2 <- foreach(i = 1:length(z2)) %do%
    {
      ctmm.select(z2[[i]], ctmm.guess(z2[[i]]), verbose = TRUE, cores = 4)
    }

  uds <- foreach(i = 1:length(z2)) %do%
    {
      akde(list(z1, z2[[i]]), list(ms1[[1]], ms2[[i]][[1]]))
    }

  bat <- map_dbl(uds, function(x) overlap(x)$CI[1, 2, 2])

  list(
    data = z1,
    sim = z2,
    ctmm_data = ms1,
    ctmm_sim = ms2,
    uds = uds,
    bat_coeff = bat
  )
}

res_ud_hr <- foreach(i = 1:10) %do%
  {
    overlap_ud(deer_mvt$stp[[i]], sim_res_hr[[i]], n_sim = 10)
  }

res_ud_hr_quad <- foreach(i = 1:10) %do%
  {
    overlap_ud(deer_mvt$stp[[i]], sim_res_hr_quad[[i]], n_sim = 10)
  }
