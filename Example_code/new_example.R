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
      res <- random_pts %>%
        # operations requires a tibble
        as_tibble() %>%
        filter(case_ == FALSE, Water == 0) %>%
        group_by(step_id_) %>%
        # sample n_pts for each step
        slice_sample(n = n_pts) %>%
        bind_rows(random_pts %>% filter(case_ == TRUE)) %>%
        ungroup() %>%
        select(-Water) %>%
        # go back to the previous structure
        structure(
          class = class(random_pts),
          crs_ = attr(random_pts, "crs_")
        )

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
Make_issf_Models <- function(data, var, formula) {
  data %>%
    mutate(
      iss = map(!!sym(var), ~ fit_mods(.x, formula)),
      tidy_res = map(iss, possibly(~ broom::tidy(.$model), otherwise = NA)),
      AIC = map(iss, possibly(~ AIC(.$model), otherwise = NA))
    )
}

#' Filter models with valid coefficients
#' @param df Dataframe with model results
#' @param season_test Season to filter
#' @param sex_test Sex to filter
#' @param mod_name Model name
#' @param max_coeff Maximum acceptable coefficient
#' @param max_se Maximum acceptable standard error
Good_coefficients <- function(
  df,
  season_test,
  sex_test,
  mod_name,
  max_coeff = 20,
  max_se = 20
) {
  # Process coefficients
  coeff_df <- df %>%
    filter(
      season == season_test,
      sex == sex_test,
      Formula == mod_name,
      !is.na(Coeffs)
    ) %>%
    select(indicator, Coeffs) %>%
    unnest(cols = Coeffs) %>%
    filter(!is.na(estimate)) %>%
    mutate(
      Significant = as.numeric(p.value < 0.005),
      Bad_mean = as.numeric(abs(estimate) > max_coeff | std.error > max_se)
    )

  # Identify good models
  summary_stats <- coeff_df %>%
    group_by(indicator) %>%
    filter(sum(Bad_mean) == 0) %>%
    ungroup() %>%
    group_by(term) %>%
    summarise(
      Number_model = n(),
      Coeff_mean = mean(estimate),
      Coeff_sd = sd(estimate),
      Se_mean = mean(std.error),
      Se_sd = sd(std.error),
      Num_significant = sum(Significant)
    ) %>%
    mutate(
      Max_Coeff = Coeff_mean + 2 * Coeff_sd,
      Min_Coeff = Coeff_mean - 2 * Coeff_sd,
      Max_Se = Se_mean + 2 * Se_sd,
      Min_Se = Se_mean - 2 * Se_sd
    )

  # Find valid models
  valid_models <- coeff_df %>%
    left_join(summary_stats, by = 'term') %>%
    mutate(
      Good_model = case_when(
        estimate < Min_Coeff | estimate > Max_Coeff ~ 0,
        std.error < Min_Se | std.error > Max_Se ~ 0,
        TRUE ~ 1
      )
    ) %>%
    filter(Good_model == 1) %>%
    pull(indicator) %>%
    unique()

  return(list(coeff_df, valid_models))
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

#' Simulate movement path
#' @param GPS_test GPS data for validation
#' @param Simulated_days Date sequence for simulation
#' @param test_ndvi NDVI raster for test year
#' @param nsims Number of simulations
#' @param Env Environmental raster
#' @param issf_train Fitted iSSF model
#' @param y_test Test year
#' @param y_train Training year
#' @param formula_name Formula identifier
Simulate_movement <- function(
  GPS_test,
  Simulated_days,
  test_ndvi,
  nsims,
  Env,
  issf_train,
  y_test,
  y_train,
  formula_name
) {
  Env_crop <- crop(
    Env,
    st_buffer(st_as_sf(GPS_test, coords = c('x1_', 'y1_'), crs = 6610), 5000)
  )

  Simulations <- map(1:nsims, function(sim) {
    simu_i <- GPS_test[1, c('x1_', 'y1_', 't1_')]
    names(simu_i) <- c('x_', 'y_', 't_')

    for (d in seq.Date(
      Simulated_days[1],
      Simulated_days[length(Simulated_days)],
      by = 1
    )) {
      new_Env <- prepare_env_for_day(Env_crop, test_ndvi, d)

      start_pt_simu <- simu_i[nrow(simu_i), ] %>%
        make_track(.x = x_, .y = y_, .t = t_) %>%
        make_start() %>%
        mutate(dt = hours(4))

      kernel <- redistribution_kernel(
        issf_train,
        map = new_Env,
        fun = function(xy, map) extract_covariates(xy, map, where = "both"),
        start = start_pt_simu,
        landscape = "discrete",
        as.rast = TRUE
      )

      Simu_result <- tryCatch(
        list(simulate_path(kernel, n = 6), 'Work'),
        warning = function(w) list(simulate_path(kernel, n = 6), 'Warning'),
        error = function(err) list(NA, 'Error')
      )

      if (!any(is.na(Simu_result[[1]]))) {
        Simu_filtered <- Simu_result[[1]] %>%
          filter(
            as.Date(t_, tz = 'America/Chicago') ==
              as.Date(d, origin = "1970-01-01")
          )
        simu_i <- bind_rows(simu_i, Simu_filtered[, -4])
      }
    }

    simu_i[-1, ] %>%
      mutate(
        sim = sim,
        y_train = y_train,
        y_test = y_test,
        Model = formula_name
      )
  })

  return(list(Simulations, 'Complete'))
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

# Join with original data
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
Formula_df <- data.frame(
  formula = "case_ ~ log(sl_) + cos(ta_) + HR_end + Coarse_scale_end * ndvi_end",
  name = 'Model_test'
)

Model_test <- Deer_mvt_random %>%
  filter(sex == 'Female', season == 'fa') %>%
  select(id, season, year, age.at.col1, sex, indicator, id_s, stp.var)

Model_res <- Model_test %>%
  select(-stp.var) %>%
  slice(rep(1:n(), times = nrow(Formula_df))) %>%
  mutate(Formula = rep(Formula_df$name, each = nrow(Model_test)))

# Fit models for each formula
for (f in 1:nrow(Formula_df)) {
  cat("Fitting formula", f, "\n")
  formula_test <- paste0(Formula_df[f, ]$formula, ' + strata(step_id_)')

  Model <- Make_issf_Models(
    data = Model_test,
    var = 'stp.var',
    formula = formula_test
  )

  Model_res[Model_res$Formula == Formula_df[f, ]$name, ][['Model']] <- Model$iss
  Model_res[Model_res$Formula == Formula_df[f, ]$name, ][[
    'Coeffs'
  ]] <- Model$tidy_res
  Model_res[Model_res$Formula == Formula_df[f, ]$name, ][['AIC']] <- Model$AIC
}

# Step 4: Check coefficient quality
cat("Checking coefficient quality...\n")
GoodCoeff_results <- Good_coefficients(
  df = Model_res,
  sex_test = 'Female',
  season_test = 'fa',
  mod_name = Formula_df[1, 'name']
)

cat("Valid models:", GoodCoeff_results[[2]], "\n")

# Step 5: Simulate movement
cat("Simulating movement...\n")
nsim <- 2
All_sims <- data.frame()

Simulation_deer <- Model_res %>%
  mutate(train = 2017, test = 2018)

for (index in 1:nrow(Simulation_deer)) {
  set.seed(index)

  row_data <- Simulation_deer[index, ]
  Issf <- row_data$Model[[1]]

  if (!any(Issf %in% c('Warning', 'Error'))) {
    GPS_test <- raw_data %>%
      filter(
        id == row_data$id,
        season == row_data$season,
        year == row_data$test
      )

    if (nrow(GPS_test) > 0) {
      # Prepare environment
      Env_test <- Env_rasters[[paste0('Env_', row_data$test)]]
      HR_center <- vect(st_as_sf(
        GPS_test[, c('x_median', 'y_median')],
        coords = c('x_median', 'y_median'),
        crs = crs(Env_test)
      ))
      Env_test$HR <- distance(Env_test$HR, HR_center) / 1000

      # Load NDVI
      test_ndvi <- unwrap(readRDS(paste0(
        "NDVI_Daily_MODIS_250M_",
        row_data$test,
        "_",
        row_data$season,
        ".Rdata"
      )))
      time(test_ndvi) <- as.Date(names(test_ndvi))

      # Reconstruct model
      name_map <- c(
        "log(sl_)" = "log(sl_)",
        "cos(ta_)" = "cos(ta_)",
        "HR_end" = "HR_end",
        "ndvi_end" = 'ndvi_end',
        "Coarse_scale_endagriculture" = "Agriculture_end",
        "Coarse_scale_endgrassland" = "Grassland_end",
        "Coarse_scale_endother" = "Other_end",
        "Coarse_scale_endagriculture:ndvi_end" = "Agriculture_end:ndvi_end",
        "Coarse_scale_endgrassland:ndvi_end" = "Grassland_end:ndvi_end",
        "Coarse_scale_endother:ndvi_end" = "Other_end:ndvi_end"
      )

      coefs <- setNames(
        as.numeric(Issf$model$coefficients),
        name_map[names(Issf$model$coefficients)]
      )
      coefs <- coefs[!is.na(coefs)]

      Model_reconstructed <- make_issf_model(
        coefs = coefs,
        sl = Issf$sl_,
        ta = Issf$ta_
      )

      # Simulate
      Simulated_days <- unique(as.Date(GPS_test$stp[[1]]$t1_, tz = 'CST6CDT'))

      Simulation <- Simulate_movement(
        GPS_test$stp[[1]],
        Simulated_days,
        test_ndvi,
        nsim,
        Env_test,
        Model_reconstructed,
        row_data$test,
        row_data$train,
        row_data$Formula
      )

      Final_sim <- bind_rows(Simulation[[1]]) %>%
        mutate(
          id = row_data$id,
          season = row_data$season,
          Formula = row_data$Formula
        ) %>%
        nest(stp = c(x_, y_, t_))

      All_sims <- bind_rows(All_sims, Final_sim)
    }
  }
}

# Visualize results
cat("Creating visualization...\n")
Mvt_5000_fa <- All_sims[1:2, ] %>%
  select(sim, stp) %>%
  unnest(stp)

Real_mvt <- raw_data %>%
  filter(id == 5000, season == 'fa', year == 2018) %>%
  unnest(stp)

ggplot() +
  geom_spatraster(data = Env_rasters$Env_2018$Coarse_scale) +
  xlim(min(Mvt_5000_fa$x_) - 1000, max(Mvt_5000_fa$x_) + 1000) +
  ylim(min(Mvt_5000_fa$y_) - 1000, max(Mvt_5000_fa$y_) + 1000) +
  geom_path(data = Real_mvt, aes(x = x1_, y = y1_), color = 'black') +
  geom_path(data = Mvt_5000_fa, aes(x = x_, y = y_, color = as.factor(sim))) +
  scale_fill_manual(
    values = c(
      'goldenrod2',
      'forestgreen',
      'lightgreen',
      'grey',
      'blue',
      'white'
    )
  )

cat("Analysis complete!\n")
