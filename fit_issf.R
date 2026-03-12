#' @description
#' Simplified deer movement analysis code
#' 1. View raw deer movement data
#' 2. Create random steps with environmental covariates
#' 3. Fit iSSF models
#' 4. Simulate movement from models
#' 5. Estimate and compare utilization distributions
#' 6. Estimate proximity

# Load packages ----------------------------------------------------------------
library(amt)
library(terra)
library(tidyverse)
library(foreach)
library(sf)
library(doParallel)
library(ctmm)

# Load data --------------------------------------------------------------------

# landscape data
env_raster <- rast("env/wiscland/wiscland2_binary.tif")

water_binary <- ifel(env_raster$wiscland == "water", 1, 0)
names(water_binary) <- "Water"

# NDVI data
ndvi_rasters <- list(
  ndvi_2017 = rast('NDVI_2017.tif'),
  ndvi_2018 = rast('NDVI_2018.tif')
)

# Deer movement data
raw_data <- readRDS('Example_code/SW_filtered_deer.RData')

# Prepare sample data
deer_mvt <- raw_data %>%
  filter(year == 2017, sex == 'Female', season == 'fa') %>%
  slice(1:10)

# helper functions
source("helper_functions.R")

# Main workflow ----------------------------------------------------------------

# Set seed for reproducibility
set.seed(1919)

# Step 1: Generate random steps
cat("Generating random steps...\n")
deer_with_random <- make_random_pt_extraction(
  data = deer_mvt,
  n_pts = 10,
  water = water_binary
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
  extract_step_variables(
    data = deer_mvt_random,
    env = env_raster,
    ndvi_list = ndvi_rasters
  )

# Step 3: Fit models
cat("Fitting models...\n")

formulas <-
  c(
    "case_ ~ (log(sl_) + cos(ta_)):tod_end_",
    "case_ ~ (log(sl_) + cos(ta_)):tod_end_ + HR_end",
    "case_ ~ (log(sl_) + cos(ta_)):tod_end_ + HR_end + days",
    "case_ ~ (log(sl_) + cos(ta_)):tod_end_ + HR_end:days",
    "case_ ~ (log(sl_) + cos(ta_)):tod_end_ + (log(sl_) + cos(ta_)):HR_end",
    "case_ ~ (log(sl_) + cos(ta_)):tod_end_:HR_end",
    "case_ ~ (log(sl_) + cos(ta_)):tod_end_:HR_end + days",
    "case_ ~ (log(sl_) + cos(ta_)):tod_end_:HR_end + HR_end:days",
    "case_ ~ (log(sl_) + cos(ta_)):tod_end_ + (log(sl_) + cos(ta_)):HR_end + 
      days",
    "case_ ~ (log(sl_) + cos(ta_)):tod_end_ + (log(sl_) + cos(ta_)):HR_end + 
      days + HR_end:days"
  )

formulas <- paste(formulas, "+ strata(step_id_)")

formula_df <- data.frame(
  formula = formulas,
  name = 1:length(formulas)
)

## Fit models for each formula
model_res <- model_res <- foreach(i = 1:nrow(formula_df)) %do%
  {
    cat("Fitting formula", i, "\n")
    fit_mods(deer_mvt_var, formula_df[i, 1])
  }

# Step 4: Simulate movement
cat("Simulating movement...\n")

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

n_models <- nrow(formula_df)
n_deer <- nrow(deer_mvt)
n_sim <- 10

sim_res <- foreach(m = 1:n_models) %do%
  {
    cat("Simulating model:", formula_df$name[m], "\n")
    run_simulations(
      deer_mvt_var,
      model_res[[m]],
      env_raster,
      ndvi_rasters,
      n_sim = n_sim,
      n_deer = nrow(deer_mvt_var)
    )
  }
names(sim_res) <- formula_df$name

# Step 5: Estimate and compare UD
cat("Estimating overlap of UDs...\n")

res_ud <- foreach(m = 1:n_models) %do%
  {
    cat("UD overlap for model:", formula_df$name[m], "\n")
    foreach(
      i = 1:n_deer,
      .packages = c("amt", "terra", "sf", "tidyverse", "foreach", "ctmm")
    ) %dopar%
      {
        overlap_ud(deer_mvt_var$stp.var[[i]], sim_res[[m]][[i]], n_sim = n_sim)
      }
  }
names(res_ud) <- formula_df$name

# Step 6: Estimate proximity between observed vs simulated paths
cat("Estimating proximity...\n")

prox <- foreach(m = 1:n_models) %do%
  {
    cat("Proximity for model:", formula_df$name[m], "\n")
    foreach(
      i = 1:n_deer,
      .combine = "rbind",
      .packages = "ctmm"
    ) %dopar%
      {
        prox_path(res_ud[[m]][[i]], n_sim = n_sim)
      }
  }
names(prox) <- formula_df$name

stopCluster(cl)

# Plots ------------------------------------------------------------------------
## Visualize results
cat("Creating visualization...\n")
library(patchwork)

plot_list <- foreach(i = c(4, 5, 8, 10)) %do%
  {
    obs_data <- deer_mvt$stp[[i]]
    sim_hr <- sim_res_hr[[i]]
    sim_quad <- sim_res_ndvi[[i]]

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
    p_ndvi <- ggplot() +
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
      labs(title = paste0("Deer ", i, " — NDVI + HR")) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "none"
      )

    list(p_lin, p_ndvi)
  }

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
