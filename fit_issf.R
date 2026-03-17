#' @description
#' Simplified deer movement analysis code
#' 1. View raw deer movement data
#' 2. Create random steps with environmental covariates
#' 3. Fit iSSF models
#' 4. Simulate movement from models
#' 5. Estimate and compare utilization distributions
#' 6. Calculate Energy Scores

# Load packages ----------------------------------------------------------------
library(amt)
library(terra)
library(tidyverse)
library(foreach)
library(sf)
library(doParallel)
library(ctmm)
library(scoringRules)

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
  filter(year %in% c(2017, 2018))

# Separate training and test splits
deer_mvt <- deer_mvt %>%
  mutate(
    stp_train = map(stp, ~ .x[1:ceiling(nrow(.x) * 0.7), ]),
    stp_test = map(stp, ~ .x[(ceiling(nrow(.x) * 0.7) + 1):nrow(.x), ])
  )

# helper functions
source("helper_functions.R")

# Main workflow ----------------------------------------------------------------

# Set seed for reproducibility
set.seed(1919)

# Step 1: Generate random steps
cat("Generating random steps for train...\n")
deer_mvt <- make_random_pt_extraction(
  data = deer_mvt,
  n_pts = 10,
  water = water_binary,
  stp_col = "stp_train",
  output_col = "stp.random.train"
)

cat("Generating random steps for test...\n")
deer_mvt <- make_random_pt_extraction(
  data = deer_mvt,
  n_pts = 10,
  water = water_binary,
  stp_col = "stp_test",
  output_col = "stp.random.test"
)

# Step 2: Extract environmental variables
cat("Extracting environmental variables for train...\n")
deer_mvt <- extract_step_variables(
  data = deer_mvt,
  env = env_raster,
  ndvi_list = ndvi_rasters,
  random_col = "stp.random.train",
  output_col = "stp.var.train"
)

cat("Extracting environmental variables for test...\n")
deer_mvt <- extract_step_variables(
  data = deer_mvt,
  env = env_raster,
  ndvi_list = ndvi_rasters,
  random_col = "stp.random.test",
  output_col = "stp.var.test"
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
res_model <- foreach(i = 1:nrow(formula_df)) %do%
  {
    cat("Fitting formula", i, "\n")
    fit_mods(deer_mvt, formula_df[i, 1])
  }

# Step 4: Simulate movement
cat("Simulating movement...\n")

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

n_models <- nrow(formula_df)
n_deer <- nrow(deer_mvt)
n_sim <- 10

res_sim <- foreach(m = 1:n_models) %do%
  {
    cat("Simulating model:", formula_df$name[m], "\n")
    run_simulations(
      deer_mvt,
      model_res[[m]],
      env_raster,
      ndvi_rasters,
      n_sim = n_sim,
      n_deer = nrow(deer_mvt)
    )
  }
names(res_sim) <- formula_df$name

saveRDS(res_sim, "results_sim.rds")

# Step 5: Estimate UD overlap
cat("Estimating overlap of UDs...\n")

res_ud <- foreach(m = 1:n_models) %do%
  {
    cat("UD overlap for model:", formula_df$name[m], "\n")
    foreach(
      i = 1:n_deer,
      .packages = c("amt", "terra", "sf", "tidyverse", "foreach", "ctmm")
    ) %dopar%
      {
        overlap_ud(
          deer_mvt$stp_test[[i]],
          res_sim[[m]][[i]],
          n_sim = n_sim
        )
        #prox <- prox_path(ud, n_sim = n_sim)

        #list(bat = ud$bat_coeff, prox = prox)
      }
  }
names(res_ud) <- formula_df$name

saveRDS(res_ud, "results_ud.rds")

# Step 6: Calculate Energy Scores
cat("Calculating Energy Scores...\n")

res_es <- foreach(m = 1:n_models, .combine = "rbind") %do%
  {
    cat("Energy Score for model:", formula_df$name[m], "\n")

    scores <- foreach(
      i = 1:n_deer,
      .combine = "c"
    ) %do%
      {
        sim_i <- res_sim[[m]][[i]]

        # Skip if simulation failed
        if (length(sim_i) == 1 && is.na(sim_i)) {
          return(NA_real_)
        }

        calc_energy_score(
          obs = deer_mvt$stp_test[[i]],
          sim = sim_i
        )
      }

    data.frame(
      model = formula_df$name[m],
      deer = deer_mvt$id,
      energy_score = scores
    )
  }

saveRDS(res_es, "results_es.rds")

stopCluster(cl)

# Plots ------------------------------------------------------------------------

library(patchwork)

# Wrangle res_ud into a data frame
res_bat <- foreach(m = names(res_ud), .combine = "rbind") %do%
  {
    scores_uds <- foreach(i = 1:n_deer, .combine = "c") %do%
      {
        if (length(res_ud[[m]][[i]]) == 1 && is.na(res_ud[[m]][[i]])) {
          return(NA_real_)
        }
        res_ud[[m]][[i]]$bat_uds
      }

    scores_ctmm <- foreach(i = 1:n_deer, .combine = "c") %do%
      {
        if (length(res_ud[[m]][[i]]) == 1 && is.na(res_ud[[m]][[i]])) {
          return(NA_real_)
        }
        res_ud[[m]][[i]]$bat_ctmm$CI[1, 2, 2]
      }

    data.frame(
      model = m,
      deer = deer_mvt$id,
      bat_uds = scores_uds,
      bat_ctmm = scores_ctmm
    )
  }

# Energy score violin plot
p_es <- ggplot(
  res_es,
  aes(x = as.factor(model), y = 1 - (energy_score / 500))
) +
  geom_violin(trim = FALSE, fill = "lightblue", alpha = 0.5) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  labs(
    x = "Model",
    y = "Energy Skill Score"
  ) +
  theme_minimal()

# Bhattacharyya coefficient violin plot
p_uds <- ggplot(res_bat, aes(x = as.factor(model), y = bat_uds)) +
  geom_violin(trim = T, fill = "lightgreen", alpha = 0.5) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  labs(
    x = "Model",
    y = "Bhattacharyya Coefficient (UD)"
  ) +
  theme_minimal()

p_ctmm <- ggplot(res_bat, aes(x = as.factor(model), y = bat_ctmm)) +
  geom_violin(trim = T, fill = "lightgreen", alpha = 0.5) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  labs(
    x = "Model",
    y = "Bhattacharyya Coefficient (CTMM)"
  ) +
  theme_minimal()

p_uds
p_ctmm
p_es


cat("Creating visualization...\n")


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
