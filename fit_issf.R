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
env_old <- rast("Example_code/Env_2017.tif")
env_raster <- terra::crop(env_raster, env_old) %>%
  terra::resample(env_old)

env_raster$ele <- env_old$ele
env_raster$east <- env_old$eastness
env_raster$east2 <- env_old$eastness2
env_raster$north <- env_old$northness
env_raster$dist <- env_old$fe.dist

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
  ) %>%
  slice(1:30)

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

saveRDS(deer_mvt, "data_deer.rds")

# Step 3: Fit models
cat("Fitting models...\n")

formulas <-
  c(
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + HR_end",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + HR_end:days",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + (log(sl_) + cos(ta_)):HR_start",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_:HR_start",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_:HR_start + HR_end:days",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + (log(sl_) + cos(ta_)):HR_start",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + (log(sl_) + cos(ta_)):HR_start + 
      HR_end:days",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + HR_end + east_end + east2_end + 
      north_end + dist_end + wiscland_end + wiscland_end:ndvi_end",
    "case_ ~ (log(sl_) + cos(ta_)):tod_start_ + east_end + east2_end + 
      north_end + dist_end + wiscland_end + wiscland_end:ndvi_end"
  )

formulas <- paste(formulas, "+ strata(step_id_)")

formula_df <- data.frame(
  formula = formulas,
  name = 1:length(formulas)
)

## Fit models for each formula
results <- list()

results$models <- foreach(i = 1:nrow(formula_df)) %do%
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

results$sim <- foreach(m = 1:n_models) %do%
  {
    cat("Simulating model:", formula_df$name[m], "\n")
    run_simulations(
      deer_mvt,
      results$models[[m]],
      env_raster,
      ndvi_rasters,
      n_sim = n_sim,
      n_deer = nrow(deer_mvt)
    )
  }
names(results$sim) <- formula_df$name

# Step 5: Estimate UD overlap
cat("Estimating overlap of UDs and CTMMs...\n")

results$ud <- foreach(m = 1:n_models) %do%
  {
    cat("UD and CTMM overlap for model:", formula_df$name[m], "\n")
    foreach(
      i = 1:n_deer,
      .packages = c("amt", "terra", "sf", "tidyverse", "foreach", "ctmm")
    ) %dopar%
      {
        overlap_ud(
          deer_mvt$stp_test[[i]],
          results$sim[[m]][[i]],
          n_sim = n_sim
        )
      }
  }
names(results$ud) <- formula_df$name

# Step 6: Calculate Energy Scores
cat("Calculating Energy Scores...\n")

results$es <- foreach(m = 1:n_models, .combine = "rbind") %do%
  {
    cat("Energy Score for model:", formula_df$name[m], "\n")

    scores <- foreach(
      i = 1:n_deer,
      .combine = "c"
    ) %do%
      {
        sim_i <- results$sim[[m]][[i]]

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

results$es <- results$es %>%
  group_by(deer) %>%
  mutate(energy_skill = 1 - (energy_score / energy_score[2]))

# Step 7: Select best models for each deer
cat("Selecting best models...\n")

# Wrangle res_ud into a data frame
res_bat <- foreach(m = names(results$ud), .combine = "rbind") %do%
  {
    scores_uds <- foreach(i = 1:n_deer, .combine = "c") %do%
      {
        if (length(results$ud[[m]][[i]]) == 1 && is.na(results$ud[[m]][[i]])) {
          return(NA_real_)
        }
        results$ud[[m]][[i]]$bat_uds
      }

    scores_ctmm <- foreach(i = 1:n_deer, .combine = "c") %do%
      {
        if (length(results$ud[[m]][[i]]) == 1 && is.na(results$ud[[m]][[i]])) {
          return(NA_real_)
        }
        results$ud[[m]][[i]]$bat_ctmm
      }

    data.frame(
      model = m,
      deer = deer_mvt$id,
      bat_uds = scores_uds,
      bat_ctmm = scores_ctmm
    )
  }

res_bat$model <- factor(res_bat$model, levels = unique(res_bat$model))
results$es$model <- factor(results$es$model, levels = unique(results$es$model))

# Filter and select best model per deer
model_selection <- res_bat %>%
  left_join(results$es, by = c("deer", "model")) %>%
  group_by(deer) %>%
  mutate(
    step1 = bat_uds > 0.8,
    step2 = step1 & bat_ctmm > 0.8,
    step3 = step2 & energy_skill == max(energy_skill[step2])
  ) %>%
  ungroup()

selected <- model_selection %>%
  filter(step3) %>%
  select(deer, model, bat_uds, bat_ctmm, energy_skill)

no_selection <- model_selection %>%
  group_by(deer) %>%
  summarize(
    passed_step1 = sum(step1),
    passed_step2 = sum(step2),
    passed_step3 = sum(step3)
  ) %>%
  filter(passed_step3 == 0)

# Step 8: Calculate proximity for selected models --------------------------------
cat("Calculating proximity for selected models...\n")

results$prox <- foreach(
  r = 1:nrow(selected),
  .packages = c("amt", "ctmm", "tidyverse", "sf", "foreach"),
  .combine = "rbind"
) %dopar%
  {
    deer_idx <- which(deer_mvt$id == selected$deer[r])
    model_name <- as.character(selected$model[r])

    obs <- deer_mvt$stp_test[[deer_idx]]
    sim <- results$sim[[model_name]][[deer_idx]]

    n_sim <- length(unique(sim$nsim))

    prox_result <- prox_path(data = obs, sim = sim, n_sim = n_sim)

    data.frame(
      deer = selected$deer[r],
      model = model_name,
      mean_prox = prox_result[1],
      prox_lt1 = prox_result[2]
    )
  }

stopCluster(cl)

# Plots ------------------------------------------------------------------------

library(patchwork)


# Energy score violin plot
es <-
  p_es <- ggplot(
    es,
    aes(x = as.factor(model), y = energy_score)
  ) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
    geom_violin(fill = "lightblue", alpha = 0.5) +
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

best_models <- res_bat %>%
  group_by(deer) %>%
  summarise(
    best_m_uds = which.max(bat_uds),
    best_m_ctmm = which.max(bat_ctmm),
    best_uds = max(bat_uds),
    best_ctmm = max(bat_ctmm)
  )

best_models <- left_join(
  best_models,
  (es %>%
    group_by(deer) %>%
    summarise(
      best_m_es = which.max(energy_score),
      best_es = max(energy_score)
    )),
  by = "deer"
)

best_models$best_m_uds %>% table()
best_models$best_m_ctmm %>% table()
best_models$best_m_es %>% table()

# Single deer plot
ggplot() +
  geom_path(
    data = deer_mvt$stp_test[[14]],
    aes(x = x1_, y = y1_),
    color = 'black',
    alpha = 0.5
  ) +
  geom_path(
    data = filter(results$sim[[11]][[14]], nsim == 7),
    aes(x = x_, y = y_),
    color = "orange",
    alpha = 0.5
  ) +
  theme_minimal()

ggplot() +
  geom_path(
    data = deer_mvt$stp_test[[14]],
    aes(x = x1_, y = y1_),
    color = 'black',
    alpha = 0.5
  ) +
  geom_path(
    data = results$sim[[11]][[14]],
    aes(x = x_, y = y_, group = n_sim),
    color = "orange",
    alpha = 0.2
  ) +
  theme_minimal()

ggplot() +
  geom_path(
    data = deer_mvt$stp_test[[4]],
    aes(x = x1_, y = y1_),
    color = 'black',
    alpha = 0.5
  ) +
  geom_path(
    data = results$sim[[2]][[4]],
    aes(x = x_, y = y_, group = n_sim),
    color = "orange",
    alpha = 0.2
  ) +
  theme_minimal()
