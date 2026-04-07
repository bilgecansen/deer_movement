library(foreach)
library(tidyverse)

# Load and combine per-deer skill results -------------------------------------

result_files <- list.files(
  "results",
  pattern = "results_skill_\\d+\\.rds",
  full.names = TRUE
)

# Sort by deer number
deer_numbers <- as.integer(gsub(
  ".*results_skill_(\\d+)\\.rds",
  "\\1",
  result_files
))
result_files <- result_files[order(deer_numbers)]
deer_numbers <- sort(deer_numbers)

cat(sprintf("Loading %d deer skill results...\n", length(result_files)))

results <- purrr::map(result_files, readRDS)

n_deer <- length(results)
n_models <- length(results[[1]]$ud)

# Load AIC from iSSF models and compute delta AIC ------------------------------
results_issf <- readRDS("results_issf_1_119.rds")
n_models <- length(results_issf)

# Build delta_aic data frame for the deer that have skill results
# deer_numbers contains the row indices of the 116 deer with results
res_aic <- purrr::map_dfr(1:n_models, function(m) {
  aic_all <- results_issf[[m]]$aic
  aic_null <- results_issf[[2]]$aic # null model (model 2)

  # delta_aic = aic_null - aic_focal; positive means focal is better
  delta <- aic_null - aic_all
  data.frame(
    model = m,
    deer = deer_numbers,
    delta_aic = delta[deer_numbers]
  )
})

# Wrangle UD results into a data frame -----------------------------------------
res_bat <- foreach(i = 1:n_deer, .combine = "rbind") %do%
  {
    scores_uds <- purrr::map_dbl(1:n_models, function(m) {
      ud <- results[[i]]$ud[[m]]
      if (length(ud) == 1 && is.na(ud)) {
        return(NA_real_)
      }
      ud$bat_uds
    })

    scores_ctmm <- purrr::map_dbl(1:n_models, function(m) {
      ud <- results[[i]]$ud[[m]]
      if (length(ud) == 1 && is.na(ud)) {
        return(NA_real_)
      }
      ud$bat_ctmm
    })

    data.frame(
      model = 1:n_models,
      deer = deer_numbers[i],
      bat_uds = scores_uds,
      bat_ctmm = scores_ctmm
    )
  }

# Combine energy scores from all deer and compute energy skill
res_es <- purrr::map2_dfr(results, deer_numbers, function(r, d) {
  r$es %>% dplyr::mutate(deer = d)
}) %>%
  dplyr::group_by(deer) %>%
  dplyr::mutate(energy_skill = 1 - (energy_score / energy_score[2])) %>%
  dplyr::ungroup()


# Ensure model columns match types before joining
res_bat$model <- as.integer(res_bat$model)
res_es$model <- as.integer(res_es$model)

# Filter and select best model per deer
model_selection <- res_bat %>%
  left_join(res_es, by = c("deer", "model")) %>%
  left_join(res_aic, by = c("deer", "model")) %>%
  group_by(deer) %>%
  mutate(
    step1 = bat_uds >= 0.8,
    step2 = step1 & bat_ctmm >= 0.8,
    step3 = step2 & (delta_aic > 2 | delta_aic == 0),
    step4 = step3 & (energy_skill >= 0)
  ) %>%
  ungroup()

selected <- model_selection %>%
  filter(step4) %>%
  group_by(deer) %>%
  filter(energy_skill == max(energy_skill)) %>%
  ungroup() %>%
  select(deer, model, bat_uds, bat_ctmm, delta_aic, energy_skill)

no_selection <- model_selection %>%
  group_by(deer) %>%
  summarize(
    passed_step1 = sum(step1),
    passed_step2 = sum(step2),
    passed_step3 = sum(step3)
  ) %>%
  filter(passed_step3 == 0)


# Plots ------------------------------------------------------------------------

library(patchwork)

# Home range overlap (BC) violin plot
p_uds <- ggplot(res_bat, aes(x = as.factor(model), y = bat_uds)) +
  geom_violin(trim = T, fill = "#FF644E", alpha = 0.5) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  labs(
    x = "Model",
    y = "Bhattacharyya Coefficient (UD)"
  ) +
  geom_hline(
    yintercept = 0.8,
    col = "red",
    alpha = 0.8,
    linetype = 2,
    linewidth = 2
  ) +
  theme_minimal()

# CTMM overlap (BC) violin plot
p_ctmm <- ggplot(
  res_bat[model_selection$step1 == T, ],
  aes(x = as.factor(model), y = bat_ctmm)
) +
  geom_violin(trim = T, fill = "#16E7CF", alpha = 0.5) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  labs(
    x = "Model",
    y = "Bhattacharyya Coefficient (CTMM)"
  ) +
  geom_hline(
    yintercept = 0.8,
    col = "red",
    alpha = 0.8,
    linetype = 2,
    linewidth = 2
  ) +
  theme_minimal()

# Delta AIC violin plot
idx_aic <- which(model_selection$step1 == T & model_selection$step2 == T)
p_aic <- ggplot(
  res_aic[idx_aic, ],
  aes(x = as.factor(model), y = delta_aic)
) +
  geom_violin(trim = T, fill = "#00A2FF", alpha = 0.5) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  labs(
    x = "Model",
    y = "Delta AIC (relative to null)"
  ) +
  geom_hline(
    yintercept = 2,
    col = "red",
    alpha = 0.8,
    linetype = 2,
    linewidth = 2
  ) +
  theme_minimal()

# Energy score violin plot
idx_es <- which(
  model_selection$step1 == T &
    model_selection$step2 == T &
    model_selection$step3 == T
)
p_es <- ggplot(
  res_es[idx_es, ],
  aes(x = as.factor(model), y = energy_skill)
) +
  geom_violin(fill = "#61D836", alpha = 0.5) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  geom_hline(
    yintercept = 0,
    col = "red",
    alpha = 0.8,
    linetype = 2,
    linewidth = 2
  ) +
  labs(
    x = "Model",
    y = "Energy Skill Score"
  ) +
  theme_minimal()

p_uds
ggsave("plots/p_uds.png", p_uds, width = 8, height = 5, dpi = 300)

p_ctmm
ggsave("plots/p_ctmm.png", p_ctmm, width = 8, height = 5, dpi = 300)

p_aic
ggsave("plots/p_aic.png", p_aic, width = 8, height = 5, dpi = 300)

p_es
ggsave("plots/p_es.png", p_es, width = 8, height = 5, dpi = 300)
