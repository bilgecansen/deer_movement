# Plots ------------------------------------------------------------------------

library(patchwork)

# Energy score violin plot
p_es <- ggplot(
  results$es,
  aes(x = as.factor(model), y = energy_skill)
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
  (results$es %>%
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
    data = deer_mvt$stp_test[[4]],
    aes(x = x1_, y = y1_),
    color = 'black',
    alpha = 0.5
  ) +
  geom_path(
    data = filter(results$sim[[9]][[4]], nsim == 1),
    aes(x = x_, y = y_),
    color = "orange",
    alpha = 0.5
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
    data = results$sim[[9]][[4]],
    aes(x = x_, y = y_, group = n_sim),
    color = "orange",
    alpha = 0.2
  ) +
  theme_minimal()

ggplot() +
  geom_path(
    data = deer_mvt$stp_test[[9]],
    aes(x = x1_, y = y1_),
    color = 'black',
    alpha = 0.5
  ) +
  geom_path(
    data = filter(results$sim[[9]][[9]], nsim == 1),
    aes(x = x_, y = y_),
    color = "orange",
    alpha = 0.5
  ) +
  theme_minimal()

ggplot() +
  geom_path(
    data = deer_mvt$stp_test[[9]],
    aes(x = x1_, y = y1_),
    color = 'black',
    alpha = 0.5
  ) +
  geom_path(
    data = results$sim[[9]][[9]],
    aes(x = x_, y = y_, group = n_sim),
    color = "orange",
    alpha = 0.2
  ) +
  theme_minimal()
