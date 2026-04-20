library(tidyverse)
source("helper_functions.R")

model_selection <- readRDS("results_sel.rds")

selected <- model_selection %>%
  filter(step4) %>%
  group_by(deer) %>%
  filter(energy_skill == max(energy_skill)) %>%
  ungroup() %>%
  select(deer, model, bat_uds, bat_ctmm, delta_aic, energy_skill)

no_selection <- model_selection %>%
  group_by(deer) %>%
  summarize(
    passed_step1 = sum(step1, na.rm = TRUE),
    passed_step2 = sum(step2, na.rm = TRUE),
    passed_step3 = sum(step3, na.rm = TRUE),
    passed_step4 = sum(step4, na.rm = TRUE)
  ) %>%
  filter(passed_step4 == 0)

# Selected deers
selected[which.max(selected$bat_uds), ]
plot_deer_paths(48)
ggsave("plots/deer_paths_48.png", width = 10, height = 8, dpi = 300)

selected[which.min(selected$bat_uds), ]
plot_deer_paths(71)
ggsave("plots/deer_paths_71.png", width = 10, height = 8, dpi = 300)

selected[which.max(selected$energy_skill), ]
plot_deer_paths(10)
ggsave("plots/deer_paths_10.png", width = 10, height = 8, dpi = 300)

# Filtered out deers
no_selection
plot_deer_paths(7)
ggsave("plots/deer_paths_7.png", width = 10, height = 8, dpi = 300)

plot_deer_paths(85)
ggsave("plots/deer_paths_85.png", width = 10, height = 8, dpi = 300)

# All deers in pdf
# Best model per selected deer (highest energy_skill among passing models)
best_selected <- selected %>%
  group_by(deer) %>%
  slice_max(energy_skill, n = 1, with_ties = FALSE) %>%
  ungroup()

pdf("plots/deer_paths_selected.pdf", width = 10, height = 8)
for (r in sort(best_selected$deer)) {
  info <- best_selected %>% filter(deer == r)
  subtitle <- sprintf(
    "Selected model: %s | BA-UDS: %.3f | BA-CTMM: %.3f | dAIC: %.2f | energy skill: %.3f",
    info$model,
    info$bat_uds,
    info$bat_ctmm,
    info$delta_aic,
    info$energy_skill
  )
  p <- plot_deer_paths(r) + labs(subtitle = subtitle)
  print(p)
}
dev.off()

pdf(
  "plots/deer_paths_no_selection.pdf",
  width = 10,
  height = 8,
  onefile = TRUE
)
for (r in sort(unique(no_selection$deer))) {
  print(plot_deer_paths(r))
}
dev.off()

# Animations
deer_mvt <- readRDS("data_deer_1_119.rds")
row_no <- 11
obs <- deer_mvt$stp_test[[row_no]]

results_sim <- readRDS(sprintf("results/results_sim_%d.rds", row_no))
sim_one <- results_sim[[8]] %>%
  dplyr::filter(nsim == 1) %>%
  dplyr::rename(x1_ = x_, y1_ = y_, t1_ = t_)

sim_two <- results_sim[[2]] %>%
  dplyr::filter(nsim == 1) %>%
  dplyr::rename(x1_ = x_, y1_ = y_, t1_ = t_)


# Compare observed vs one simulated replicate from a chosen model
animate_deer_path(
  path = obs,
  path2 = sim_one,
  path3 = sim_two,
  step_duration = 0.7,
  wake_length = 0.01,
  file = "plots/deer_11-(2,8)_obs_vs_sim.mp4"
)

# Single Deer
animate_deer_path(obs, file = "plots/deer_48.mp4", step_duration = 0.5)
