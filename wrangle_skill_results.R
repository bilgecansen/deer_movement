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

# Deer movement data — needed for test step length distributions (exceedance)
deer_mvt <- readRDS("data_deer_1_119.rds")

# Wrangle UD and CTMM overlap into a long data frame --------------------------
# Each deer may have a different subset of models (the names of results$ud are
# the actual model numbers that were simulated for that deer — null plus any
# models beating null by > 3 delta_logp).
res_bat <- purrr::map2_dfr(results, deer_numbers, function(r, d) {
  model_nums <- as.integer(names(r$ud))
  purrr::map_dfr(model_nums, function(m) {
    ud <- r$ud[[as.character(m)]]
    if (length(ud) == 1 && is.na(ud)) {
      return(data.frame(
        deer = d,
        model = m,
        bat_uds = NA_real_,
        bat_ctmm = NA_real_
      ))
    }
    data.frame(
      deer = d,
      model = m,
      bat_uds = ud$bat_uds,
      bat_ctmm = ud$bat_ctmm
    )
  })
})

# Energy scores (raw, per-model, per-deer) + exceedance probability -----------
# p_excd = P(observed test step length > energy_score), i.e. fraction of the
# deer's actual test steps that are longer than the energy score value.
res_es <- purrr::map2_dfr(results, deer_numbers, function(r, d) {
  sl_test <- deer_mvt$stp_test[[d]]$sl_
  r$es %>%
    dplyr::mutate(
      deer = d,
      p_excd = purrr::map_dbl(energy_score, function(es) {
        if (is.na(es)) {
          return(NA_real_)
        }
        mean(sl_test > es, na.rm = TRUE)
      })
    )
})

# Ensure model column types match before joining
res_bat$model <- as.integer(res_bat$model)
res_es$model <- as.integer(res_es$model)

# Filter: UD overlap (step1), CTMM overlap (step2), ES exceedance (step3) -----
model_selection <- res_bat %>%
  left_join(
    res_es %>% dplyr::select(deer, model, energy_score, p_excd),
    by = c("deer", "model")
  ) %>%
  group_by(deer) %>%
  mutate(
    step1 = !is.na(bat_uds) & bat_uds >= 0.8,
    step2 = step1 & !is.na(bat_ctmm) & bat_ctmm >= 0.8,
    step3 = step2 & !is.na(p_excd) & p_excd >= 0.5
  ) %>%
  ungroup()

null_model <- 2L

selected <- model_selection %>%
  filter(step3) %>%
  select(deer, model, bat_uds, bat_ctmm, energy_score, p_excd) %>%
  group_by(deer) %>%
  # Drop the null model for deer that have at least one alternative selected.
  # Keep it only for deer where the null is the sole passing model.
  filter(!(model == null_model & any(model != null_model))) %>%
  ungroup()

no_selection <- model_selection %>%
  group_by(deer) %>%
  summarize(
    passed_step1 = sum(step1),
    passed_step2 = sum(step2),
    passed_step3 = sum(step3)
  ) %>%
  filter(passed_step3 == 0)

saveRDS(model_selection, "results_sel.rds")

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

# CTMM overlap (BC) violin plot — restricted to models passing step 1
p_ctmm <- ggplot(
  model_selection %>% filter(step1),
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

# ES exceedance violin plot — for models passing step 2
p_excd <- ggplot(
  model_selection %>% filter(step2),
  aes(x = as.factor(model), y = p_excd)
) +
  geom_violin(fill = "#BF5AF2", alpha = 0.5) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  geom_hline(
    yintercept = 0.5,
    col = "red",
    alpha = 0.8,
    linetype = 2,
    linewidth = 2
  ) +
  labs(
    x = "Model",
    y = "P(observed step length > energy score)"
  ) +
  theme_minimal()

# Energy score violin plot — for models passing step 3 (final selection)
p_es <- ggplot(
  model_selection %>% filter(step3),
  aes(x = as.factor(model), y = energy_score)
) +
  geom_violin(fill = "#61D836", alpha = 0.5) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
  labs(
    x = "Model",
    y = "Energy Score"
  ) +
  theme_minimal()

p_uds
ggsave("plots/p_uds.png", p_uds, width = 8, height = 5, dpi = 300)

p_ctmm
ggsave("plots/p_ctmm.png", p_ctmm, width = 8, height = 5, dpi = 300)

p_excd
ggsave("plots/p_excd.png", p_excd, width = 8, height = 5, dpi = 300)

p_es
ggsave("plots/p_es.png", p_es, width = 8, height = 5, dpi = 300)
