# Step 7: Select best models for each deer
cat("Selecting best models...\n")

# Wrangle res_ud into a data frame
res_bat <- foreach(m = names(results$ud), .combine = "rbind") %do%
  {
    scores_uds <- map_dbl(i = 1:n_deer, function(i) {
      if (length(results$ud[[m]][[i]]) == 1 && is.na(results$ud[[m]][[i]])) {
        return(NA_real_)
      }
      results$ud[[m]][[i]]$bat_uds
    })

    scores_ctmm <- map_dbl(i = 1:n_deer, function(i) {
      if (length(results$ud[[m]][[i]]) == 1 && is.na(results$ud[[m]][[i]])) {
        return(NA_real_)
      }
      results$ud[[m]][[i]]$bat_ctmm
    })

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
