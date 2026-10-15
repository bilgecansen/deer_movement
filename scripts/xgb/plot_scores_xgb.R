#' @description
#' Violin plots of the path-simulation scores: one panel per score, one
#' violin per season, one point per deer-year. A deer-year's value is that
#' of whichever of its season's deer models scored better on that score
#' (higher is better for all three), so each score can come from a
#' different model.
#'
#' The energy score enters as p_excd, as in the amt and GAM filters: the
#' share of the deer-year's real step lengths longer than the energy score,
#' rather than in metres, which would mostly reflect how far each deer
#' moves.
#'
#' Input:  filters/xgb/scores_<set>/scores_xgb_<key>.rds (score_sims_xgb.R);
#'         read per deer-year, so a run still in progress can be plotted
#' Output: plots/scores_xgb_<set>.png and .pdf, and a summary printed
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
# The deer-years to plot: "test" or "train"
SET <- "test"
# Dashed reference line in each panel
REFERENCE <- c(hr_overlap = 0.8, svf_agree = 0.8, p_excd = 0.45)
WIDTH <- 12
HEIGHT <- 5.5

# Load packages ---------------------------------------------------------------
library(tidyverse)

files <- list.files(sprintf("filters/xgb/scores_%s", SET),
                    pattern = "^scores_xgb_.*[.]rds$", full.names = TRUE)
if (!length(files)) {
  stop(sprintf("No scores for set '%s'; run score_sims_xgb.R first", SET))
}
scores <- purrr::map_dfr(files, readRDS) |> filter(deer)

SEASON <- c(br = "breeding", fa = "fawning", nb = "non-breeding",
            pf = "post-fawning")
METRIC <- c(hr_overlap = "HR overlap (1 is best)",
            svf_agree = "svf agreement (1 is best)",
            p_excd = "P(real step > energy score) (1 is best)")
LOWER_BETTER <- character(0)

# The better of the deer models, per deer-year and score
long <- scores |>
  pivot_longer(all_of(names(METRIC)), names_to = "metric",
               values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(signed = if_else(metric %in% LOWER_BETTER, -value, value))
better <- long |>
  group_by(key, season, metric) |>
  slice_max(signed, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(winner = if_else(simplest, "simplest", "best"))

cat(sprintf("=== %s: %d deer-years scored ===\n", SET,
            dplyr::n_distinct(scores$key)))
cat(sprintf("model runs with a failed ctmm fit: %d of %d\n",
            sum(!is.na(scores$ud_error)), nrow(scores)))
print(as.data.frame(
  better |>
    group_by(metric, season) |>
    summarise(n = n(), median = median(value),
              q25 = quantile(value, 0.25), q75 = quantile(value, 0.75),
              at_reference = mean(value >= REFERENCE[metric[1]]),
              simplest_wins = mean(winner == "simplest"), .groups = "drop") |>
    arrange(match(metric, names(METRIC)), season)
), row.names = FALSE, digits = 3)

n_per <- better |> distinct(key, season) |> count(season)
xlab <- stats::setNames(sprintf("%s\nn = %d", SEASON[n_per$season], n_per$n),
                        n_per$season)
plot_data <- better |>
  mutate(metric = factor(METRIC[metric], levels = METRIC),
         season = factor(season, levels = names(SEASON)))

reference <- tibble(metric = factor(METRIC[names(REFERENCE)], levels = METRIC),
                    value = unname(REFERENCE))

p <- ggplot(plot_data, aes(season, value)) +
  geom_hline(data = reference, aes(yintercept = value), linetype = "22",
             linewidth = 0.5, colour = "#52514e") +
  geom_violin(fill = "#2a78d6", alpha = 0.18, colour = "#2a78d6",
              linewidth = 0.5, scale = "width", trim = TRUE) +
  geom_jitter(width = 0.12, height = 0, size = 0.9, alpha = 0.55,
              colour = "#52514e") +
  stat_summary(fun = median, geom = "crossbar", width = 0.45,
               linewidth = 0.4, colour = "#0b0b0b") +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_x_discrete(labels = xlab) +
  labs(
    title = sprintf(
      "Simulated against real paths, %s deer-years, the better deer model",
      if (SET == "test") "2021 test" else "training"
    ),
    subtitle = paste0(
      "One point per deer-year: the value of whichever of its season's two ",
      "models with deer (simplest, best) scored better on that score. ",
      "Black bar: median. Dashed line: reference value.\n10 simulated ",
      "paths per model, compared with the ",
      "real path at the real fix times. Energy score: the share of the ",
      "deer-year's real steps longer than it (p_excd, as for amt and GAM)."
    ),
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "#fcfcfb", colour = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "#e1e0d9", linewidth = 0.3),
    strip.text = element_text(face = "bold", colour = "#0b0b0b",
                              hjust = 0),
    axis.text = element_text(colour = "#52514e"),
    plot.title = element_text(face = "bold", colour = "#0b0b0b"),
    plot.subtitle = element_text(colour = "#52514e", size = 9)
  )

dir.create("plots", showWarnings = FALSE)
out <- sprintf("plots/scores_xgb_%s", SET)
ggsave(paste0(out, ".png"), p, width = WIDTH, height = HEIGHT, dpi = 150,
       bg = "#fcfcfb")
ggsave(paste0(out, ".pdf"), p, width = WIDTH, height = HEIGHT,
       bg = "#fcfcfb")
cat(sprintf("-> %s.png, %s.pdf\n", out, out))
