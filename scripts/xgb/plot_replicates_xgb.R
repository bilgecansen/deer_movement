#' @description
#' Plot the replicate importance from replicate_importance_xgb.R: a point at
#' each variable's median and a bar across the middle of its runs.
#'
#' A variable whose bar clears zero carries signal; one straddling zero
#' scores positive or negative depending on the draw, which is the same
#' practical verdict as negligible.
#'
#' The x axis is linear. Distance to the home-range centre is an order of
#' magnitude above the rest and would flatten everything, so EXCLUDE drops it
#' from the plot and the subtitle reports its numbers instead.
#'
#' Input:  results/xgb/replicates_xgb_<season>_<year>.rds
#' Output: plots/importance_replicates_xgb_<season>_<year>.png
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
SEASON <- "fa"
YEAR <- 2021L
# Variables left out of the plot (still named in the subtitle)
EXCLUDE <- c("HR_center_end")
# Which importance column to summarise: "out_of_bag" or "in_sample"
METRIC <- "out_of_bag"
# Scale: "deer" is the headline (reads against the delta_logp >= 3 gate),
# "steps100" for comparing season-years, "total" for the raw drop
PER <- "deer"
# Interval drawn across the replicates
LOWER <- 0.1
UPPER <- 0.9
WIDTH <- 9.5
HEIGHT <- 5.5

# Load packages ---------------------------------------------------------------
library(tidyverse)
library(scales)

# helper functions (xgb_scale_value / xgb_scale_label)
source("scripts/helper_functions.R")

key <- sprintf("%s_%d", SEASON, YEAR)
reps <- readRDS(sprintf("results/xgb/replicates_xgb_%s.rds", key))
if (!METRIC %in% names(reps)) {
  stop(sprintf("No '%s' column in the replicate file", METRIC))
}
if (!all(c("n_deer", "n_steps") %in% names(reps))) {
  stop("Replicate file lacks n_deer / n_steps; rerun the replicates")
}
reps$value <- xgb_scale_value(reps[[METRIC]], reps$n_deer, reps$n_steps,
                              PER)
out_path <- sprintf("plots/importance_replicates_xgb_%s.png", key)
dir.create("plots", showWarnings = FALSE)

# Display names; anything not listed keeps its column name.
PRETTY <- c(
  HR_center_end = "HR centre distance",
  ndvi_end = "NDVI",
  landcover = "Landcover",
  oak_mast_end = "Oak mast",
  oak_dist_end = "Oak distance",
  forest_edge_end = "Forest edge distance",
  elevation_end = "Elevation",
  northness_end = "Northness",
  eastness_end = "Eastness",
  famd1_end = "FAMD1",
  famd2_end = "FAMD2",
  famd3_end = "FAMD3",
  famd4_end = "FAMD4",
  famd5_end = "FAMD5"
)

summ <- reps |>
  group_by(variable) |>
  summarise(
    median = median(value),
    lower = quantile(value, LOWER),
    upper = quantile(value, UPPER),
    share_positive = mean(value > 0),
    .groups = "drop"
  ) |>
  mutate(label_var = dplyr::coalesce(PRETTY[variable], variable))

n_reps <- dplyr::n_distinct(reps$rep)
dropped <- summ |>
  filter(variable %in% EXCLUDE) |>
  mutate(text = sprintf("%s (median %.1f, %.1f to %.1f)", label_var,
                        median, lower, upper))

df <- summ |>
  filter(!variable %in% EXCLUDE) |>
  mutate(label_var = fct_reorder(label_var, median))

p <- ggplot(df, aes(y = label_var)) +
  geom_vline(xintercept = 0, colour = "#c3c2b7", linewidth = 0.5) +
  geom_linerange(aes(xmin = lower, xmax = upper), colour = "#eb6834",
                 linewidth = 1.8) +
  geom_point(aes(x = median), colour = "#eb6834", size = 4.5) +
  scale_x_continuous(
    labels = scales::label_comma(),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  labs(
    title = sprintf("Variable importance over %d runs, pooled %s %d",
                    n_reps, SEASON, YEAR),
    subtitle = paste0(
      "Out-of-bag drop in conditional log-likelihood. Point: median. Bar: ",
      100 * LOWER, "th-", 100 * UPPER, "th percentile across runs.\n",
      "Each run redraws the random points and refits with a new seed; the ",
      "deer and their steps are fixed.",
      if (nrow(dropped)) {
        sprintf("\nOff the scale and left out: %s.",
                paste(dropped$text, collapse = ", "))
      } else {
        ""
      }
    ),
    x = xgb_scale_label(PER),
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "#fcfcfb", colour = NA),
    panel.background = element_rect(fill = "#fcfcfb", colour = NA),
    panel.grid.major.x = element_line(colour = "#e1e0d9", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(colour = "#0b0b0b", size = 11),
    axis.text.x = element_text(colour = "#898781", size = 9),
    axis.title.x = element_text(colour = "#52514e", size = 10),
    plot.title = element_text(colour = "#0b0b0b", face = "bold"),
    plot.subtitle = element_text(colour = "#52514e", size = 9)
  )

ggsave(out_path, p, width = WIDTH, height = HEIGHT, dpi = 150,
       bg = "#fcfcfb")
cat(sprintf("-> %s\n", out_path))
