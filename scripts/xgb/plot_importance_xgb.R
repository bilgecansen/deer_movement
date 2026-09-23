#' @description
#' Plot the variable importance from fit_xgb.R: one horizontal line per
#' variable ending in a point, ordered by importance.
#'
#' The x axis is linear, so the bars are read as they are. Distance to the
#' home-range centre is an order of magnitude above every habitat variable
#' and would flatten the rest against zero, so EXCLUDE drops it from the
#' plot; the subtitle still reports what it scored.
#'
#' Input:  results/xgb/importance_xgb_<season>_<year>.rds
#' Output: plots/importance_xgb_<season>_<year>.png
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
SEASON <- "fa"
YEAR <- 2021L
# Variables left out of the plot (still named in the subtitle)
EXCLUDE <- c("HR_center_end")
# Which importance column to plot: "out_of_bag" or "in_sample"
METRIC <- "out_of_bag"
# Scale: "deer" is the headline (reads against the delta_logp >= 3 gate),
# "steps100" for comparing season-years, "total" for the raw drop
PER <- "deer"
# Plot size in inches
WIDTH <- 9
HEIGHT <- 5

# Load packages ---------------------------------------------------------------
library(tidyverse)
library(scales)

# helper functions (xgb_scale_value / xgb_scale_label)
source("scripts/helper_functions.R")

key <- sprintf("%s_%d", SEASON, YEAR)
imp <- readRDS(sprintf("results/xgb/importance_xgb_%s.rds", key))
if (!METRIC %in% names(imp)) {
  stop(sprintf("No '%s' column; refit with BAG_FRAC < 1 for out-of-bag",
               METRIC))
}
if (!all(c("n_deer", "n_steps") %in% names(imp))) {
  stop("Importance file lacks n_deer / n_steps; refit with fit_xgb.R")
}
imp$importance <- xgb_scale_value(imp[[METRIC]], imp$n_deer, imp$n_steps,
                                  PER)
# Excluding a variable changes the axis, so the two versions get their own
# files rather than overwriting each other.
out_path <- sprintf("plots/importance_xgb_%s%s.png", key,
                    if (length(EXCLUDE)) "_trimmed" else "")
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

pretty_name <- function(v) dplyr::coalesce(PRETTY[v], v)

dropped <- imp |>
  filter(variable %in% EXCLUDE) |>
  mutate(text = sprintf("%s (%s)", pretty_name(variable),
                        formatC(importance, format = "f", digits = 0,
                                big.mark = ",")))

df <- imp |>
  filter(!variable %in% EXCLUDE) |>
  mutate(
    label_var = pretty_name(variable),
    label_var = fct_reorder(label_var, importance),
    value = formatC(importance, format = "f", digits = 2, big.mark = ",")
  )

# Value labels sit a fixed distance from each point, in plot units.
GAP <- 0.02 * diff(range(c(0, df$importance)))
df <- df |>
  mutate(
    x_label = importance + ifelse(importance >= 0, GAP, -GAP),
    label_hjust = ifelse(importance >= 0, 0, 1)
  )

p <- ggplot(df, aes(x = importance, y = label_var)) +
  geom_vline(xintercept = 0, colour = "#c3c2b7", linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = importance, yend = label_var),
               colour = "#eb6834", linewidth = 1.1, lineend = "round") +
  geom_point(colour = "#eb6834", size = 4.5) +
  geom_text(aes(x = x_label, label = value, hjust = label_hjust),
            colour = "#0b0b0b", size = 3.6) +
  scale_x_continuous(
    labels = scales::label_comma(),
    # Room for the value labels: they sit left of a negative point and right
    # of a positive one.
    expand = expansion(
      mult = c(if (any(df$importance < 0)) 0.14 else 0.02, 0.12)
    )
  ) +
  labs(
    title = sprintf("Variable importance, pooled %s %d", SEASON, YEAR),
    subtitle = paste0(
      "Drop in conditional log-likelihood when the variable is shuffled ",
      "within strata (FAMD: among forest points).\n",
      if (METRIC == "out_of_bag") {
        paste("Out-of-bag: each step scored by the trees that did not",
              "learn from it.")
      } else {
        "In-sample, so this ranks variables within one fit."
      },
      if (nrow(dropped)) {
        sprintf(" Off the scale and left out: %s.",
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
    plot.subtitle = element_text(colour = "#52514e", size = 9.5)
  )

ggsave(out_path, p, width = WIDTH, height = HEIGHT, dpi = 150,
       bg = "#fcfcfb")
cat(sprintf("-> %s\n", out_path))
