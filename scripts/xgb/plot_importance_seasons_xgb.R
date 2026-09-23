#' @description
#' Plot the variable importance that run_all_xgb.R has fitted: one figure per
#' season, one panel per year within it.
#'
#' Each panel has its own x axis: a season-year with more steps gives larger
#' drops, so the panels are read within themselves, not against each other.
#' The panel strip carries the deer and step counts.
#'
#' Each season shows only the variables its own models used, so NDVI does not
#' appear in the winter figure.
#'
#' Input:  results/xgb/importance_xgb_<season>_<year>.rds
#' Output: plots/importance_xgb_season_<season>.png
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
# Variables left out of the plot (still named in the subtitle)
EXCLUDE <- c("HR_center_end")
# Which importance column to plot: "out_of_bag" or "in_sample"
METRIC <- "out_of_bag"
# Scale: "steps100" for cross-season figures (deer differ in track length
# between seasons), "deer" for the headline per-deer number, "total" for raw
PER <- "steps100"
# Panel columns, and the size of one panel in inches
N_COLS <- 3
PANEL_WIDTH <- 4.2
PANEL_HEIGHT <- 3.6

# Load packages ---------------------------------------------------------------
library(tidyverse)
library(scales)

# helper functions (xgb_scale_value / xgb_scale_label)
source("scripts/helper_functions.R")

files <- list.files("results/xgb", pattern = "^importance_xgb_.*[.]rds$",
                    full.names = TRUE)
if (!length(files)) {
  stop("No importance files in results/xgb/; run run_all_xgb.R first")
}
imp <- purrr::map_dfr(files, readRDS)
if (!all(c("season", "year") %in% names(imp))) {
  stop("Importance files lack season / year; refit with run_all_xgb.R")
}
imp$value <- xgb_scale_value(imp[[METRIC]], imp$n_deer, imp$n_steps, PER)

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

dir.create("plots", showWarnings = FALSE)

plot_season <- function(season_code) {
  season_imp <- imp |> filter(season == season_code)
  dropped <- season_imp |>
    filter(variable %in% EXCLUDE) |>
    summarise(text = sprintf("%s (%.1f to %.1f across years)",
                             dplyr::coalesce(PRETTY[variable[1]],
                                             variable[1]),
                             min(value), max(value)))

  df <- season_imp |>
    filter(!variable %in% EXCLUDE) |>
    mutate(
      label_var = dplyr::coalesce(PRETTY[variable], variable),
      panel = sprintf("%s %d  ·  %d deer, %s steps", season, year, n_deer,
                      formatC(n_steps, format = "d", big.mark = ","))
    )
  # One variable order for the season's panels, by the median across its
  # years. Only variables this season's models used appear.
  order_tab <- df |>
    group_by(label_var) |>
    summarise(m = median(value), .groups = "drop") |>
    arrange(m)
  df <- df |>
    mutate(
      label_var = factor(label_var, levels = order_tab$label_var),
      panel = factor(panel, levels = unique(panel[order(year)]))
    )

  n_panels <- dplyr::n_distinct(df$panel)
  n_cols <- min(N_COLS, n_panels)
  n_rows <- ceiling(n_panels / n_cols)

  p <- ggplot(df, aes(x = value, y = label_var)) +
    geom_vline(xintercept = 0, colour = "#c3c2b7", linewidth = 0.4) +
    geom_point(colour = "#eb6834", size = 3) +
    # One x axis for the season: per-100-step values are comparable between
    # years, which is the reason for scaling them.
    facet_wrap(~panel, ncol = n_cols,
               scales = if (PER == "total") "free_x" else "fixed") +
    scale_x_continuous(labels = scales::label_comma(),
                       expand = expansion(mult = c(0.1, 0.1))) +
    labs(
      title = sprintf("Drop in total log score by variable, pooled %s",
                      season_code),
      subtitle = paste0(
        "Out-of-bag drop in the total log score when a variable is ",
        "shuffled within strata (FAMD: among forest points).\n",
        switch(
          PER,
          steps100 = paste("Scored per 100 observed steps, so years with",
                           "longer tracks are not inflated."),
          deer = paste("Scored per deer, the same units as delta_logp.",
                       "Deer with longer tracks contribute more, so winter",
                       "years run higher."),
          total = "Raw totals over the whole pool."
        ),
        "\nOff the scale and left out: ", dropped$text, "."
      ),
      x = xgb_scale_label(PER),
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.background = element_rect(fill = "#fcfcfb", colour = NA),
      panel.background = element_rect(fill = "#fcfcfb", colour = NA),
      panel.grid.major.x = element_line(colour = "#e1e0d9", linewidth = 0.3),
      panel.grid.major.y = element_line(colour = "#e1e0d9", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(colour = "#0b0b0b", size = 9),
      axis.text.x = element_text(colour = "#898781", size = 8),
      axis.title.x = element_text(colour = "#52514e", size = 10),
      strip.text = element_text(colour = "#0b0b0b", hjust = 0, size = 9.5),
      plot.title = element_text(colour = "#0b0b0b", face = "bold"),
      plot.subtitle = element_text(colour = "#52514e", size = 9),
      panel.spacing.x = grid::unit(14, "pt"),
      panel.spacing.y = grid::unit(12, "pt")
    )

  # The scale is in the filename, so a per-deer and a per-100-step version
  # can sit side by side.
  out <- sprintf("plots/importance_xgb_season_%s_%s.png", season_code, PER)
  ggsave(out, p, width = n_cols * PANEL_WIDTH,
         height = n_rows * PANEL_HEIGHT + 1.4, dpi = 150, bg = "#fcfcfb")
  cat(sprintf("-> %s\n", out))
}

for (s in sort(unique(imp$season))) {
  plot_season(s)
}
