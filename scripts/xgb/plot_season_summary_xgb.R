#' @description
#' One figure per season: each variable gets a point at its median and a bar
#' across the spread, with the variable names written once.
#'
#' The spread combines two things that cannot be separated here — how much a
#' variable's importance differs between years, and how much it moves when the
#' procedure is rerun (fresh random points, fresh bagging and column sampling,
#' fresh permutations). A bar clear of zero means the variable mattered in
#' every year and every run; a bar straddling zero means it did not.
#'
#' Values are per 100 observed steps, because deer in different seasons
#' contribute very different numbers of steps each (winter tracks are two to
#' four times longer than autumn ones). Set PER to "deer" for the per-deer
#' figure, which reads against the delta_logp >= 3 gate.
#'
#' Input:  results/xgb/replicates_xgb_<season>_<year>.rds
#'         (falls back to importance_xgb_<season>_<year>.rds, one fit per
#'         year, if the replicates have not been run)
#' Output: plots/importance_xgb_summary_<season>_<per>.png
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
# Variables left out of the plot (still named in the subtitle)
EXCLUDE <- c("HR_center_end")
# Which importance column to summarise: "out_of_bag" or "in_sample"
METRIC <- "out_of_bag"
# "steps100" compares seasons fairly; "deer" reads against the gate
PER <- "steps100"
# Bar range. With fewer than 10 values per variable the full range is drawn
# instead, since percentiles of a handful of numbers say little.
LOWER <- 0.1
UPPER <- 0.9
WIDTH <- 9
HEIGHT <- 5.5

# Load packages ---------------------------------------------------------------
library(tidyverse)
library(scales)

# helper functions (xgb_scale_value / xgb_scale_label)
source("scripts/helper_functions.R")

files <- list.files("results/xgb", pattern = "^replicates_xgb_.*[.]rds$",
                    full.names = TRUE)
source_label <- "years x runs"
if (!length(files)) {
  files <- list.files("results/xgb", pattern = "^importance_xgb_.*[.]rds$",
                      full.names = TRUE)
  source_label <- "years"
}
if (!length(files)) {
  stop("No importance or replicate files in results/xgb/")
}
imp <- purrr::map_dfr(files, readRDS)
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
  season_imp <- imp |>
    filter(season == season_code) |>
    mutate(label_var = dplyr::coalesce(PRETTY[variable], variable))

  summarise_spread <- function(d) {
    d |>
      group_by(label_var) |>
      summarise(
        median = median(value),
        lo = if (dplyr::n() >= 10) quantile(value, LOWER) else min(value),
        hi = if (dplyr::n() >= 10) quantile(value, UPPER) else max(value),
        .groups = "drop"
      )
  }

  dropped <- summarise_spread(season_imp |> filter(variable %in% EXCLUDE)) |>
    mutate(text = sprintf("%s (median %.1f, %.1f to %.1f)", label_var,
                          median, lo, hi))

  df <- summarise_spread(season_imp |> filter(!variable %in% EXCLUDE)) |>
    arrange(median) |>
    mutate(label_var = factor(label_var, levels = label_var))

  # Pale stripe behind every second row, to carry the eye across.
  bands <- tibble(row = seq(2, nrow(df), by = 2))
  # The stripes are a geom, and a theme's gridlines are always drawn beneath
  # geoms, so the lines are drawn here instead, after the stripes, at the
  # same breaks the axis labels use.
  breaks <- scales::breaks_extended(6)(range(c(0, df$lo, df$hi)))
  n_years <- dplyr::n_distinct(season_imp$year)
  n_vals <- nrow(season_imp) / dplyr::n_distinct(season_imp$variable)
  # Deer counts are per season-year, and a deer tracked in two years counts
  # in both, so this is deer-years rather than distinct animals.
  per_year <- season_imp |>
    dplyr::distinct(year, n_deer, n_steps)
  n_deer_years <- sum(per_year$n_deer)
  n_steps_total <- sum(per_year$n_steps)

  p <- ggplot(df, aes(y = label_var)) +
    geom_rect(data = bands, inherit.aes = FALSE,
              aes(ymin = row - 0.5, ymax = row + 0.5),
              xmin = -Inf, xmax = Inf, fill = "#edece5") +
    geom_vline(xintercept = breaks, colour = "#d8d6cd", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "#c3c2b7", linewidth = 0.5) +
    geom_linerange(aes(xmin = lo, xmax = hi), colour = "#eb6834",
                   linewidth = 1.6) +
    geom_point(aes(x = median), colour = "#eb6834", size = 4) +
    scale_x_continuous(breaks = breaks, labels = scales::label_comma(),
                       expand = expansion(mult = c(0.05, 0.05))) +
    labs(
      title = sprintf("Drop in total log score by variable, pooled %s",
                      season_code),
      subtitle = paste0(
        n_deer_years, " deer-years across ", n_years, " years (",
        min(per_year$n_deer), "-", max(per_year$n_deer), " per year), ",
        formatC(n_steps_total, format = "d", big.mark = ","),
        " observed steps.\n",
        "Out-of-bag drop when a variable is shuffled within strata (FAMD: ",
        "among forest points).\nPoint: median over ", n_years, " years",
        if (source_label == "years x runs") {
          sprintf(" x %d runs", round(n_vals / n_years))
        } else {
          ""
        },
        ". Bar: ",
        if (n_vals >= 10) {
          sprintf("%gth-%gth percentile", 100 * LOWER, 100 * UPPER)
        } else {
          "full range"
        },
        ".\nOff the scale and left out: ", dropped$text, "."
      ),
      x = xgb_scale_label(PER),
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.background = element_rect(fill = "#fcfcfb", colour = NA),
      panel.background = element_rect(fill = "#fcfcfb", colour = NA),
      panel.grid = element_blank(),
      axis.text.y = element_text(colour = "#0b0b0b", size = 10),
      axis.text.x = element_text(colour = "#898781", size = 9),
      axis.title.x = element_text(colour = "#52514e", size = 10),
      plot.title = element_text(colour = "#0b0b0b", face = "bold"),
      plot.subtitle = element_text(colour = "#52514e", size = 9)
    )

  out <- sprintf("plots/importance_xgb_summary_%s_%s.png", season_code, PER)
  ggsave(out, p, width = WIDTH, height = HEIGHT, dpi = 150, bg = "#fcfcfb")
  cat(sprintf("-> %s\n", out))
}

for (s in sort(unique(imp$season))) {
  plot_season(s)
}
