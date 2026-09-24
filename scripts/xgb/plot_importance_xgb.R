#' @description
#' Plot the variable importance that fit_season_xgb.R has scored: one figure
#' per season, one row per variable.
#'
#' The noise column is drawn as a row of its own and as a dashed reference
#' line, because a variable earns its place by beating the noise, not by
#' clearing zero. On the pooled fits it lands within a few hundredths of
#' zero in every season, so the two nearly coincide — that agreement is the
#' point, and it is worth being able to see it rather than assume it.
#'
#' Nothing is left off the scale: distance to the home-range centre is in
#' the nuisance block, so it is not ranked and cannot dominate the axis.
#'
#' Values are per deer-year by default, the scale that reads against the
#' pipeline's delta_logp >= 3 gate. Set PER to "steps100" for the per-100-
#' step version. Neither makes seasons strictly comparable — a season with
#' fewer steps over-fits more and leans harder on every variable, which
#' inflates its whole column. Compare ranks between seasons, not values.
#'
#' Input:  results/xgb/season_xgb_<season>.rds
#' Output: plots/importance_xgb_<season>_<per>.png
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
# "deer" reads against the gate; "steps100" puts seasons on one axis
PER <- "deer"
WIDTH <- 9
HEIGHT <- 5.5

# Load packages ---------------------------------------------------------------
library(tidyverse)
library(scales)

# helper functions (xgb_scale_value / xgb_scale_label)
source("scripts/helper_functions.R")

files <- list.files("results/xgb", pattern = "^season_xgb_.*[.]rds$",
                    full.names = TRUE)
if (!length(files)) {
  stop("No season files in results/xgb/; run fit_season_xgb.R first")
}

PRETTY <- c(
  ndvi_end = "NDVI",
  landcover = "Landcover",
  forest_edge_end = "Forest edge distance",
  elevation_end = "Elevation",
  northness_end = "Northness",
  eastness_end = "Eastness",
  oak_mast_end = "Oak mast",
  oak_dist_end = "Oak distance",
  shadow_gauss = "Noise column",
  famd1_end = "FAMD1",
  famd2_end = "FAMD2",
  famd3_end = "FAMD3",
  famd4_end = "FAMD4",
  famd5_end = "FAMD5"
)

dir.create("plots", showWarnings = FALSE)

plot_season <- function(path) {
  r <- readRDS(path)
  df <- r$importance |>
    mutate(value = xgb_scale_value(cv, n_deer, n_steps, PER),
           label = dplyr::coalesce(PRETTY[variable], variable),
           is_noise = variable == "shadow_gauss") |>
    arrange(value) |>
    mutate(label = factor(label, levels = label))

  noise <- df$value[df$is_noise]
  # Pale stripe behind every second row, to carry the eye across.
  bands <- tibble(row = seq(2, nrow(df), by = 2))
  # The stripes are a geom, and a theme's gridlines are always drawn beneath
  # geoms, so the lines are drawn here instead, after the stripes, at the
  # same breaks the axis labels use.
  breaks <- scales::breaks_extended(6)(range(c(0, df$value)))

  p <- ggplot(df, aes(y = label, x = value)) +
    geom_rect(data = bands, inherit.aes = FALSE,
              aes(ymin = row - 0.5, ymax = row + 0.5),
              xmin = -Inf, xmax = Inf, fill = "#edece5") +
    geom_vline(xintercept = breaks, colour = "#d8d6cd", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "#c3c2b7", linewidth = 0.5) +
    geom_vline(xintercept = noise, colour = "#2a78d6", linewidth = 0.5,
               linetype = "22") +
    geom_segment(aes(x = 0, xend = value, yend = label, colour = is_noise),
                 linewidth = 0.9) +
    geom_point(aes(colour = is_noise), size = 4) +
    scale_colour_manual(values = c(`FALSE` = "#eb6834",
                                   `TRUE` = "#2a78d6"),
                        guide = "none") +
    scale_x_continuous(breaks = breaks, labels = scales::label_comma(),
                       expand = expansion(mult = c(0.05, 0.05))) +
    labs(
      title = sprintf("Variable importance, pooled %s, all years",
                      r$season),
      subtitle = paste0(
        r$n_deer_years, " deer-years from ", r$n_animals, " animals over ",
        nrow(r$per_year), " years (",
        formatC(r$n_steps, format = "d", big.mark = ","),
        " observed steps).\nDrop in held-out log score when a variable is ",
        "shuffled within strata (FAMD: among forest points), over ",
        "five random\nfolds of whole steps, above a movement + HR-centre ",
        "null. Trees may combine variables; depth 2; 500 rounds.\n",
        "The dashed line is a column of pure noise carried through the ",
        "whole fit."
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
      plot.subtitle = element_text(colour = "#52514e", size = 8.5)
    )

  out <- sprintf("plots/importance_xgb_%s_%s.png", r$season, PER)
  ggsave(out, p, width = WIDTH, height = HEIGHT, dpi = 150, bg = "#fcfcfb")
  cat(sprintf("-> %s\n", out))
}

for (f in files) {
  plot_season(f)
}
