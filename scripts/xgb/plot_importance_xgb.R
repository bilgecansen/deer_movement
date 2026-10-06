#' @description
#' Plot the variable importance of each season's SELECTED model type: one
#' figure per season, one row per variable.
#'
#' Which type that is comes from the same rule plot_models_xgb.R applies —
#' the simplest one within THRESHOLD log units per deer of the best — so
#' the two scripts always describe the same model. Keep the two
#' configuration blocks in step.
#'
#' Each ranked block carries a noise column, and they are drawn in blue.
#' A variable earns its place by beating its block's noise column, not by
#' clearing zero: the habitat one has come out slightly negative while the
#' FAMD one sits clearly above zero, so reading the FAMD axes against zero
#' overstates them.
#'
#' Two nulls appear on one axis when the selected model carries the
#' start-of-step block. End-point columns are shuffled within a stratum —
#' does this predict which of these points was chosen. Start columns are
#' shuffled between strata a whole step at a time — does the movement
#' kernel depend on where the animal stood. Same units, same fit, different
#' questions: compare ranks, not an equivalence. The within-stratum score
#' of a start column is exactly zero by construction and is left off.
#'
#' Nothing is off the scale: distance to the home-range centre is in the
#' nuisance block, so it is not ranked and cannot dominate the axis — not
#' even in rsf_hr, where the habitat trees may also use it.
#'
#' Values are per deer-year by default, the scale that reads against the
#' pipeline's delta_logp >= 3 gate; on that scale a dashed line marks
#' THRESHOLD, the bar a variable has to clear to be carried into the
#' individual-deer models. Neither scale makes seasons strictly
#' comparable — a season with fewer steps over-fits more and leans harder
#' on every variable, which inflates its whole column. Compare ranks
#' between seasons, not values.
#'
#' Input:  results/xgb/compare_<season>_<type>.rds
#' Output: plots/importance_xgb_<season>_<type>.png
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
# Must match plot_models_xgb.R, or this figure describes a model the
# selection did not pick.
THRESHOLD <- 3
COMPLEXITY <- c("main", "rsf", "rsf_hr", "start", "full")
# "deer" reads against the gate; "steps100" puts seasons on one axis
PER <- "deer"
WIDTH <- 9
HEIGHT <- 5.5

# Load packages ---------------------------------------------------------------
library(tidyverse)
library(scales)

# helper functions (xgb_scale_value / xgb_scale_label)
source("scripts/helper_functions.R")

files <- list.files("results/xgb", pattern = "^compare_.*[.]rds$",
                    full.names = TRUE)
if (!length(files)) {
  stop("No model files in results/xgb/; run run_models_xgb.R first")
}
scores <- purrr::map_dfr(files, function(f) {
  x <- readRDS(f)
  tibble(season = x$season, config = x$config, held_out = x$ll_cv,
         n_deer = x$n_deer_years)
})
selected <- scores |>
  group_by(season) |>
  mutate(per_deer = (held_out - max(held_out)) / n_deer) |>
  filter(per_deer > -THRESHOLD) |>
  slice_min(match(config, COMPLEXITY), n = 1) |>
  ungroup()

PRETTY <- c(
  ndvi_end = "NDVI",
  landcover = "Landcover",
  forest_edge_end = "Forest edge distance",
  elevation_end = "Elevation",
  northness_end = "Northness",
  eastness_end = "Eastness",
  oak_mast_end = "Oak mast",
  oak_dist_end = "Oak distance",
  shadow_gauss = "Noise (habitat)",
  shadow_famd = "Noise (FAMD)",
  shadow_start = "Noise (start, per step)",
  famd1_end = "FAMD1",
  famd2_end = "FAMD2",
  famd3_end = "FAMD3",
  famd4_end = "FAMD4",
  famd5_end = "FAMD5",
  ndvi_start = "NDVI at start",
  wiscland_start = "Cover at start",
  forest_edge_start = "Forest edge at start"
)

dir.create("plots", showWarnings = FALSE)

plot_season <- function(season, config) {
  r <- readRDS(sprintf("results/xgb/compare_%s_%s.rds", season, config))
  # A start column scored the within-stratum way is exactly zero by
  # construction, so it is dropped rather than drawn as a real zero.
  df <- r$importance |>
    filter(!(booster == "modifier" & scheme == "within")) |>
    mutate(value = xgb_scale_value(cv, r$n_deer_years, r$n_steps, PER),
           label = dplyr::coalesce(PRETTY[variable], variable),
           kind = dplyr::case_when(
             grepl("^shadow", variable) ~ "noise",
             scheme == "across" ~ "start",
             TRUE ~ "end"
           )) |>
    arrange(value) |>
    mutate(label = factor(label, levels = label))

  # Pale stripe behind every second row, to carry the eye across.
  bands <- tibble(row = seq(2, nrow(df), by = 2))
  # The stripes are a geom, and a theme's gridlines are always drawn
  # beneath geoms, so the lines are drawn here instead, after the stripes,
  # at the same breaks the axis labels use.
  # The threshold is per deer, so it is drawn only on that scale.
  show_line <- PER == "deer"
  line_at <- if (show_line) THRESHOLD else 0
  threshold_line <- if (show_line) {
    geom_vline(xintercept = THRESHOLD, colour = "#2a78d6", linewidth = 0.5,
               linetype = "22")
  }
  breaks <- scales::breaks_extended(6)(range(c(0, line_at, df$value)))
  has_start <- any(df$kind == "start")

  p <- ggplot(df, aes(y = label, x = value, colour = kind)) +
    geom_rect(data = bands, inherit.aes = FALSE,
              aes(ymin = row - 0.5, ymax = row + 0.5),
              xmin = -Inf, xmax = Inf, fill = "#edece5") +
    geom_vline(xintercept = breaks, colour = "#d8d6cd", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "#c3c2b7", linewidth = 0.5) +
    threshold_line +
    geom_segment(aes(x = 0, xend = value, yend = label), linewidth = 0.9) +
    geom_point(size = 4) +
    scale_colour_manual(
      values = c(end = "#eb6834", start = "#2a78d6", noise = "#898781"),
      breaks = c("end", "start", "noise"),
      labels = c(end = "end point (selection)",
                 start = "start (movement kernel)",
                 noise = "noise column"),
      name = NULL
    ) +
    scale_x_continuous(breaks = breaks, labels = scales::label_comma(),
                       expand = expansion(mult = c(0.05, 0.05))) +
    expand_limits(x = line_at) +
    labs(
      title = sprintf("Variable importance, pooled %s (%s model)",
                      r$season, r$config),
      subtitle = paste0(
        r$n_deer_years, " deer-years",
        # Fits written before n_animals was recorded simply omit it.
        if (!is.null(r$n_animals)) {
          paste0(" from ", r$n_animals, " animals")
        } else {
          ""
        },
        " (", formatC(r$n_steps, format = "d", big.mark = ","),
        " observed steps). Drop in held-out log score when a variable ",
        "is\nshuffled, over five random folds of whole steps, above a ",
        "movement + HR-centre null.\n",
        if (has_start) {
          paste0("End-point columns are shuffled within a stratum, start ",
                 "columns between strata.\n")
        } else {
          ""
        },
        "Each ranked block carries its own noise column: a variable earns ",
        "its place by beating that, not zero.",
        if (show_line) {
          paste0(" Dashed line: ", THRESHOLD, " per deer.")
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
      panel.grid = element_blank(),
      legend.position = "top",
      legend.justification = "left",
      axis.text.y = element_text(colour = "#0b0b0b", size = 10),
      axis.text.x = element_text(colour = "#898781", size = 9),
      axis.title.x = element_text(colour = "#52514e", size = 10),
      plot.title = element_text(colour = "#0b0b0b", face = "bold"),
      plot.subtitle = element_text(colour = "#52514e", size = 8.5)
    )

  out <- sprintf("plots/importance_xgb_%s_%s.png", season, config)
  ggsave(out, p, width = WIDTH, height = HEIGHT, dpi = 150, bg = "#fcfcfb")
  cat(sprintf("-> %s\n", out))
}

for (i in seq_len(nrow(selected))) {
  plot_season(selected$season[i], selected$config[i])
}
