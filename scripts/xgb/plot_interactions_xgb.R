#' @description
#' Which variable the selected model split on under which, and what each
#' combination is worth on steps the model never saw — for the seasons
#' whose selected model type has any pairs.
#'
#' A tree using two variables is not yet an interaction. The unit is a
#' root-to-leaf path: the root splits on A, a child splits on B, so the
#' tree says how B matters depends on which side of A's threshold you are.
#' At max_depth 2 that parent-child pair is the whole of it.
#'
#' Each cell is the held-out gain of the child splits on B under a root on
#' A: the held-out log-likelihood they add over leaving their node unsplit,
#' summed over the five folds (every step is held out once) and divided by
#' deer-years. The units are the importance plot's. A cell can be
#' negative: splits that fit noise in training make held-out predictions
#' worse. See xgb_heldout_split_gain().
#'
#' Read it with two things in hand:
#'
#'   * A variable under ITSELF is a second cut on one curve, not an
#'     interaction. Read the off-diagonal.
#'   * Gain under a parent is not all interaction: a child split also
#'     carries the child variable's own effect within that side of the
#'     parent. Read a pair against what the block's noise column carries.
#'
#' How often each pair was chosen, and its in-sample gain, are printed with
#' the tables, not drawn.
#'
#' A season whose selected model is "main" has single-variable trees and so
#' no pairs; it is skipped, and saying so is part of the output.
#'
#' Input:  results/xgb/compare_<season>_<type>.rds
#' Output: plots/interactions_xgb_<season>_<type>.png and .pdf
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
# Must match plot_models_xgb.R, or the figures describe a model the
# selection did not pick.
THRESHOLD <- 3
COMPLEXITY <- c("main", "rsf", "rsf_hr", "start", "full")
# Cross-variable pairs listed per block, by held-out gain
TOP_N <- 10
WIDTH <- 12
HEIGHT <- 6.5

# Load packages ---------------------------------------------------------------
library(tidyverse)

# helper functions
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
  ndvi_end = "NDVI", landcover = "Landcover",
  forest_edge_end = "Forest edge", elevation_end = "Elevation",
  northness_end = "Northness", eastness_end = "Eastness",
  HR_center_end = "HR center", shadow_gauss = "Noise",
  shadow_famd = "Noise", shadow_start = "Noise at start",
  famd1_end = "FAMD1", famd2_end = "FAMD2", famd3_end = "FAMD3",
  famd4_end = "FAMD4", famd5_end = "FAMD5", ndvi_start = "NDVI at start",
  wiscland_start = "Cover at start",
  forest_edge_start = "Forest edge at start", sl_ = "Step length",
  cos_ta = "cos(turn angle)"
)
BLOCK <- c(hab = "Habitat", famd = "FAMD", modifier = "Movement modifier")

dir.create("plots", showWarnings = FALSE)

for (i in seq_len(nrow(selected))) {
  season <- selected$season[i]
  config <- selected$config[i]
  if (config == "main") {
    cat(sprintf("%s: selected model is 'main', single-variable trees, ",
                season))
    cat("no pairs to show\n")
    next
  }
  r <- readRDS(sprintf("results/xgb/compare_%s_%s.rds", season, config))
  if (is.null(r$split_gain)) {
    stop(sprintf("%s / %s has no held-out split gain; refit it with ",
                 season, config),
         "fit_model_xgb.R")
  }
  n_deer <- r$n_deer_years
  pairs <- r$split_gain$pairs |>
    mutate(per_deer = heldout / n_deer)
  roots <- r$split_gain$roots |>
    mutate(per_deer = heldout / n_deer)

  cat(sprintf("\n########## %s / %s (per deer-year) ##########\n", season,
              config))
  for (b in unique(pairs$booster)) {
    cat(sprintf("\n--- %s: root splits by variable ---\n", b))
    print(as.data.frame(
      roots |>
        filter(booster == b) |>
        mutate(insample_share = insample / sum(insample)) |>
        select(variable, splits, per_deer, insample_share) |>
        arrange(desc(per_deer))
    ), row.names = FALSE, digits = 3)

    cat(sprintf("\n--- %s: top %d cross-variable pairs by held-out gain ---\n",
                b, TOP_N))
    print(head(as.data.frame(
      pairs |>
        filter(booster == b) |>
        mutate(insample_share = insample / sum(insample)) |>
        filter(parent != child) |>
        arrange(desc(per_deer)) |>
        select(parent, child, splits, per_deer, insample_share)
    ), TOP_N), row.names = FALSE, digits = 3)
  }

  hm <- pairs |>
    mutate(block = BLOCK[booster],
           parent = dplyr::coalesce(PRETTY[parent], parent),
           child = dplyr::coalesce(PRETTY[child], child),
           label = sub("^-(0[.]00)$", "\\1", sprintf("%.2f", per_deer)))
  lim <- max(abs(hm$per_deer))

  p <- ggplot(hm, aes(x = child, y = parent, fill = per_deer)) +
    geom_tile(colour = "#fcfcfb", linewidth = 0.6) +
    geom_text(aes(label = label), size = 3, colour = "#0b0b0b") +
    facet_wrap(~block, scales = "free") +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f2f1ea",
                         high = "#eb6834", midpoint = 0,
                         limits = c(-lim, lim),
                         name = "held-out gain\nper deer-year") +
    labs(
      title = sprintf("Which splits the model put under which, %s (%s)",
                      season, config),
      subtitle = paste0(
        "Root variable on the vertical, the variable its child splits on ",
        "along the horizontal. Each cell is the held-out log score\nper ",
        "deer-year those child splits add over leaving their node ",
        "unsplit, summed over five folds of whole steps. Orange helps on\n",
        "new steps, blue hurts.\n\nA variable under itself is a second ",
        "cut on one curve, not an interaction, and a child split also ",
        "carries the child's own effect\non that side of its parent. Read ",
        "each pair against the noise column's."
      ),
      x = "child split", y = "root split"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.background = element_rect(fill = "#fcfcfb", colour = NA),
      panel.background = element_rect(fill = "#fcfcfb", colour = NA),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 40, hjust = 1, colour = "#0b0b0b"),
      axis.text.y = element_text(colour = "#0b0b0b"),
      strip.text = element_text(colour = "#0b0b0b", face = "bold",
                                hjust = 0),
      plot.title = element_text(colour = "#0b0b0b", face = "bold"),
      plot.subtitle = element_text(colour = "#52514e", size = 8)
    )

  out <- sprintf("plots/interactions_xgb_%s_%s", season, config)
  ggsave(paste0(out, ".png"), p, width = WIDTH, height = HEIGHT,
         dpi = 150, bg = "#fcfcfb")
  ggsave(paste0(out, ".pdf"), p, width = WIDTH, height = HEIGHT,
         bg = "#fcfcfb")
  cat(sprintf("\n-> %s.png, %s.pdf\n", out, out))
}
