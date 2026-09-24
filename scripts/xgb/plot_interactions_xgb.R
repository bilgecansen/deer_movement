#' @description
#' Show which variable the model split on under which — the interaction
#' structure of the fitted trees, for the record that there is little there.
#'
#' A tree using two variables is not yet an interaction. The unit is a
#' root-to-leaf path: the root splits on A, a child splits on B, so the tree
#' says how B matters depends on which side of A's threshold you are. At
#' max_depth 2 that parent-child pair is the whole of it.
#'
#' Cells are observed pairs over what independence between the parent's and
#' the child's overall split frequencies would give. 1.00 is what an
#' additive model looks like once it is made to grow depth-2 trees.
#'
#' Two things have to be read alongside the numbers:
#'
#'   * The noise column is in the grid. On every fit so far the single
#'     largest ratio anywhere is noise pairing with itself, so the ratio has
#'     no usable null — a greedy learner always finds a second split. Where
#'     the real pairs separate is gain share, printed below the figure: the
#'     noise column takes a couple of percent of gain, NDVI with forest edge
#'     and NDVI with landcover take ten or more.
#'   * A variable splitting on itself is refinement of one curve, not an
#'     interaction, and it is the most enriched pattern throughout. Read the
#'     off-diagonal.
#'
#' Nothing here says a combination holds up out of sample; the gains in the
#' dump are in-sample. The held-out test is to refit with the trees confined
#' to one variable and compare the folds' score — on fa 2021 that bought 2.4
#' log units against 186 from the same trees without pairing.
#'
#' Input:  results/xgb/season_xgb_<season>.rds
#' Output: plots/interactions_xgb_<season>.png
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
# Pairs listed in the printed table, most enriched first
TOP_N <- 12
WIDTH <- 12
HEIGHT <- 6.5

# Load packages ---------------------------------------------------------------
library(tidyverse)

# helper functions (xgb_pair_lift)
source("scripts/helper_functions.R")

files <- list.files("results/xgb", pattern = "^season_xgb_.*[.]rds$",
                    full.names = TRUE)
if (!length(files)) {
  stop("No season files in results/xgb/; run fit_season_xgb.R first")
}

PRETTY <- c(
  ndvi_end = "NDVI", landcover = "Landcover",
  forest_edge_end = "Forest edge", elevation_end = "Elevation",
  northness_end = "Northness", eastness_end = "Eastness",
  oak_mast_end = "Oak mast", oak_dist_end = "Oak distance",
  shadow_gauss = "Noise", famd1_end = "FAMD1", famd2_end = "FAMD2",
  famd3_end = "FAMD3", famd4_end = "FAMD4", famd5_end = "FAMD5"
)

dir.create("plots", showWarnings = FALSE)

plot_season <- function(path) {
  r <- readRDS(path)
  lift <- purrr::imap_dfr(r$structure, function(st, nm) {
    xgb_pair_lift(st$pairs) |>
      mutate(booster = ifelse(nm == "hab", "Habitat", "FAMD"))
  })

  cat(sprintf("\n=== %s: splits per variable ===\n", r$season))
  print(as.data.frame(
    purrr::imap_dfr(r$structure, function(st, nm) {
      st$splits |>
        group_by(variable) |>
        summarise(splits = n(), as_root = sum(depth == 0, na.rm = TRUE),
                  gain = sum(gain, na.rm = TRUE), .groups = "drop") |>
        mutate(booster = nm, gain_share = gain / sum(gain)) |>
        select(booster, variable, splits, as_root, gain_share)
    }) |>
      arrange(booster, desc(splits))
  ), row.names = FALSE, digits = 3)

  cat(sprintf("\n=== %s: top %d pairs by enrichment ===\n", r$season,
              TOP_N))
  print(head(as.data.frame(lift |> arrange(desc(ratio))), TOP_N),
        row.names = FALSE, digits = 3)

  df <- lift |>
    mutate(parent = dplyr::coalesce(PRETTY[parent], parent),
           child = dplyr::coalesce(PRETTY[child], child))

  p <- ggplot(df, aes(x = child, y = parent, fill = log2(ratio))) +
    geom_tile(colour = "#fcfcfb", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f", ratio)), size = 2.7,
              colour = "#0b0b0b") +
    facet_wrap(~booster, scales = "free") +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f2f1ea",
                         high = "#eb6834", midpoint = 0,
                         name = "log2 obs/exp") +
    labs(
      title = sprintf("Which splits the model put under which, pooled %s",
                      r$season),
      subtitle = paste0(
        "Root variable on the vertical, the variable its child splits on ",
        "along the horizontal. Values are how often that\ncombination was ",
        "chosen against what each variable's overall taste for splitting ",
        "would give. 1.00 is what an\nadditive model looks like once it ",
        "is made to grow depth-2 trees. A variable under itself is ",
        "refinement of one\ncurve, not an interaction. Counts only, from ",
        "in-sample gains - nothing here says a pair holds up out of ",
        "sample."
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

  out <- sprintf("plots/interactions_xgb_%s.png", r$season)
  ggsave(out, p, width = WIDTH, height = HEIGHT, dpi = 150, bg = "#fcfcfb")
  cat(sprintf("\n-> %s\n", out))
}

for (f in files) {
  plot_season(f)
}
