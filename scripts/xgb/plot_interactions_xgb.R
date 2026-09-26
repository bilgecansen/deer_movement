#' @description
#' Which variable the selected shape split on under which — the interaction
#' structure of its fitted trees, for the seasons whose selected shape has
#' any.
#'
#' A tree using two variables is not yet an interaction. The unit is a
#' root-to-leaf path: the root splits on A, a child splits on B, so the
#' tree says how B matters depends on which side of A's threshold you are.
#' At max_depth 2 that parent-child pair is the whole of it. Cells are
#' observed pairs over what independence between the parent's and the
#' child's overall split frequencies would give, so 1.00 is what an
#' additive model looks like once it is made to grow depth-2 trees.
#'
#' Read it with three things in hand:
#'
#'   * A variable under ITSELF is a second cut on one curve, not an
#'     interaction, and it is usually the most enriched cell on the grid.
#'     Read the off-diagonal.
#'   * The ratio has no null of its own — on every fit so far the largest
#'     ratio anywhere is a noise column pairing with itself — so the pairs
#'     that matter are the ones carrying gain, printed below the figure,
#'     against what that block's noise column carries.
#'   * The gains in the dump are in-sample. What says interactions matter
#'     at all is the held-out shape comparison in plot_shapes_xgb.R; this
#'     figure only allocates that between pairs.
#'
#' A season whose selected shape is "main" has single-variable trees and so
#' no pairs; it is skipped, and saying so is part of the output.
#'
#' Input:  results/xgb/compare_<season>_<shape>.rds
#' Output: plots/interactions_xgb_<season>_<shape>.png
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
# Must match plot_shapes_xgb.R, or the figures describe a shape the
# selection did not pick.
THRESHOLD <- 3
COMPLEXITY <- c("main", "rsf", "start", "full")
# Cross-variable pairs listed per block, by share of the block's gain
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
  stop("No shape files in results/xgb/; run run_shapes_xgb.R first")
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
  shadow_gauss = "Noise", shadow_famd = "Noise",
  shadow_start = "Noise at start", famd1_end = "FAMD1",
  famd2_end = "FAMD2", famd3_end = "FAMD3", famd4_end = "FAMD4",
  famd5_end = "FAMD5", ndvi_start = "NDVI at start",
  wiscland_start = "Cover at start",
  forest_edge_start = "Forest edge at start"
)
BLOCK <- c(hab = "Habitat", famd = "FAMD", modifier = "Movement modifier")

dir.create("plots", showWarnings = FALSE)

for (i in seq_len(nrow(selected))) {
  season <- selected$season[i]
  config <- selected$config[i]
  if (config == "main") {
    cat(sprintf("%s: selected shape is 'main', single-variable trees, ",
                season))
    cat("no pairs to show\n")
    next
  }
  r <- readRDS(sprintf("results/xgb/compare_%s_%s.rds", season, config))

  cat(sprintf("\n########## %s / %s ##########\n", season, config))
  for (nm in names(r$structure)) {
    st <- r$structure[[nm]]
    cat(sprintf("\n--- %s: splits per variable ---\n", nm))
    print(as.data.frame(
      st$splits |>
        group_by(variable) |>
        summarise(splits = n(), as_root = sum(depth == 0, na.rm = TRUE),
                  gain = sum(gain, na.rm = TRUE), .groups = "drop") |>
        mutate(gain_share = gain / sum(gain)) |>
        select(-gain) |>
        arrange(desc(gain_share))
    ), row.names = FALSE, digits = 3)

    cat(sprintf("\n--- %s: top %d cross-variable pairs by gain ---\n", nm,
                TOP_N))
    print(head(as.data.frame(
      xgb_pair_lift(st$pairs) |>
        filter(parent != child) |>
        arrange(desc(gain_share)) |>
        select(parent, child, observed, ratio, gain_share)
    ), TOP_N), row.names = FALSE, digits = 3)
  }

  hm <- purrr::imap_dfr(r$structure, function(st, nm) {
    xgb_pair_lift(st$pairs) |> mutate(block = BLOCK[[nm]])
  }) |>
    mutate(parent = dplyr::coalesce(PRETTY[parent], parent),
           child = dplyr::coalesce(PRETTY[child], child))

  p <- ggplot(hm, aes(x = child, y = parent, fill = log2(ratio))) +
    geom_tile(colour = "#fcfcfb", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f\n%.0f%%", ratio,
                                  100 * gain_share)),
              size = 2.5, lineheight = 0.95, colour = "#0b0b0b") +
    facet_wrap(~block, scales = "free") +
    scale_fill_gradient2(low = "#2a78d6", mid = "#f2f1ea",
                         high = "#eb6834", midpoint = 0,
                         name = "log2 obs/exp") +
    labs(
      title = sprintf("Which splits the model put under which, %s (%s)",
                      season, config),
      subtitle = paste0(
        "Root variable on the vertical, the variable its child splits on ",
        "along the horizontal. Top number is how often\nthat combination ",
        "was chosen against independence; below it, that pair's share of ",
        "the block's gain.\n\nA variable under itself is a second cut on ",
        "one curve, not an interaction. The ratio has no null of its own, ",
        "so\nread the gain share, and read it against the noise column's."
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

  out <- sprintf("plots/interactions_xgb_%s_%s.png", season, config)
  ggsave(out, p, width = WIDTH, height = HEIGHT, dpi = 150,
         bg = "#fcfcfb")
  cat(sprintf("\n-> %s\n", out))
}
