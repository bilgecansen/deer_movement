#' @description
#' Plot the interaction log score of each season's SELECTED model type: one
#' figure per season, one row per pair of variables, laid out like the
#' importance plot and read against the same threshold.
#'
#' A tree using two variables is not yet an interaction. The unit is a
#' root-to-leaf path: the root splits on A, a child splits on B, so the
#' tree says how B matters depends on which side of A's threshold you are.
#' At max_depth 2 that parent-child pair is the whole of it.
#'
#' A pair's interaction log score is the held-out log score its child
#' splits gain over leaving their node unsplit (xgb_heldout_split_gain()),
#' with B under A and A under B added together — both are the same pair —
#' summed over the five folds, so every step counts once, and divided by
#' deer-years. The units are the importance plot's, which measures the
#' drop when a variable is shuffled; this measures the gain a pair's splits
#' bring.
#'
#' A child split also carries the child variable's own effect on that side
#' of its parent, so the score is an upper bound on the interaction, not
#' the interaction itself. A pair has to clear THRESHOLD per deer on that
#' upper bound to be worth carrying into the individual-deer models. The
#' model comparison asks whether interactions are worth looking at at all;
#' this asks whether any one of them is.
#'
#' A variable under itself is a second cut on one curve, not an
#' interaction, and is left off. Pairs with a noise column mark where
#' "carries nothing" sits.
#'
#' A season whose selected model is "main" has single-variable trees and so
#' no pairs, and one whose selected model is the null type has no ranked
#' trees at all; either is skipped, and saying so is part of the output.
#'
#' Input:  results/xgb/compare_<season>_<type>.rds
#' Output: plots/interactions_xgb_<season>_<type>.png and .pdf
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
# Must match plot_models_xgb.R, or the figures describe a model the
# selection did not pick.
THRESHOLD <- 3
COMPLEXITY <- c("null", "main", "rsf", "rsf_hr", "start", "full")
WIDTH <- 9
# Height grows with the number of pairs
HEIGHT_BASE <- 2.2
HEIGHT_PER_ROW <- 0.2

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
  cos_ta = "cos(turn angle)", day_of_season = "Day of season"
)
BLOCK <- c(hab = "Habitat", famd = "FAMD", modifier = "Movement modifier")

dir.create("plots", showWarnings = FALSE)

plot_season <- function(season, config) {
  r <- readRDS(sprintf("results/xgb/compare_%s_%s.rds", season, config))
  if (is.null(r$split_gain)) {
    stop(sprintf("%s / %s has no held-out split gain; refit it with ",
                 season, config),
         "fit_model_xgb.R")
  }
  # B under A and A under B are one pair: name it by its two variables in
  # a fixed order and add the two directions.
  df <- r$split_gain$pairs |>
    filter(parent != child) |>
    mutate(a = pmin(parent, child), b = pmax(parent, child)) |>
    group_by(booster, a, b) |>
    summarise(value = sum(heldout) / r$n_deer_years, .groups = "drop") |>
    mutate(
      block = factor(BLOCK[booster], levels = BLOCK),
      label = paste(dplyr::coalesce(PRETTY[a], a), "×",
                    dplyr::coalesce(PRETTY[b], b)),
      kind = dplyr::case_when(
        grepl("^shadow", a) | grepl("^shadow", b) ~ "noise",
        booster == "modifier" ~ "start",
        TRUE ~ "pair"
      )
    ) |>
    arrange(value)
  stopifnot(!anyDuplicated(df$label))
  df$label <- factor(df$label, levels = df$label)

  cat(sprintf("\n########## %s / %s: interaction log score per deer ######",
              season, config))
  cat("####\n")
  print(as.data.frame(df |> arrange(block, desc(value)) |>
                        select(block, label, value)),
        row.names = FALSE, digits = 3)

  # Pale stripe behind every second row of each block, to carry the eye
  # across, with the gridlines drawn over the stripes at the axis breaks.
  bands <- df |>
    group_by(block) |>
    mutate(row = rank(value, ties.method = "first")) |>
    filter(row %% 2 == 0) |>
    ungroup()
  breaks <- scales::breaks_extended(6)(range(c(0, THRESHOLD, df$value)))

  p <- ggplot(df, aes(y = label, x = value, colour = kind)) +
    geom_rect(data = bands, inherit.aes = FALSE,
              aes(ymin = row - 0.5, ymax = row + 0.5),
              xmin = -Inf, xmax = Inf, fill = "#edece5") +
    geom_vline(xintercept = breaks, colour = "#d8d6cd", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "#c3c2b7", linewidth = 0.5) +
    geom_vline(xintercept = THRESHOLD, colour = "#2a78d6",
               linewidth = 0.5, linetype = "22") +
    geom_segment(aes(x = 0, xend = value, yend = label), linewidth = 0.9) +
    geom_point(size = 4) +
    facet_grid(block ~ ., scales = "free_y", space = "free_y") +
    scale_colour_manual(
      values = c(pair = "#eb6834", start = "#2a78d6", noise = "#898781"),
      breaks = c("pair", "start", "noise"),
      labels = c(pair = "pair of variables",
                 start = "movement x start of step",
                 noise = "pair with a noise column"),
      name = NULL
    ) +
    scale_x_continuous(breaks = breaks, labels = scales::label_comma(),
                       expand = expansion(mult = c(0.05, 0.05))) +
    expand_limits(x = THRESHOLD) +
    labs(
      title = sprintf("Interaction log score, pooled %s (%s model)",
                      season, config),
      subtitle = paste0(
        r$n_deer_years, " deer-years (",
        formatC(r$n_steps, format = "d", big.mark = ","),
        " observed steps). Held-out log score the trees gain by splitting ",
        "on one variable\nunder the other, both orders added, over five ",
        "random folds of whole steps. Each split also carries its own\n",
        "variable's effect on that side of its parent, so this is an upper ",
        "bound on the interaction. Dashed line: ", THRESHOLD, " per deer."
      ),
      x = "Gain in total log score per deer",
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
      strip.text.y = element_text(colour = "#0b0b0b", face = "bold",
                                  angle = 0, hjust = 0),
      plot.title = element_text(colour = "#0b0b0b", face = "bold"),
      plot.subtitle = element_text(colour = "#52514e", size = 8.5),
      panel.spacing.y = grid::unit(10, "pt")
    )

  height <- HEIGHT_BASE + HEIGHT_PER_ROW * nrow(df)
  out <- sprintf("plots/interactions_xgb_%s_%s", season, config)
  ggsave(paste0(out, ".png"), p, width = WIDTH, height = height, dpi = 150,
         bg = "#fcfcfb")
  ggsave(paste0(out, ".pdf"), p, width = WIDTH, height = height,
         bg = "#fcfcfb")
  cat(sprintf("-> %s.png, %s.pdf\n", out, out))
}

for (i in seq_len(nrow(selected))) {
  if (selected$config[i] %in% c("main", "null")) {
    cat(sprintf("%s: selected model is '%s', no pairs to show\n",
                selected$season[i], selected$config[i]))
    next
  }
  plot_season(selected$season[i], selected$config[i])
}
