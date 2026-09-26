#' @description
#' Compare the four model shapes within each season and pick one.
#'
#' The rule: take the SIMPLEST shape whose held-out score is within
#' THRESHOLD log units per deer of the best shape in that season. The
#' threshold is the pipeline's own gate 3, where a model must beat the null
#' by 3 log units for a deer.
#'
#' Complexity runs main < rsf < start < full. The ranking between rsf and
#' start is arguable — fewer trees but interactions, against more trees but
#' additive — and it has not mattered yet, since the within-threshold set
#' has been either all four shapes or just full and rsf.
#'
#' Per-deer favours seasons with longer tracks: nb and pf carry about 400
#' steps per deer against fa's 164, so the same per-step effect reads
#' roughly 2.5x larger there. The per-100-step column is printed alongside
#' so that is visible rather than buried.
#'
#' Input:  results/xgb/compare_<season>_<shape>.rds
#' Output: plots/compare_shapes.png, and the selection printed
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
# Log units per deer. A shape within this much of the best is a candidate;
# the simplest candidate wins.
THRESHOLD <- 3
# Simplest first
COMPLEXITY <- c("main", "rsf", "start", "full")
WIDTH <- 9.5
HEIGHT <- 8

# Load packages ---------------------------------------------------------------
library(tidyverse)

# helper functions
source("scripts/helper_functions.R")

files <- list.files("results/xgb", pattern = "^compare_.*[.]rds$",
                    full.names = TRUE)
if (!length(files)) {
  stop("No shape files in results/xgb/; run run_shapes_xgb.R first")
}
r <- purrr::map_dfr(files, function(f) {
  x <- readRDS(f)
  tibble(season = x$season, config = x$config, trees = x$n_trees,
         in_sample = x$loglik_in, held_out = x$ll_cv,
         n_deer = x$n_deer_years, n_steps = x$n_steps)
})

r <- r |>
  group_by(season) |>
  mutate(delta = held_out - max(held_out),
         per_deer = delta / n_deer,
         per_100 = 100 * delta / n_steps,
         within = per_deer > -THRESHOLD,
         selected = config == COMPLEXITY[min(match(config[within],
                                                   COMPLEXITY))]) |>
  ungroup()

cat("=== held-out score by shape ===\n")
print(as.data.frame(
  r |>
    mutate(config = factor(config, levels = COMPLEXITY)) |>
    arrange(season, config) |>
    select(season, config, trees, in_sample, held_out, per_deer, per_100,
           selected)
), row.names = FALSE, digits = 5)

cat("\n=== selection ===\n")
print(as.data.frame(
  r |>
    filter(selected) |>
    select(season, config, per_deer)
), row.names = FALSE, digits = 4)

# Each effect read twice, which is what the 2x2 is for.
w <- r |>
  select(season, config, held_out, n_deer) |>
  pivot_wider(names_from = config, values_from = held_out)
if (all(c("full", "rsf", "start", "main") %in% names(w))) {
  cat("\n=== each effect read twice, per deer ===\n")
  print(as.data.frame(w |> transmute(
    season,
    ranked_interactions_a = (full - start) / n_deer,
    ranked_interactions_b = (rsf - main) / n_deer,
    start_block_a = (full - rsf) / n_deer,
    start_block_b = (start - main) / n_deer
  )), row.names = FALSE, digits = 4)
}

# Figure ----------------------------------------------------------------------
LABEL <- c(
  full = "Full\nstart block + RSF interactions",
  rsf = "RSF interaction\nRSF interactions only",
  start = "Start only\nstart block only",
  main = "Main effects\nneither"
)
df <- r |>
  mutate(config = factor(config, levels = COMPLEXITY),
         label = factor(LABEL[as.character(config)],
                        levels = LABEL[COMPLEXITY]),
         season = factor(season, levels = c("fa", "nb", "pf")))

per_season <- df |>
  distinct(season, n_deer, n_steps) |>
  mutate(spd = round(n_steps / n_deer))
strip <- stats::setNames(
  sprintf("%s  -  %d deer-years, %s steps (%d per deer)",
          per_season$season, per_season$n_deer,
          formatC(per_season$n_steps, format = "d", big.mark = ","),
          per_season$spd),
  as.character(per_season$season)
)

p <- ggplot(df, aes(x = per_deer, y = label)) +
  geom_vline(xintercept = 0, colour = "#c3c2b7", linewidth = 0.5) +
  geom_vline(xintercept = -THRESHOLD, colour = "#2a78d6",
             linewidth = 0.5, linetype = "22") +
  geom_segment(aes(x = 0, xend = per_deer, yend = label),
               colour = "#eb6834", linewidth = 0.9) +
  geom_point(aes(colour = selected), size = 4) +
  geom_text(aes(label = sprintf("%.2f", per_deer)), hjust = 1.35,
            size = 3, colour = "#52514e") +
  scale_colour_manual(values = c(`FALSE` = "#d9c4bb", `TRUE` = "#eb6834"),
                      guide = "none") +
  facet_wrap(~season, ncol = 1, scales = "free_x",
             labeller = labeller(season = strip)) +
  scale_x_continuous(expand = expansion(mult = c(0.22, 0.06))) +
  labs(
    title = sprintf(
      "Simplest shape within %g log score per deer of the best", THRESHOLD
    ),
    subtitle = paste0(
      "Five random folds of whole steps, 500 rounds per booster. The best ",
      "shape sits at 0; the others show what they\ngive up, in held-out ",
      "log score per deer-year. Anything right of the dashed line is ",
      "within ", THRESHOLD, " of the best,\nand the simplest of those is ",
      "the selected shape, drawn solid.\n\n",
      "Per deer-year favours seasons with longer tracks, so the ",
      "per-100-step gaps are printed with the table."
    ),
    x = "Held-out log score per deer, against the best shape", y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "#fcfcfb", colour = NA),
    panel.background = element_rect(fill = "#fcfcfb", colour = NA),
    panel.grid.major.x = element_line(colour = "#e1e0d9", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(colour = "#0b0b0b", size = 9),
    axis.text.x = element_text(colour = "#898781", size = 9),
    axis.title.x = element_text(colour = "#52514e", size = 10),
    strip.text = element_text(colour = "#0b0b0b", face = "bold", hjust = 0),
    plot.title = element_text(colour = "#0b0b0b", face = "bold"),
    plot.subtitle = element_text(colour = "#52514e", size = 8.5),
    panel.spacing.y = grid::unit(12, "pt")
  )

dir.create("plots", showWarnings = FALSE)
ggsave("plots/compare_shapes.png", p, width = WIDTH, height = HEIGHT,
       dpi = 150, bg = "#fcfcfb")
cat("\n-> plots/compare_shapes.png\n")
