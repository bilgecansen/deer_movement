#' @description
#' Compare the seven deer model types (the animal id in every block, see
#' make_xgb_specs()) within each season and pick one.
#'
#' The rule: take the SIMPLEST model type whose held-out score is within
#' THRESHOLD log units per deer of the season's best type. The threshold is
#' the pipeline's own gate 3, where a model must beat the null by 3 log
#' units for a deer.
#'
#' Complexity runs null < main < rsf < rsf_hr < start < full < full_hr. The
#' null type is movement and the home-range centre alone, so if no other
#' type beats it by the threshold, the season selects it: no habitat model
#' is worth it there. rsf_hr is rsf with one more
#' variable allowed to interact, so it sits just above rsf, and full_hr
#' stands to full as rsf_hr to rsf, so it comes last. rsf_hr's place
#' against start and full, which it is not nested in, is a judgment call,
#' as is rsf against start — fewer trees but interactions, against more
#' trees but additive.
#'
#' Per-deer favours seasons with longer tracks: nb carries about 390 steps
#' per deer-year against 140-200 in the other seasons, so the same per-step
#' effect reads roughly twice as large there. The per-100-step column is
#' printed alongside so that is visible rather than buried.
#'
#' Input:  results/xgb/compare_<season>_<type>_deer.rds
#' Output: plots/compare_models_deer.png and .pdf, and the selection printed
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
# Log units per deer. A model type within this much of the season's best is
# a candidate; the simplest candidate wins.
THRESHOLD <- 3
# Simplest first
COMPLEXITY <- c("null", "main", "rsf", "rsf_hr", "start", "full", "full_hr")
WIDTH <- 15
HEIGHT <- 9

# Load packages ---------------------------------------------------------------
library(tidyverse)

# helper functions
source("scripts/helper_functions.R")

files <- list.files("results/xgb", pattern = "^compare_.*_deer[.]rds$",
                    full.names = TRUE)
if (!length(files)) {
  stop("No model files in results/xgb/; run run_models_xgb.R first")
}
r <- purrr::map_dfr(files, function(f) {
  x <- readRDS(f)
  tibble(season = x$season, config = x$config, trees = x$n_trees,
         in_sample = x$loglik_in, held_out = x$ll_cv,
         n_deer = x$n_deer_years, n_steps = x$n_steps)
})

r <- r |>
  mutate(type = sub("_deer$", "", config)) |>
  group_by(season) |>
  mutate(delta = held_out - max(held_out),
         per_deer = delta / n_deer,
         per_100 = 100 * delta / n_steps,
         within = per_deer > -THRESHOLD,
         selected = type == COMPLEXITY[min(match(type[within],
                                                 COMPLEXITY))]) |>
  ungroup()

cat("=== held-out score by model type, against the season's best ===\n")
print(as.data.frame(
  r |>
    arrange(season, match(type, COMPLEXITY)) |>
    select(season, config, trees, in_sample, held_out, per_deer, per_100,
           selected)
), row.names = FALSE, digits = 5)

cat("\n=== selection ===\n")
print(as.data.frame(
  r |>
    filter(selected) |>
    arrange(season) |>
    select(season, config, per_deer)
), row.names = FALSE, digits = 4)

# Each effect read twice, which is what the 2x2 is for
w <- r |>
  select(season, type, held_out, n_deer) |>
  pivot_wider(names_from = type, values_from = held_out)
if (all(c("full", "rsf", "start", "main") %in% names(w))) {
  cat("\n=== each effect read twice, per deer ===\n")
  print(as.data.frame(w |> transmute(
    season,
    ranked_interactions_a = (full - start) / n_deer,
    ranked_interactions_b = (rsf - main) / n_deer,
    start_block_a = (full - rsf) / n_deer,
    start_block_b = (start - main) / n_deer
  ) |> arrange(season)), row.names = FALSE, digits = 4)
}

# Figures ---------------------------------------------------------------------
LABEL <- c(
  full_hr = "Movement, HR and RSF Interactions",
  full = "Movement and RSF Interactions",
  rsf = "RSF Interactions",
  rsf_hr = "RSF and HR Interactions",
  start = "Movement Interactions",
  main = "Main effects",
  null = "Null (movement + HR)"
)

per_season <- r |>
  distinct(season, n_deer, n_steps) |>
  mutate(spd = round(n_steps / n_deer))
strip <- stats::setNames(
  sprintf("%s  -  %d deer-years, %s steps (%d per deer)",
          per_season$season, per_season$n_deer,
          formatC(per_season$n_steps, format = "d", big.mark = ","),
          per_season$spd),
  as.character(per_season$season)
)

dir.create("plots", showWarnings = FALSE)

df <- r |>
  mutate(label = factor(LABEL[type], levels = LABEL[COMPLEXITY]),
         season = factor(season, levels = c("br", "fa", "nb", "pf")))

p <- ggplot(df, aes(x = per_deer, y = label)) +
  geom_vline(xintercept = 0, colour = "#c3c2b7", linewidth = 0.5) +
  geom_vline(xintercept = -THRESHOLD, colour = "#2a78d6",
             linewidth = 0.5, linetype = "22") +
  geom_segment(aes(x = 0, xend = per_deer, yend = label,
                   colour = selected),
               linewidth = 0.9) +
  geom_point(aes(colour = selected), size = 4) +
  geom_text(aes(label = sprintf("%.2f", per_deer)), hjust = 1.35,
            size = 3, colour = "#52514e") +
  scale_colour_manual(values = c(`FALSE` = "#d9c4bb", `TRUE` = "#eb6834"),
                      guide = "none") +
  facet_wrap(~season, ncol = 2, scales = "free_x",
             labeller = labeller(season = strip)) +
  scale_x_continuous(expand = expansion(mult = c(0.22, 0.06))) +
  labs(
    title = sprintf(
      paste("Simplest model within %g log score per deer of the best,",
            "with the deer id (the animal id in every block, trees 3 deep)"),
      THRESHOLD
    ),
    subtitle = paste0(
      "Five random folds of whole steps, 500 rounds per booster. The best ",
      "model sits at 0; the others show what they\ngive up, in held-out ",
      "log score per deer-year. Anything right of the dashed line is ",
      "within ", THRESHOLD, " of the best,\nand the simplest of those is ",
      "the selected model, drawn solid.\n\n",
      "Per deer-year favours seasons with longer tracks, so the ",
      "per-100-step gaps are printed with the table."
    ),
    x = "Held-out log score per deer, against the season's best model",
    y = NULL
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
    strip.text = element_text(colour = "#0b0b0b", face = "bold",
                              hjust = 0),
    plot.title = element_text(colour = "#0b0b0b", face = "bold"),
    plot.subtitle = element_text(colour = "#52514e", size = 8.5),
    panel.spacing.y = grid::unit(12, "pt")
  )

out <- "plots/compare_models_deer"
ggsave(paste0(out, ".png"), p, width = WIDTH, height = HEIGHT, dpi = 150,
       bg = "#fcfcfb")
ggsave(paste0(out, ".pdf"), p, width = WIDTH, height = HEIGHT,
       bg = "#fcfcfb")
cat(sprintf("\n-> %s.png, %s.pdf\n", out, out))
