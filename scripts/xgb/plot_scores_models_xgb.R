#' @description
#' Violin plots of the two path filters' scores, HR overlap and svf
#' agreement, model by model: the season's two xgb deer picks, the GAM null
#' and the Gilbertson et al. (2025) model, one point per deer-year. Where
#' plot_scores_xgb.R keeps the better xgb pick per deer-year, this one keeps
#' every model apart, so the habitat models can be set beside the null. The
#' number under each violin is the deer-years at or above the filter's
#' threshold; a deer-year whose ctmm fit failed has no HR overlap and is
#' left out of that violin, its count printed.
#'
#' Inputs:  filters/xgb/scores_xgb_<tag>.rds and scores_gam_<tag>.rds
#'          (score_sims_xgb.R, xgb and gam), <tag> = xgb_scores_tag(SET,
#'          REACH)
#' Output:  plots/scores_models_<tag>.png, and a summary printed
#'
#' Configuration: edit the block below before running. The set can also be
#' given on the command line: Rscript scripts/xgb/plot_scores_models_xgb.R
#' test

# Configuration ---------------------------------------------------------------
SET <- "train"
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1) {
  SET <- .args[1]
}
stopifnot(SET %in% c("test", "train"))
REACH <- "gamma99"
# The filters' thresholds, drawn as a dashed line
THRESHOLD <- c(hr_overlap = 0.8, svf_agree = 0.8)

# Load packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
})

# helper functions
source("scripts/helper_functions.R")

tag <- xgb_scores_tag(SET, REACH)
xgb <- readRDS(sprintf("filters/xgb/scores_xgb_%s.rds", tag))
gam <- readRDS(sprintf("filters/xgb/scores_gam_%s.rds", tag))
models <- c("xgb best", "xgb simplest", "GAM null", "Gilbertson")
# A type that is both a season's best and its simplest pick counts as both
scores <- dplyr::bind_rows(
  dplyr::mutate(dplyr::filter(xgb, best), model = "xgb best"),
  dplyr::mutate(dplyr::filter(xgb, simplest), model = "xgb simplest"),
  dplyr::mutate(gam, model = dplyr::recode(config,
                                           null_gam = "GAM null",
                                           gilbertson_gam = "Gilbertson"))
) |>
  dplyr::select(key, season, model, hr_overlap, svf_agree) |>
  tidyr::pivot_longer(c(hr_overlap, svf_agree), names_to = "score",
                      values_to = "value") |>
  dplyr::mutate(model = factor(model, models),
                threshold = THRESHOLD[score])
stopifnot(!anyNA(scores$model))

summ <- scores |>
  dplyr::group_by(score, season, model) |>
  dplyr::summarise(n = dplyr::n(), missing = sum(is.na(value)),
                   median = stats::median(value, na.rm = TRUE),
                   pass = sum(value >= threshold, na.rm = TRUE),
                   .groups = "drop")
cat(sprintf("=== %s: median score and deer-years at or above 0.8 ===\n",
            SET))
print(as.data.frame(summ |>
                      dplyr::mutate(cell = sprintf("%.3f  %d/%d%s", median,
                                                   pass, n,
                                                   ifelse(missing > 0,
                                                          sprintf(" (%d NA)",
                                                                  missing),
                                                          ""))) |>
                      dplyr::select(score, season, model, cell) |>
                      tidyr::pivot_wider(names_from = model,
                                         values_from = cell)),
      row.names = FALSE)
cat("\nall seasons:\n")
print(as.data.frame(scores |>
                      dplyr::group_by(score, model) |>
                      dplyr::summarise(
                        median = stats::median(value, na.rm = TRUE),
                        pass = sum(value >= threshold, na.rm = TRUE),
                        n = dplyr::n(), missing = sum(is.na(value)),
                        .groups = "drop"
                      )),
      digits = 3, row.names = FALSE)

# Plot -----------------------------------------------------------------------
cols <- c("xgb best" = "#1B9E77", "xgb simplest" = "#66A61E",
          "GAM null" = "#7570B3", "Gilbertson" = "#D95F02")
score_lab <- c(hr_overlap = "HR overlap", svf_agree = "svf agreement")
season_lab <- c(br = "Breeding", fa = "Fawning", nb = "Non-breeding",
                pf = "Post-fawning")
labels <- summ |>
  dplyr::mutate(label = sprintf("%d of %d", pass, n))
g <- ggplot(dplyr::filter(scores, !is.na(value)),
            aes(model, value, fill = model)) +
  geom_hline(aes(yintercept = threshold), linetype = 2, colour = "grey30") +
  geom_violin(scale = "width", width = 0.85, linewidth = 0.3,
              colour = "grey25") +
  geom_jitter(width = 0.12, height = 0, size = 0.5, alpha = 0.35) +
  stat_summary(fun = stats::median, geom = "point", shape = 23, size = 2.2,
               fill = "white", colour = "black") +
  geom_text(data = labels, aes(model, -0.06, label = label),
            inherit.aes = FALSE, size = 2.9) +
  scale_fill_manual(NULL, values = cols, guide = "none") +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(-0.08, 1)) +
  facet_grid(score ~ season,
             labeller = labeller(score = score_lab, season = season_lab)) +
  labs(x = NULL, y = "Score per deer-year (1 is best)",
       title = sprintf(paste("HR overlap and svf agreement of the simulated",
                             "paths, model by model, %s"), SET),
       subtitle = paste("Dashed line: the filter's threshold, 0.8. White",
                        "diamond: median. Numbers: deer-years at or above",
                        "the threshold."),
       caption = paste("A deer-year whose ctmm fit failed has no HR overlap",
                       "and is left out of the violin; it still counts in",
                       "the total under it.")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
out <- sprintf("plots/scores_models_%s.png", tag)
ggsave(out, g, width = 14, height = 8, dpi = 140)
cat(sprintf("-> %s\n", out))
