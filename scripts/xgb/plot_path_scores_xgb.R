#' @description
#' Box plots of one path score per deer-year (score_sims_xgb.R) for the
#' season's two final xgb deer models, the Gilbertson et al. (2025)
#' habitat model (gam_hab) and the GAM null (gam_hrc: movement and the home
#' range center). SCORE picks the score:
#'   "hr_overlap"  the Bhattacharyya coefficient between the real deer's
#'                 AKDE and the mean AKDE of the model's simulated paths
#'   "svf_agree"   the agreement of the real path's variogram with the
#'                 mean of the simulated paths' variograms
#' One figure per set, one panel per season. The dashed line is the score's
#' filter threshold (filter_paths_xgb.R); under each box, the deer-years at
#' or above it. A deer-year without a value (a failed ctmm fit) is left out
#' of its model's box and count.
#'
#' Inputs:  filters/xgb/scores_xgb_<tag>.rds and scores_gam_<tag>.rds
#'          (score_sims_xgb.R, xgb and gam); <tag> = xgb_scores_tag(set,
#'          REACH)
#' Output:  plots/<SCORE>_<tag>.png, one per set
#'
#' Configuration: edit the block below before running. The score can also
#' be given on the command line:
#'   Rscript scripts/xgb/plot_path_scores_xgb.R svf_agree

# Configuration ---------------------------------------------------------------
SCORE <- "hr_overlap"
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1) {
  SCORE <- .args[1]
}
REACH <- "gamma99"
# Per score: the filter threshold, the title and the y-axis title
SCORES <- list(
  hr_overlap = list(
    threshold = 0.8,
    title = "Do observed and simulated home ranges overlap?",
    y = "HR overlap\n(Bhattacharyya coefficient)"
  ),
  svf_agree = list(
    threshold = 0.8,
    title = "Do observed and simulated variograms agree?",
    y = "svf agreement"
  )
)
stopifnot(SCORE %in% names(SCORES))
spec <- SCORES[[SCORE]]
# The share of the most extreme values of each panel left outside the axis
TRIM <- 0.01

# Load packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# helper functions
source("scripts/helper_functions.R")

models <- c("xgb best", "xgb simplest", "gam_hab", "gam_hrc")
read_set <- function(set) {
  tag <- xgb_scores_tag(set, REACH)
  xgb <- readRDS(sprintf("filters/xgb/scores_xgb_%s.rds", tag))
  gam <- readRDS(sprintf("filters/xgb/scores_gam_%s.rds", tag))
  # A type that is both a season's best and its simplest pick counts as both
  dplyr::bind_rows(
    dplyr::mutate(dplyr::filter(xgb, best), model = "xgb best"),
    dplyr::mutate(dplyr::filter(xgb, simplest), model = "xgb simplest"),
    dplyr::mutate(gam, model = dplyr::recode(config,
                                             null_gam = "gam_hrc",
                                             gilbertson_gam = "gam_hab"))
  ) |>
    dplyr::transmute(set = set, key, season,
                     model = factor(model, models), value = .data[[SCORE]])
}
d <- dplyr::bind_rows(read_set("train"), read_set("test")) |>
  dplyr::filter(!is.na(value))
stopifnot(!anyNA(d$model))

cat(sprintf("=== %s per deer-year, median and deer-years >= %g ===\n",
            SCORE, spec$threshold))
print(as.data.frame(d |>
                      dplyr::group_by(set, season, model) |>
                      dplyr::summarise(deer_years = dplyr::n(),
                                       median = stats::median(value),
                                       at_least = sum(value >=
                                                        spec$threshold),
                                       .groups = "drop")),
      digits = 3, row.names = FALSE)

# Plots ----------------------------------------------------------------------
# One figure per set, one panel per season, each on its own axis
cols <- c("xgb best" = "#1B9E77", "xgb simplest" = "#66A61E",
          "gam_hab" = "#D95F02", "gam_hrc" = "#7570B3")
season_lab <- c(br = "Breeding", fa = "Fawning", nb = "Non-breeding",
                pf = "Post-fawning")
panel <- function(s, season) {
  ds <- d[d$set == s & d$season == season, ]
  lims <- range(c(spec$threshold, stats::quantile(ds$value,
                                                  c(TRIM, 1 - TRIM))))
  # Each model's name, with its deer-years at or above the threshold
  counts <- ds |>
    dplyr::group_by(model) |>
    dplyr::summarise(label = sprintf("%s\n%d of %d ≥ %g", model[1],
                                     sum(value >= spec$threshold),
                                     dplyr::n(), spec$threshold),
                     .groups = "drop")
  ggplot(ds, aes(model, value, fill = model)) +
    geom_hline(yintercept = spec$threshold, linetype = 2,
               colour = "grey30") +
    geom_boxplot(width = 0.6, linewidth = 0.3, outliers = FALSE) +
    geom_jitter(width = 0.15, height = 0, size = 0.6, alpha = 0.35) +
    scale_x_discrete(labels = stats::setNames(counts$label,
                                              as.character(counts$model))) +
    scale_fill_manual(NULL, values = cols, guide = "none") +
    coord_cartesian(ylim = lims) +
    labs(x = NULL, y = NULL, title = season_lab[[season]]) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 17),
          axis.text = element_text(size = 12),
          axis.title.y = element_text(size = 14))
}
for (s in c("train", "test")) {
  panels <- lapply(names(season_lab), panel, s = s)
  for (i in c(1, 3)) {
    panels[[i]] <- panels[[i]] + labs(y = spec$y)
  }
  g <- patchwork::wrap_plots(panels, nrow = 2) +
    plot_annotation(
      title = spec$title,
      subtitle = c(train = "Nowcasting", test = "Forecasting")[[s]],
      theme = theme(plot.title = element_text(size = 21, face = "bold"),
                    plot.subtitle = element_text(size = 18))
    )
  plot_out <- sprintf("plots/%s_%s.png", SCORE, xgb_scores_tag(s, REACH))
  ggsave(plot_out, g, width = 12, height = 11, dpi = 140)
  cat(sprintf("-> %s\n", plot_out))
}
