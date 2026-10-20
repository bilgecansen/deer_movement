#' @description
#' Box plots of the one-step log score against the GAM null, gam_hrc
#' (movement and the pull of the home-range centre), per deer-year and per
#' step: the season's two final xgb deer models and the Gilbertson et al.
#' (2025) habitat model, gam_hab, for the training years and the test year.
#'
#' A deer-year's difference is a model's total log score over the
#' deer-year's real steps minus gam_hrc's, both summed over the steps all
#' four models could score (status "ok" for every one), so the three
#' differences share one step set. The figure shows it per step as
#' exp(difference / steps) - 1 in percent: how much more probability, per
#' step on average (the geometric mean of the ratio), the model gives the
#' cell the deer moved to than gam_hrc does. Above 0: the model gives the
#' deer's real steps more probability than gam_hrc. The per-step scores are
#' the ones logscore_xgb.R and logscore_gam_xgb.R saved.
#'
#' Inputs:  per deer-year, filters/xgb/logscore_<tag>/logscore_xgb_<key>.rds,
#'          filters/xgb/logscore_null_gam_<tag>/logscore_null_gam_<key>.rds
#'          and filters/xgb/logscore_gilbertson_gam_<tag>/
#'          logscore_gilbertson_gam_<key>.rds; results/xgb/final_selection.rds;
#'          <tag> = xgb_scores_tag(set, REACH)
#' Output:  filters/xgb/logscore_vs_gam_hrc_<REACH>.rds, the total
#'          differences and step counts per set, deer-year and model;
#'          plots/logscore_vs_gam_hrc_<tag>.png, one per set, a panel per
#'          season; and a summary printed
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
REACH <- "gamma99"
# The share of the most extreme values of each panel left outside the axis
TRIM <- 0.01
# Deer-years left out of the figures (kept in the saved table)
EXCLUDE <- c("7745_nb_2021")

# Load packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# helper functions
source("scripts/helper_functions.R")

models <- c("xgb best", "xgb simplest", "gam_hab")
sel <- readRDS("results/xgb/final_selection.rds")

# Each model's total over the steps all four could score, less gam_hrc's
read_set <- function(set) {
  tag <- xgb_scores_tag(set, REACH)
  dir_of <- function(what) sprintf("filters/xgb/logscore_%s%s", what, tag)
  keys <- sub("^logscore_xgb_(.*)[.]rds$", "\\1",
              list.files(dir_of(""), pattern = "^logscore_xgb_.*[.]rds$"))
  purrr::map_dfr(keys, function(key) {
    x <- readRDS(file.path(dir_of(""), sprintf("logscore_xgb_%s.rds", key)))
    hrc <- readRDS(file.path(dir_of("null_gam_"),
                             sprintf("logscore_null_gam_%s.rds", key)))
    hab <- readRDS(file.path(dir_of("gilbertson_gam_"),
                             sprintf("logscore_gilbertson_gam_%s.rds", key)))
    picks <- sel[sel$season == x$season, ]
    label <- c(stats::setNames(ifelse(picks$best, "xgb best", "xgb simplest"),
                               picks$config),
               null_gam = "gam_hrc", gilbertson_gam = "gam_hab")
    steps <- dplyr::bind_rows(x$steps[x$steps$config %in% picks$config, ],
                              hrc$steps, hab$steps) |>
      dplyr::mutate(model = label[config])
    stopifnot(!anyNA(steps$model))
    common <- steps |>
      dplyr::group_by(step) |>
      dplyr::summarise(all_ok = dplyr::n() == 4 && all(status == "ok"),
                       .groups = "drop") |>
      dplyr::filter(all_ok) |>
      dplyr::pull(step)
    totals <- steps |>
      dplyr::filter(step %in% common) |>
      dplyr::group_by(model) |>
      dplyr::summarise(logp = sum(logp), .groups = "drop")
    hrc_total <- totals$logp[totals$model == "gam_hrc"]
    totals |>
      dplyr::filter(model != "gam_hrc") |>
      dplyr::transmute(set = set, key = key, season = x$season, model,
                       diff = logp - hrc_total, n_steps = length(common))
  })
}

d <- dplyr::bind_rows(read_set("train"), read_set("test")) |>
  dplyr::mutate(model = factor(model, models),
                set = factor(set, c("train", "test")),
                per_step = 100 * (exp(diff / n_steps) - 1))
out <- sprintf("filters/xgb/logscore_vs_gam_hrc_%s.rds", REACH)
saveRDS(dplyr::select(d, -per_step), out)

cat("=== per-step probability over gam_hrc (%), per deer-year ===\n")
print(as.data.frame(d |>
                      dplyr::group_by(set, season, model) |>
                      dplyr::summarise(
                        deer_years = dplyr::n(),
                        median = stats::median(per_step),
                        above_0 = sum(per_step > 0),
                        median_steps = stats::median(n_steps),
                        .groups = "drop"
                      )),
      digits = 3, row.names = FALSE)

# Plots ----------------------------------------------------------------------
# One figure per set, one panel per season, each on its own axis
cols <- c("xgb best" = "#1B9E77", "xgb simplest" = "#66A61E",
          "gam_hab" = "#D95F02")
season_lab <- c(br = "Breeding", fa = "Fawning", nb = "Non-breeding",
                pf = "Post-fawning")
panel <- function(s, season) {
  ds <- d[d$set == s & d$season == season & !d$key %in% EXCLUDE, ]
  lims <- range(c(0, stats::quantile(ds$per_step, c(TRIM, 1 - TRIM))))
  # Each model's name, with its deer-years above 0 under it
  counts <- ds |>
    dplyr::group_by(model) |>
    dplyr::summarise(label = sprintf("%s\n%d of %d above 0", model[1],
                                     sum(per_step > 0), dplyr::n()),
                     .groups = "drop")
  ggplot(ds, aes(model, per_step, fill = model)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey30") +
    geom_boxplot(width = 0.6, linewidth = 0.3, outliers = FALSE) +
    geom_jitter(width = 0.15, height = 0, size = 0.6, alpha = 0.35) +
    scale_y_continuous(labels = function(b) paste0(b, "%")) +
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
    panels[[i]] <- panels[[i]] +
      labs(y = "Probability of average step\nover gam_hrc")
  }
  g <- patchwork::wrap_plots(panels, nrow = 2) +
    plot_annotation(
      title = paste("Does habitat make steps more likely than movement",
                    "and home range center alone?"),
      subtitle = c(train = "Nowcasting", test = "Forecasting")[[s]],
      theme = theme(plot.title = element_text(size = 21, face = "bold"),
                    plot.subtitle = element_text(size = 18))
    )
  plot_out <- sprintf("plots/logscore_vs_gam_hrc_%s.png",
                      xgb_scores_tag(s, REACH))
  ggsave(plot_out, g, width = 12, height = 11, dpi = 140)
  cat(sprintf("-> %s\n", plot_out))
}
cat(sprintf("-> %s\n", out))
