#' @description
#' Path-level used-habitat calibration (UHC; Fieberg et al. 2018) for
#' landcover: do the simulated paths spend their time in the same landcover
#' classes as the real deer? For every deer-year of one set it sets the
#' landcover at the real fixes against the landcover at the simulated
#' positions at the same times, for each model simulated:
#'   the season's two xgb deer picks (sim_paths_xgb.R),
#'   the GAM null and the Gilbertson et al. (2025) model
#'   (sim_paths_gam_xgb.R).
#' Every model's paths start at the deer's first fix and run free on the
#' same clock, so the first fix is left out of both. Each simulated path is
#' read at the clock time of every later real fix (xgb_path_clock() puts a
#' tick within TOLERANCE_MINS of each); a path that stopped early has no
#' position from then on, and those times are left out of that path only.
#' Landcover is read from the deer-year's own annual raster, the one the
#' models read, open water included.
#'
#' Four summaries, printed and plotted:
#'   * per season, the share of each class over all the season's real fixes
#'     against the same share over the simulated positions, path by path
#'     (path r of every deer-year pooled), so the 10 paths give 10 pooled
#'     shares and their range is the envelope;
#'   * per deer-year, the total variation distance between the real class
#'     shares and each model's mean share over its paths (half the sum of
#'     the absolute differences: 0 the same mix, 1 no class in common), and
#'     which models come closer than the GAM null;
#'   * per deer-year and class, the model's mean share over its paths minus
#'     the real share, its spread over the deer-years showing how far
#'     individual deer are missed;
#'   * per deer-year, the energy score of the real class mix under the
#'     model's paths (scoringRules::es_sample(): the vectors of class shares,
#'     the real one against the 10 simulated ones). It is the mean distance
#'     from the paths to the real mix less half the mean distance between
#'     the paths, so it rewards paths close to the real deer and a spread no
#'     wider than needed; lower is better, 0 every path exactly the real
#'     mix. Each model is set against the GAM null deer-year by deer-year.
#'
#' Inputs:  data/tracks/data_<key>.rds, sims_xgb_<key>.rds and
#'          sims_gam_<key>.rds in xgb_sims_dir(SET, REACH), the landcover
#'          rasters
#' Output:  filters/xgb/uhc_landcover_<tag>.rds — class counts per
#'          deer-year, source (real or model) and path; and
#'          plots/uhc_landcover_<tag>.png, plots/uhc_landcover_tvd_<tag>.png,
#'          plots/uhc_landcover_deer_error_<tag>.png and
#'          plots/uhc_landcover_energy_<tag>.png; <tag> is
#'          xgb_scores_tag(SET, REACH)
#'
#' Configuration: edit the block below before running. The set can also be
#' given on the command line: Rscript scripts/xgb/uhc_landcover_xgb.R train

# Configuration ---------------------------------------------------------------
SET <- "test"
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1) {
  SET <- .args[1]
}
stopifnot(SET %in% c("test", "train"))
REACH <- "gamma99"
TOLERANCE_MINS <- 30
# FALSE: reuse the saved class counts if they exist, and only summarise
RECOUNT <- FALSE
N_PAR <- parallel::detectCores() - 1L

# Load packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(furrr)
})

# helper functions
source("scripts/helper_functions.R")

tag <- xgb_scores_tag(SET, REACH)
sims_dir <- xgb_sims_dir(SET, REACH)
deer <- canonical_deer(readRDS("library/SW_filtered_deer.RData"))
deer <- deer[deer$set == SET, ]
keys <- sprintf("%s_%s_%d", deer$id, deer$season, as.integer(deer$year))
have <- file.exists(file.path(sims_dir, sprintf("sims_xgb_%s.rds", keys))) &
  file.exists(file.path(sims_dir, sprintf("sims_gam_%s.rds", keys)))
if (any(!have)) {
  cat(sprintf("%d deer-years lack xgb or GAM paths and are left out\n",
              sum(!have)))
}
keys <- keys[have]

count_one <- function(key) {
  source("scripts/helper_functions.R")
  track <- readRDS(sprintf("data/tracks/data_%s.rds", key))
  stp <- track$stp[[1]]
  xs <- readRDS(file.path(sims_dir, sprintf("sims_xgb_%s.rds", key)))
  gs <- readRDS(file.path(sims_dir, sprintf("sims_gam_%s.rds", key)))
  stopifnot(identical(xs$times, gs$times))
  times <- xs$times

  # The real fixes after the first
  last_of_burst <- !duplicated(stp$burst_, fromLast = TRUE)
  real <- dplyr::bind_rows(
    tibble::tibble(t_ = stp$t1_, x_ = stp$x1_, y_ = stp$y1_),
    tibble::tibble(t_ = stp$t2_[last_of_burst], x_ = stp$x2_[last_of_burst],
                   y_ = stp$y2_[last_of_burst])
  ) |>
    dplyr::distinct(t_, .keep_all = TRUE) |>
    dplyr::arrange(t_)
  real <- real[-1, ]
  # The clock tick of each real fix
  gap <- abs(outer(as.numeric(real$t_), as.numeric(times), "-"))
  tick <- max.col(-gap, ties.method = "first")
  stopifnot(gap[cbind(seq_along(tick), tick)] <= TOLERANCE_MINS * 60)

  lc <- load_landcover(track$year, track$season)[["wiscland"]]
  class_at <- function(x, y) {
    v <- terra::extract(lc, cbind(x, y))[[1]]
    as.character(v)
  }
  counts <- function(cls) {
    tibble::tibble(class = cls[!is.na(cls)]) |> dplyr::count(class)
  }

  real_n <- counts(class_at(real$x_, real$y_)) |>
    dplyr::mutate(source = "real", path = NA_integer_)
  picks <- xs$picks
  pick_label <- ifelse(picks$best, "xgb best", "xgb simplest")
  names(pick_label) <- picks$config
  model_paths <- c(
    stats::setNames(xs$paths, pick_label[names(xs$paths)]),
    list(`GAM null` = gs$paths$null_gam,
         `Gilbertson` = gs$paths$gilbertson_gam)
  )
  sim_n <- purrr::imap_dfr(model_paths, function(p, model) {
    purrr::map_dfr(sort(unique(p$nsim)), function(r) {
      pr <- p[p$nsim == r, ]
      stopifnot(identical(pr$t_, times))
      pos <- pr[tick, ]
      counts(class_at(pos$x_, pos$y_)) |>
        dplyr::mutate(source = model, path = r)
    })
  })
  dplyr::bind_rows(real_n, sim_n) |>
    dplyr::mutate(key = key, season = track$season, n_fixes = nrow(real),
                  config_best = picks$config[picks$best],
                  config_simplest = picks$config[picks$simplest])
}

out_path <- sprintf("filters/xgb/uhc_landcover_%s.rds", tag)
if (RECOUNT || !file.exists(out_path)) {
  cat(sprintf("reading landcover for %d %s deer-years (%s)\n", length(keys),
              tag, format(Sys.time(), "%H:%M")))
  future::plan(future::multisession, workers = min(N_PAR, length(keys)))
  tab <- furrr::future_map_dfr(keys, count_one,
                               .options = furrr::furrr_options(
                                 seed = TRUE, chunk_size = 1L
                               ))
  future::plan(future::sequential)
  saveRDS(tab, out_path)
  cat(sprintf("-> %s (%s)\n", out_path, format(Sys.time(), "%H:%M")))
} else {
  tab <- readRDS(out_path)
  stopifnot(setequal(unique(tab$key), keys))
  cat(sprintf("reusing %s\n", out_path))
}

# Summaries ------------------------------------------------------------------
models <- c("xgb best", "xgb simplest", "GAM null", "Gilbertson")
classes <- sort(unique(tab$class))
full <- tab |>
  tidyr::complete(tidyr::nesting(key, season, source, path), class = classes,
                  fill = list(n = 0)) |>
  dplyr::group_by(key, source, path) |>
  dplyr::mutate(share = n / sum(n)) |>
  dplyr::ungroup()

# Per season: real pooled over deer-years, each path r pooled over them
season_share <- full |>
  dplyr::group_by(season, source, path, class) |>
  dplyr::summarise(n = sum(n), .groups = "drop") |>
  dplyr::group_by(season, source, path) |>
  dplyr::mutate(share = n / sum(n)) |>
  dplyr::ungroup()
env <- season_share |>
  dplyr::filter(source != "real") |>
  dplyr::group_by(season, source, class) |>
  dplyr::summarise(lo = min(share), hi = max(share), mean = mean(share),
                   .groups = "drop") |>
  dplyr::left_join(season_share |>
                     dplyr::filter(source == "real") |>
                     dplyr::select(season, class, real = share),
                   by = c("season", "class")) |>
  dplyr::mutate(inside = real >= lo & real <= hi,
                source = factor(source, models))
cat("\n=== real share inside the 10-path range, per season and model",
    "(classes with >= 1% of real fixes) ===\n")
print(as.data.frame(env |>
                      dplyr::filter(real >= 0.01) |>
                      dplyr::group_by(season, source) |>
                      dplyr::summarise(classes = n(), inside = sum(inside),
                                       .groups = "drop") |>
                      tidyr::pivot_wider(names_from = source,
                                         values_from = inside)),
      row.names = FALSE)

# Per deer-year: total variation distance to the real mix
tvd <- full |>
  dplyr::filter(source != "real") |>
  dplyr::group_by(key, season, source, class) |>
  dplyr::summarise(share = mean(share), .groups = "drop") |>
  dplyr::left_join(full |>
                     dplyr::filter(source == "real") |>
                     dplyr::select(key, class, real = share),
                   by = c("key", "class")) |>
  dplyr::group_by(key, season, source) |>
  dplyr::summarise(tvd = sum(abs(share - real)) / 2, .groups = "drop") |>
  dplyr::mutate(source = factor(source, models))
cat("\n=== total variation distance to the real landcover mix, median over",
    "deer-years ===\n")
print(as.data.frame(tvd |>
                      dplyr::group_by(season, source) |>
                      dplyr::summarise(tvd = median(tvd), .groups = "drop") |>
                      tidyr::pivot_wider(names_from = source,
                                         values_from = tvd)),
      digits = 3, row.names = FALSE)
wide <- tvd |> tidyr::pivot_wider(names_from = source, values_from = tvd)
cat("\n=== deer-years where a model comes closer than the GAM null ===\n")
for (m in c("xgb best", "xgb simplest", "Gilbertson")) {
  d <- wide[[m]] - wide[["GAM null"]]
  cat(sprintf("%-13s closer in %d of %d (median difference %.3f)\n", m,
              sum(d < 0), length(d), median(d)))
}
d <- wide[["xgb best"]] - wide[["Gilbertson"]]
cat(sprintf("%-13s closer than Gilbertson in %d of %d\n", "xgb best",
            sum(d < 0), length(d)))

# Plots ----------------------------------------------------------------------
cols <- c("xgb best" = "#1B9E77", "xgb simplest" = "#66A61E",
          "GAM null" = "#7570B3", "Gilbertson" = "#D95F02")
lab_class <- function(x) gsub("_", " ", x)
keep_cls <- env |>
  dplyr::group_by(class) |>
  dplyr::summarise(m = max(real, mean), .groups = "drop") |>
  dplyr::filter(m >= 0.005) |>
  dplyr::pull(class)
g1 <- ggplot(dplyr::filter(env, class %in% keep_cls),
             aes(x = lab_class(class))) +
  geom_linerange(aes(ymin = lo, ymax = hi, colour = source),
                 position = position_dodge(width = 0.7), linewidth = 2.2) +
  geom_point(aes(y = real), shape = 21, fill = "black", colour = "white",
             size = 2.4) +
  scale_colour_manual("Simulated paths\n(range over 10 paths)",
                      values = cols) +
  facet_wrap(~season, ncol = 2,
             labeller = labeller(season = c(br = "Breeding",
                                            fa = "Fawning",
                                            nb = "Non-breeding",
                                            pf = "Post-fawning"))) +
  labs(x = NULL, y = "Share of fixes in the class",
       title = sprintf("Landcover at real fixes (black) vs simulated paths, %s",
                       SET),
       caption = paste("Each bar: the range over the 10 simulated paths of",
                       "the share of positions in the class, every deer-year",
                       "of the season pooled.")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(sprintf("plots/uhc_landcover_%s.png", tag), g1, width = 13,
       height = 8, dpi = 140)

lims <- range(c(wide[["GAM null"]], wide[["xgb best"]], wide[["Gilbertson"]]))
pd <- dplyr::bind_rows(
  tibble::tibble(season = wide$season, null = wide[["GAM null"]],
                 other = wide[["xgb best"]], model = "xgb best"),
  tibble::tibble(season = wide$season, null = wide[["GAM null"]],
                 other = wide[["Gilbertson"]], model = "Gilbertson")
)
g2 <- ggplot(pd, aes(null, other)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
  geom_point(aes(colour = model), alpha = 0.6, size = 1.6) +
  scale_colour_manual(NULL, values = cols) +
  facet_grid(model ~ season) +
  coord_equal(xlim = lims, ylim = lims) +
  labs(x = "GAM null: distance to the real landcover mix",
       y = "Model: distance to the real landcover mix",
       title = sprintf(paste("Per deer-year landcover mismatch, %s (below the",
                             "line: closer to the real deer than the null)"),
                       SET)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")
ggsave(sprintf("plots/uhc_landcover_tvd_%s.png", tag), g2, width = 13,
       height = 7, dpi = 140)
cat(sprintf("-> plots/uhc_landcover_%s.png, plots/uhc_landcover_tvd_%s.png\n",
            tag, tag))

# Per deer-year, by class ----------------------------------------------------
# Each model's mean share over its paths, against the real share
per_deer <- full |>
  dplyr::filter(source != "real") |>
  dplyr::group_by(key, season, source, class) |>
  dplyr::summarise(mean = mean(share), hi = max(share), .groups = "drop") |>
  dplyr::left_join(full |>
                     dplyr::filter(source == "real") |>
                     dplyr::select(key, class, real = share),
                   by = c("key", "class")) |>
  # A class neither the deer nor any path of the model used says nothing
  dplyr::filter(real > 0 | hi > 0) |>
  dplyr::mutate(error = mean - real, source = factor(source, models))
cat("\n=== per deer-year and class: median absolute error ===\n")
print(as.data.frame(per_deer |>
                      dplyr::group_by(source) |>
                      dplyr::summarise(abs_error = median(abs(error)),
                                       .groups = "drop")),
      digits = 3, row.names = FALSE)

season_lab <- c(br = "Breeding", fa = "Fawning", nb = "Non-breeding",
                pf = "Post-fawning")
pd_deer <- dplyr::filter(per_deer, class %in% keep_cls) |>
  dplyr::mutate(class = factor(lab_class(class),
                               lab_class(keep_cls[order(-tapply(
                                 env$real, env$class, mean
                               )[keep_cls])])))
ylim_err <- c(-0.4, 0.4)
n_out <- pd_deer |>
  dplyr::filter(error < ylim_err[1] | error > ylim_err[2]) |>
  dplyr::count(source)
g3 <- ggplot(pd_deer, aes(class, error, fill = source)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey30") +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7,
               linewidth = 0.3, outlier.size = 0.6, outlier.alpha = 0.5) +
  scale_fill_manual("Simulated paths", values = cols) +
  facet_wrap(~season, ncol = 2, labeller = labeller(season = season_lab)) +
  coord_cartesian(ylim = ylim_err) +
  labs(x = NULL,
       y = "Simulated share minus real share, per deer-year",
       title = sprintf(paste("Landcover per deer-year: how far each model's",
                             "paths are from the real deer, %s"), SET),
       caption = paste0(
         "One value per deer-year and class: the model's mean share over ",
         "its 10 paths minus the real deer's share. 0 = the same share. ",
         "Values beyond ±0.4 not drawn: ",
         paste(sprintf("%s %d", n_out$source, n_out$n), collapse = ", "),
         "."
       )) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(sprintf("plots/uhc_landcover_deer_error_%s.png", tag), g3,
       width = 14, height = 8.5, dpi = 140)

# Energy score per deer-year -------------------------------------------------
# The real class mix against the model's 10 path mixes; lower is better
real_mix <- full |>
  dplyr::filter(source == "real") |>
  dplyr::select(key, class, real = share)
es <- full |>
  dplyr::filter(source != "real") |>
  dplyr::left_join(real_mix, by = c("key", "class")) |>
  dplyr::group_by(key, season, source) |>
  dplyr::group_modify(function(d, g) {
    # Classes in rows, paths in columns
    m <- tapply(d$share, list(d$class, d$path), sum)
    y <- as.numeric(tapply(d$real, d$class, `[`, 1)[rownames(m)])
    stopifnot(!anyNA(m), !anyNA(y))
    tibble::tibble(es = scoringRules::es_sample(y, dat = m))
  }) |>
  dplyr::ungroup() |>
  dplyr::mutate(source = factor(source, models))
es_wide <- es |> tidyr::pivot_wider(names_from = source, values_from = es)
cat("\n=== energy score of the real landcover mix, median over deer-years",
    "(lower is better) ===\n")
print(as.data.frame(es |>
                      dplyr::group_by(season, source) |>
                      dplyr::summarise(es = median(es), .groups = "drop") |>
                      tidyr::pivot_wider(names_from = source,
                                         values_from = es)),
      digits = 3, row.names = FALSE)
cat("\n=== deer-years where a model scores better than the GAM null ===\n")
for (m in c("xgb best", "xgb simplest", "Gilbertson")) {
  d <- es_wide[[m]] - es_wide[["GAM null"]]
  cat(sprintf("%-13s better in %d of %d (median difference %.4f)\n", m,
              sum(d < 0), length(d), median(d)))
}
d_es <- es_wide |>
  tidyr::pivot_longer(c("xgb best", "xgb simplest", "Gilbertson"),
                      names_to = "model", values_to = "es") |>
  dplyr::mutate(diff = es - .data[["GAM null"]],
                model = factor(model, models))
n_better <- d_es |>
  dplyr::group_by(season, model) |>
  dplyr::summarise(label = sprintf("%d of %d", sum(diff < 0), dplyr::n()),
                   .groups = "drop")
ylim_es <- stats::quantile(d_es$diff, c(0.01, 0.99))
g4 <- ggplot(d_es, aes(model, diff, fill = model)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey30") +
  geom_boxplot(width = 0.6, linewidth = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, size = 0.8, alpha = 0.45) +
  geom_text(data = n_better, aes(model, ylim_es[1], label = label),
            inherit.aes = FALSE, vjust = 0, size = 3.3) +
  scale_fill_manual(NULL, values = cols, guide = "none") +
  facet_wrap(~season, ncol = 4, labeller = labeller(season = season_lab)) +
  coord_cartesian(ylim = ylim_es) +
  labs(x = NULL,
       y = "Energy score minus the GAM null's, per deer-year",
       title = sprintf(paste("Landcover mix per deer-year: energy score",
                             "against the GAM null, %s"), SET),
       subtitle = paste("Below 0: the model's paths match the real deer's",
                        "landcover mix better than the null's. Numbers:",
                        "deer-years below 0."),
       caption = paste("Energy score of the real deer's class shares under",
                       "the model's 10 simulated paths; the 1% most extreme",
                       "differences fall outside the axis.")) +
  theme_bw(base_size = 11)
ggsave(sprintf("plots/uhc_landcover_energy_%s.png", tag), g4,
       width = 14, height = 6, dpi = 140)
cat(sprintf(paste0("-> plots/uhc_landcover_deer_error_%s.png, ",
                   "plots/uhc_landcover_energy_%s.png\n"), tag, tag))
