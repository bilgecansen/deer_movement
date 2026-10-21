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
#' The summary, printed and plotted: per deer-year and coarse class
#' (COARSE: forest, agricultural, grass, other, as Gilbertson et al. 2025
#' coarsen theirs), the CRPS of the real share of the class under the
#' model's 10 path shares (scoringRules::crps_sample()), each class on its
#' own; lower is better. A class neither the deer nor any model used in a
#' deer-year is left out. Each model is set against the GAM null: the
#' null's CRPS minus the model's, above 0 where the model does better.
#'
#' MODE "onestep" runs the same summary on the one-step draws of
#' sim_onestep_xgb.R instead of the paths: at every real step with a turning
#' angle, 10 endpoints from each model's kernel at the real start, draw r of
#' every step standing in for path r, against the real ends of those steps
#' in place of the real fixes.
#'
#' Inputs:  data/tracks/data_<key>.rds, sims_xgb_<key>.rds and
#'          sims_gam_<key>.rds in xgb_sims_dir(SET, REACH) (paths), or
#'          onestep_<key>.rds in xgb_onestep_dir(SET, REACH) (onestep), the
#'          landcover rasters
#' Output:  filters/xgb/uhc_landcover_<tag>.rds — class counts per
#'          deer-year, source (real or model) and path; and
#'          plots/uhc_landcover_crps_<tag>.png; <tag> is
#'          xgb_scores_tag(SET, REACH), with "_onestep" appended in that mode
#'
#' Configuration: edit the block below before running. The set and the mode
#' can also be given on the command line:
#'   Rscript scripts/xgb/uhc_landcover_xgb.R train onestep

# Configuration ---------------------------------------------------------------
SET <- "test"
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1) {
  SET <- .args[1]
}
stopifnot(SET %in% c("test", "train"))
# "paths": the free-running paths; "onestep": the one-step draws
MODE <- "paths"
if (length(.args) >= 2) {
  MODE <- .args[2]
}
stopifnot(MODE %in% c("paths", "onestep"))
REACH <- "gamma99"
TOLERANCE_MINS <- 30
# The coarse classes of the CRPS summary, as Gilbertson et al. (2025)
# coarsen theirs; open water, where no model steps, goes with other
COARSE <- c(forest = "forest", wetland_forested = "forest",
            corn = "agricultural", soybeans = "agricultural",
            alfalfa_hay = "agricultural", small_grains = "agricultural",
            other_ag = "agricultural", grassland = "grass",
            developed = "other", wetland_open = "other",
            open_water = "other")
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
if (MODE == "onestep") {
  tag <- paste0(tag, "_onestep")
  sims_dir <- xgb_onestep_dir(SET, REACH)
} else {
  sims_dir <- xgb_sims_dir(SET, REACH)
}
deer <- canonical_deer(readRDS("library/SW_filtered_deer.RData"))
deer <- deer[deer$set == SET, ]
keys <- sprintf("%s_%s_%d", deer$id, deer$season, as.integer(deer$year))
have <- if (MODE == "onestep") {
  file.exists(file.path(sims_dir, sprintf("onestep_%s.rds", keys)))
} else {
  file.exists(file.path(sims_dir, sprintf("sims_xgb_%s.rds", keys))) &
    file.exists(file.path(sims_dir, sprintf("sims_gam_%s.rds", keys)))
}
if (any(!have)) {
  cat(sprintf("%d deer-years lack xgb or GAM paths and are left out\n",
              sum(!have)))
}
keys <- keys[have]

count_one <- function(key) {
  source("scripts/helper_functions.R")
  track <- readRDS(sprintf("data/tracks/data_%s.rds", key))
  stp <- track$stp[[1]]
  if (MODE == "onestep") {
    os <- readRDS(file.path(sims_dir, sprintf("onestep_%s.rds", key)))
    times <- os$times
    # The ends of the real steps drawn from, each with its own draws
    real <- os$real
    tick <- seq_along(times)
    picks <- os$picks
    pick_label <- ifelse(picks$best, "xgb best", "xgb simplest")
    names(pick_label) <- picks$config
    model_paths <- c(
      stats::setNames(os$paths[picks$config], pick_label[picks$config]),
      list(`GAM null` = os$paths$null_gam,
           `Gilbertson` = os$paths$gilbertson_gam)
    )
  } else {
    xs <- readRDS(file.path(sims_dir, sprintf("sims_xgb_%s.rds", key)))
    gs <- readRDS(file.path(sims_dir, sprintf("sims_gam_%s.rds", key)))
    stopifnot(identical(xs$times, gs$times))
    times <- xs$times

    # The real fixes after the first
    last_of_burst <- !duplicated(stp$burst_, fromLast = TRUE)
    real <- dplyr::bind_rows(
      tibble::tibble(t_ = stp$t1_, x_ = stp$x1_, y_ = stp$y1_),
      tibble::tibble(t_ = stp$t2_[last_of_burst],
                     x_ = stp$x2_[last_of_burst],
                     y_ = stp$y2_[last_of_burst])
    ) |>
      dplyr::distinct(t_, .keep_all = TRUE) |>
      dplyr::arrange(t_)
    real <- real[-1, ]
    # The clock tick of each real fix
    gap <- abs(outer(as.numeric(real$t_), as.numeric(times), "-"))
    tick <- max.col(-gap, ties.method = "first")
    stopifnot(gap[cbind(seq_along(tick), tick)] <= TOLERANCE_MINS * 60)

    picks <- xs$picks
    pick_label <- ifelse(picks$best, "xgb best", "xgb simplest")
    names(pick_label) <- picks$config
    model_paths <- c(
      stats::setNames(xs$paths, pick_label[names(xs$paths)]),
      list(`GAM null` = gs$paths$null_gam,
           `Gilbertson` = gs$paths$gilbertson_gam)
    )
  }

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

# Class shares ---------------------------------------------------------------
models <- c("xgb best", "xgb simplest", "GAM null", "Gilbertson")
classes <- sort(unique(tab$class))
full <- tab |>
  tidyr::complete(tidyr::nesting(key, season, source, path), class = classes,
                  fill = list(n = 0)) |>
  dplyr::group_by(key, source, path) |>
  dplyr::mutate(share = n / sum(n)) |>
  dplyr::ungroup()

# Words for the figure: the models' positions
if (MODE == "onestep") {
  set_lab <- paste(SET, "(one step)")
  sims_word <- "draws"
} else {
  set_lab <- SET
  sims_word <- "paths"
}
cols <- c("xgb best" = "#1B9E77", "xgb simplest" = "#66A61E",
          "gam_hab" = "#D95F02")
season_lab <- c(br = "Breeding", fa = "Fawning", nb = "Non-breeding",
                pf = "Post-fawning")

# CRPS per coarse class ------------------------------------------------------
# The real share of each coarse class against the model's 10 path shares,
# each class on its own; lower is better
stopifnot(all(classes %in% names(COARSE)))
coarse_lab <- c(forest = "Forest", agricultural = "Agricultural",
                grass = "Grass", other = "Other")
coarse <- full |>
  dplyr::mutate(coarse = factor(COARSE[class], names(coarse_lab))) |>
  dplyr::group_by(key, season, source, path, coarse) |>
  dplyr::summarise(share = sum(share), .groups = "drop")
real_coarse <- coarse |>
  dplyr::filter(source == "real") |>
  dplyr::select(key, coarse, real = share)
crps <- coarse |>
  dplyr::filter(source != "real") |>
  dplyr::left_join(real_coarse, by = c("key", "coarse")) |>
  dplyr::group_by(key, season, source, coarse) |>
  dplyr::summarise(crps = scoringRules::crps_sample(real[1], dat = share),
                   any_use = real[1] > 0 | any(share > 0),
                   .groups = "drop") |>
  # A class neither the deer nor any model used in the deer-year says
  # nothing
  dplyr::group_by(key, coarse) |>
  dplyr::filter(any(any_use)) |>
  dplyr::ungroup() |>
  dplyr::mutate(source = factor(source, models))
crps_wide <- crps |>
  dplyr::select(key, season, coarse, source, crps) |>
  tidyr::pivot_wider(names_from = source, values_from = crps)
cat("\n=== CRPS of the real share per coarse class, median over deer-years",
    "(lower is better) ===\n")
print(as.data.frame(crps |>
                      dplyr::group_by(coarse, season, source) |>
                      dplyr::summarise(crps = median(crps), .groups = "drop") |>
                      tidyr::pivot_wider(names_from = source,
                                         values_from = crps)),
      digits = 3, row.names = FALSE)
cat("\n=== deer-years where a model scores better than the GAM null, per",
    "coarse class (difference: the null's CRPS minus the model's) ===\n")
for (cl in names(coarse_lab)) {
  w <- crps_wide[crps_wide$coarse == cl, ]
  for (m in c("xgb best", "xgb simplest", "Gilbertson")) {
    d <- w[["GAM null"]] - w[[m]]
    cat(sprintf(paste("%-12s %-13s better in %d of %d, tied %d (median",
                      "difference %.4f)\n"),
                cl, m, sum(d > 0), length(d), sum(d == 0), median(d)))
  }
}
d_crps <- crps_wide |>
  tidyr::pivot_longer(c("xgb best", "xgb simplest", "Gilbertson"),
                      names_to = "model", values_to = "crps") |>
  dplyr::mutate(diff = .data[["GAM null"]] - crps,
                model = factor(dplyr::recode(model, Gilbertson = "gam_hab"),
                               names(cols)))
n_better_crps <- d_crps |>
  dplyr::group_by(coarse, season, model) |>
  dplyr::summarise(label = sprintf("%d of %d", sum(diff > 0), dplyr::n()),
                   .groups = "drop")
# The 1% most extreme differences of each class are not drawn as points
d_pts <- d_crps |>
  dplyr::group_by(coarse) |>
  dplyr::filter(diff >= stats::quantile(diff, 0.01),
                diff <= stats::quantile(diff, 0.99)) |>
  dplyr::ungroup()
g5 <- ggplot(d_crps, aes(model, diff, fill = model)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey30") +
  geom_boxplot(width = 0.6, linewidth = 0.3, outliers = FALSE) +
  geom_jitter(data = d_pts, width = 0.15, height = 0, size = 0.6,
              alpha = 0.4) +
  geom_text(data = n_better_crps, aes(model, -Inf, label = label),
            inherit.aes = FALSE, vjust = -0.6, size = 3) +
  scale_fill_manual(NULL, values = cols, guide = "none") +
  scale_y_continuous(breaks = scales::breaks_width(0.05),
                     expand = expansion(mult = c(0.15, 0.05))) +
  facet_grid(coarse ~ season,
             labeller = labeller(coarse = coarse_lab,
                                 season = season_lab)) +
  labs(x = NULL,
       y = "gam_hrc's CRPS minus the model's, per deer-year",
       title = sprintf(paste("Share of each coarse landcover class per",
                             "deer-year: CRPS against gam_hrc, %s"),
                       set_lab),
       subtitle = paste("Above 0: the model's", sims_word, "put the real",
                        "deer's share of the class closer than gam_hrc's.",
                        "Numbers: deer-years above 0, of those where the",
                        "deer or a model used the class.")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(sprintf("plots/uhc_landcover_crps_%s.png", tag), g5,
       width = 14, height = 10, dpi = 140)
cat(sprintf("-> plots/uhc_landcover_crps_%s.png\n", tag))
