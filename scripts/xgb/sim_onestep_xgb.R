#' @description
#' One-step draws for the one-step used-habitat calibration: at every real
#' step of every deer-year of one set, N_DRAWS endpoints from each model's
#' kernel built at the real start, for
#'   the season's two final xgb deer models (xgb_onestep_draw()),
#'   the GAM null and the Gilbertson et al. (2025) model
#'   (xgb_gam_onestep_draw()).
#' The kernels are the ones the log scores score (logscore_xgb.R,
#' logscore_gam_xgb.R): the same reach rule, candidate cells, heading from
#' the step before, hour and NDVI layers. The GAM drawing for a deer-year is
#' the one fitted to the training deer-year closest in time
#' (xgb_closest_train_key()), with that year's gamma and von Mises. Before
#' any draw, xgb_path_check() and xgb_gam_check() confirm the kernels read
#' the deer-year's own steps as the models were fitted on them.
#'
#' Every step with a turning angle is drawn from, whether or not the deer's
#' own step ended within the reach, so the real use is every such step's
#' end. Draw r of every step makes up "path" r: the UHC scripts
#' (uhc_landcover_xgb.R, uhc_ndvi_xgb.R, uhc_forest_edge_xgb.R, mode
#' "onestep") read the draws as they read the free-running paths, at the
#' real steps' ends instead of the real fixes.
#'
#' Every model's draws start from the same seed, SEED, with the
#' Mersenne-Twister generator whatever the session's, so a rerun draws the
#' same endpoints.
#'
#' Inputs:  data/tracks/data_<key>.rds, results/xgb/final_selection.rds,
#'          results/xgb/final_<season>_<type>.rds, results/xgb/null_gam/,
#'          results/xgb/gilbertson_gam/, the rasters the wrangle read
#' Output:  onestep_<key>.rds in xgb_onestep_dir(SET, REACH), one per
#'          deer-year: the real steps' ends (`real`) and their times
#'          (`times`), the draws of each model (`paths`, one tibble per
#'          model as the path files hold them), the picks, the reach, and
#'          per model the log probability of every real end (`logp_end`)
#'
#' Configuration: edit the block below before running. The set can also be
#' given on the command line: Rscript scripts/xgb/sim_onestep_xgb.R test

# Configuration ---------------------------------------------------------------
# The deer-years to draw for: "test" or "train"
SET <- "test"
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1) {
  SET <- .args[1]
}
stopifnot(SET %in% c("test", "train"))
GAMS <- c("null_gam", "gilbertson_gam")
N_DRAWS <- 10L
# How far a step can reach: "longest" or "gamma99" (see sim_paths_xgb.R)
REACH <- "gamma99"
SEED <- 1L
# Deer-years to draw for: NULL for every one of the set not yet drawn
KEYS <- NULL
N_PAR <- parallel::detectCores() - 1L

# Load packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(furrr)
})

# helper functions
source("scripts/helper_functions.R")

out_dir <- xgb_onestep_dir(SET, REACH)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (is.null(KEYS)) {
  deer <- canonical_deer(readRDS("library/SW_filtered_deer.RData"))
  deer <- deer[deer$set == SET, ]
  KEYS <- sprintf("%s_%s_%d", deer$id, deer$season, as.integer(deer$year))
  # The longest deer-years first
  n_steps <- vapply(deer$stp, nrow, integer(1))
  KEYS <- KEYS[order(-n_steps)]
  KEYS <- KEYS[!file.exists(file.path(out_dir,
                                      sprintf("onestep_%s.rds", KEYS)))]
}
cat(sprintf("drawing for %d %s deer-years on %d workers (%s)\n",
            length(KEYS), SET, N_PAR, format(Sys.time(), "%H:%M")))

with_draw_seed <- function(expr) {
  withr::with_seed(SEED, expr, .rng_kind = "Mersenne-Twister",
                   .rng_normal_kind = "Inversion",
                   .rng_sample_kind = "Rejection")
}

draw_one <- function(key) {
  source("scripts/helper_functions.R")
  t_all <- Sys.time()
  track <- readRDS(sprintf("data/tracks/data_%s.rds", key))
  stopifnot(track$set == SET)
  season <- track$season
  stp <- track$stp[[1]]
  reach <- xgb_path_reach_rule(track, REACH)
  land <- xgb_path_landscape(track, load_landcover(track$year, season),
                             load_ndvi(track$year), load_topo())

  # The season's two xgb picks, read and checked as logscore_xgb.R does
  picks <- readRDS("results/xgb/final_selection.rds")
  picks <- picks[picks$season == season, ]
  models <- lapply(sprintf("results/xgb/final_%s_%s.rds", season,
                           picks$config),
                   xgb_path_model, nthread = 1L)
  names(models) <- picks$config
  deer_ids <- lapply(models, function(m) {
    i <- match(track$id, m$deer_levels)
    stopifnot(!is.na(i))
    i - 1L
  })
  v <- track$stp.var.nonp[[1]]
  layers <- stats::na.omit(unique(xgb_path_ndvi_layer(
    c(stp$t1_, stp$t2_, v$t2_), land$ndvi_times
  )))
  statics <- Map(function(m, id) xgb_path_static(m, land, layers, id),
                 models, deer_ids)
  xgb_path_check(track, land, models, statics, deer_ids)

  paths <- list()
  logp_end <- list()
  seconds <- c()
  for (cf in names(models)) {
    t0 <- Sys.time()
    d <- with_draw_seed(xgb_onestep_draw(models[[cf]], land, statics[[cf]],
                                         stp, reach, deer_ids[[cf]],
                                         N_DRAWS))
    logp_end[[cf]] <- attr(d, "logp_end")
    paths[[cf]] <- d
    seconds[cf] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }

  gam_key <- xgb_closest_train_key(track$id, season, as.integer(track$year))
  for (gam in GAMS) {
    m <- readRDS(sprintf("results/xgb/%s/%s_%s.rds", gam, gam, gam_key))
    end_vars <- xgb_gam_end_vars(m$fit)
    xgb_gam_check(xgb_gam_data(track, 25L, 1L, end_vars), land, end_vars)
    t0 <- Sys.time()
    d <- with_draw_seed(xgb_gam_onestep_draw(m$fit, land, stp, reach,
                                             m$sl_distr, m$ta_distr,
                                             N_DRAWS))
    logp_end[[gam]] <- attr(d, "logp_end")
    paths[[gam]] <- d
    seconds[gam] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }

  # The real use: the end of every step drawn from
  steps <- which(!is.na(stp$ta_))
  real <- tibble::tibble(step = steps, t_ = stp$t2_[steps],
                         x_ = stp$x2_[steps], y_ = stp$y2_[steps])
  for (nm in names(paths)) {
    stopifnot(identical(paths[[nm]]$step[paths[[nm]]$nsim == 1], steps))
  }
  saveRDS(list(key = key, season = season, set = SET, picks = picks,
               gam_key = gam_key, reach = reach, reach_rule = REACH,
               n_draws = N_DRAWS, real = real, times = real$t_,
               paths = paths, logp_end = logp_end, seconds = seconds,
               minutes = as.numeric(difftime(Sys.time(), t_all,
                                             units = "mins"))),
          file.path(out_dir, sprintf("onestep_%s.rds", key)))
  "ok"
}

if (length(KEYS)) {
  future::plan(future::multisession, workers = min(N_PAR, length(KEYS)))
  status <- furrr::future_map_chr(
    KEYS,
    function(key) {
      tryCatch(draw_one(key),
               error = function(e) paste("failed:", conditionMessage(e)))
    },
    .options = furrr::furrr_options(seed = TRUE, stdout = FALSE,
                                    chunk_size = 1L)
  )
  future::plan(future::sequential)
  failed <- status != "ok"
  if (any(failed)) {
    cat(sprintf("%d deer-years failed:\n", sum(failed)))
    cat(sprintf("  %s: %s\n", KEYS[failed], status[failed]), sep = "")
  }
}

files <- list.files(out_dir, pattern = "^onestep_.*[.]rds$",
                    full.names = TRUE)
no_draw <- vapply(files, function(f) {
  x <- readRDS(f)
  sum(vapply(x$paths, function(p) sum(is.na(p$x_)), integer(1)))
}, integer(1))
cat(sprintf(paste0("-> %s: %d deer-years with one-step draws, %d draws ",
                   "with no cell to go to (%s)\n"),
            out_dir, length(files), sum(no_draw),
            format(Sys.time(), "%H:%M")))
