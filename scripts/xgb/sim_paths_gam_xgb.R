#' @description
#' Simulate every deer-year of one set from the per-deer-year GAMs the xgb
#' filters set the xgb models against: the GAM null (fit_null_gam_xgb.R)
#' and the Gilbertson et al. (2025) habitat model
#' (fit_gilbertson_gam_xgb.R), N_PATHS paths each. They are simulated as
#' sim_paths_xgb.R simulates the xgb models, so the paths can be compared:
#'   * the same clock: a step every STEP_HOURS from the deer's first fix to
#'     its last, through the gaps, stretched where the collar's schedule
#'     shifts (xgb_path_clock());
#'   * the same start, the deer's first fix, with a random heading;
#'   * the same reach (xgb_path_reach_rule());
#'   * the deer-year's own landscape, home-range centre included.
#' The GAM simulating a deer-year is the one fitted to the training
#' deer-year closest in time (xgb_closest_train_key()), with that year's
#' gamma and von Mises: the deer-year itself in training, and for a test
#' deer-year a training year of the same animal and season. The kernel is
#' xgb_gam_simulate_paths()'s, the one the log score scores
#' (xgb_gam_onestep_logscore()). Before any path is drawn, xgb_gam_check()
#' confirms it reads every column the GAMs read as the wrangle stored it for
#' the deer-year's own random points.
#'
#' Every GAM's paths are drawn from the same seed, SEED, with the
#' Mersenne-Twister generator whatever the session's, so a rerun draws the
#' same paths.
#'
#' Inputs:  data/tracks/data_<key>.rds, results/xgb/null_gam/,
#'          results/xgb/gilbertson_gam/, the rasters the wrangle read
#' Output:  sims_gam_<key>.rds in xgb_sims_dir(SET, REACH), next to the xgb
#'          paths (sims_xgb_<key>.rds): the paths of each GAM, the fits
#'          they came from and the run time
#'
#' Configuration: edit the block below before running. The set can also be
#' given on the command line: Rscript scripts/xgb/sim_paths_gam_xgb.R test
#' A second argument names the GAMs to (re)simulate, comma-separated, e.g.
#'   Rscript scripts/xgb/sim_paths_gam_xgb.R test gilbertson_gam
#' Then every deer-year of the set is simulated from those GAMs, and their
#' paths replace those in the deer-year's file; the other GAMs' paths in it
#' are kept.

# Configuration ---------------------------------------------------------------
# The deer-years to simulate: "test" or "train"
SET <- "test"
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1) {
  SET <- .args[1]
}
stopifnot(SET %in% c("test", "train"))
GAMS <- c("null_gam", "gilbertson_gam")
# TRUE: simulate GAMS for every deer-year, keeping the other GAMs' paths
UPDATE <- FALSE
if (length(.args) >= 2) {
  GAMS <- strsplit(.args[2], ",")[[1]]
  UPDATE <- TRUE
}
stopifnot(all(GAMS %in% c("null_gam", "gilbertson_gam")))
N_PATHS <- 10L
STEP_HOURS <- 4
# How far a step can reach: "longest" or "gamma99" (see sim_paths_xgb.R)
REACH <- "gamma99"
SEED <- 1L
# Deer-years to simulate: NULL for every one of the set not yet simulated
KEYS <- NULL
N_PAR <- parallel::detectCores() - 1L

# Load packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(furrr)
})

# helper functions
source("scripts/helper_functions.R")

out_dir <- xgb_sims_dir(SET, REACH)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (is.null(KEYS)) {
  deer <- canonical_deer(readRDS("library/SW_filtered_deer.RData"))
  deer <- deer[deer$set == SET, ]
  KEYS <- sprintf("%s_%s_%d", deer$id, deer$season, as.integer(deer$year))
  # The longest deer-years first
  hours <- vapply(deer$stp, function(s) {
    as.numeric(difftime(max(s$t2_), min(s$t1_), units = "hours"))
  }, numeric(1))
  KEYS <- KEYS[order(-hours)]
  if (!UPDATE) {
    KEYS <- KEYS[!file.exists(file.path(out_dir,
                                        sprintf("sims_gam_%s.rds", KEYS)))]
  }
}
cat(sprintf("simulating %d %s deer-years from %s on %d workers (%s)\n",
            length(KEYS), SET, paste(GAMS, collapse = ", "), N_PAR,
            format(Sys.time(), "%H:%M")))

sim_one <- function(key) {
  source("scripts/helper_functions.R")
  t_all <- Sys.time()
  track <- readRDS(sprintf("data/tracks/data_%s.rds", key))
  stopifnot(track$set == SET)
  season <- track$season
  stp <- track$stp[[1]]

  # The clock of sim_paths_xgb.R
  last_of_burst <- !duplicated(stp$burst_, fromLast = TRUE)
  fixes <- sort(unique(c(stp$t1_, stp$t2_[last_of_burst])))
  times <- xgb_path_clock(fixes, STEP_HOURS)
  start <- c(stp$x1_[which.min(stp$t1_)], stp$y1_[which.min(stp$t1_)])
  reach <- xgb_path_reach_rule(track, REACH)
  land <- xgb_path_landscape(track, load_landcover(track$year, season),
                             load_ndvi(track$year), load_topo())
  gam_key <- xgb_closest_train_key(track$id, season, as.integer(track$year))

  out_path <- file.path(out_dir, sprintf("sims_gam_%s.rds", key))
  paths <- list()
  seconds <- c()
  if (UPDATE && file.exists(out_path)) {
    old <- readRDS(out_path)
    stopifnot(identical(old$times, times), old$reach == reach,
              old$gam_key == gam_key)
    paths <- old$paths
    seconds <- old$seconds
  }
  for (gam in GAMS) {
    m <- readRDS(sprintf("results/xgb/%s/%s_%s.rds", gam, gam, gam_key))
    end_vars <- xgb_gam_end_vars(m$fit)
    xgb_gam_check(xgb_gam_data(track, 25L, 1L, end_vars), land, end_vars)
    t0 <- Sys.time()
    paths[[gam]] <- withr::with_seed(
      SEED,
      xgb_gam_simulate_paths(m$fit, land, start, times, reach, N_PATHS,
                             m$sl_distr, m$ta_distr),
      .rng_kind = "Mersenne-Twister",
      .rng_normal_kind = "Inversion",
      .rng_sample_kind = "Rejection"
    )
    seconds[gam] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  saveRDS(list(key = key, season = season, gam_key = gam_key, reach = reach,
               reach_rule = REACH, times = times, paths = paths,
               seconds = seconds,
               minutes = as.numeric(difftime(Sys.time(), t_all,
                                             units = "mins"))),
          out_path)
  "ok"
}

if (length(KEYS)) {
  future::plan(future::multisession, workers = min(N_PAR, length(KEYS)))
  status <- furrr::future_map_chr(
    KEYS,
    function(key) {
      tryCatch(sim_one(key),
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

files <- list.files(out_dir, pattern = "^sims_gam_.*[.]rds$",
                    full.names = TRUE)
stopped <- vapply(files, function(f) {
  x <- readRDS(f)
  sum(vapply(x$paths, function(p) length(attr(p, "stopped")), integer(1)))
}, integer(1))
cat(sprintf(paste0("-> %s: %d deer-years with GAM paths, %d paths stopped ",
                   "early (%s)\n"),
            out_dir, length(files), sum(stopped),
            format(Sys.time(), "%H:%M")))
