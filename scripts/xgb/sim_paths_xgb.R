#' @description
#' Simulate one deer-year's path from each of its season's final models
#' (results/xgb/final_selection.rds: the best and the simplest type,
#' without and with deer), N_PATHS paths per model. DEER_ONLY keeps to the
#' two deer types. The deer-year may be a test year or a training year;
#' either is simulated the same way, from its own HR centre, rasters and
#' dates.
#'
#' A path starts at the deer's first fix with a random heading and steps
#' every STEP_HOURS to its last fix, through the gaps in the fixes, so it is
#' never put back on the real deer. Where the collar's schedule shifts, the
#' step across the gap before the shifted fix is stretched to end at it
#' (xgb_path_clock()), so every real fix has a simulated position at its
#' time. A step reaches at most as far as the
#' animal's longest step with a turning angle in its training years of the
#' season (xgb_path_reach()). The kernel is xgb_simulate_paths()'s.
#'
#' Before any path is drawn, xgb_path_check() confirms the kernel reads the
#' deer's own observed and random steps as the models were fitted on them.
#'
#' Every model's paths start from the same random seed.
#'
#' Inputs:  data/tracks/data_<key>.rds, results/xgb/final_selection.rds,
#'          results/xgb/final_<season>_<type>.rds, the rasters the wrangle
#'          read (load_landcover(), load_ndvi(), load_topo(), data/HR/)
#' Output:  sims/xgb_paths/sims_xgb_<key>.rds for a test year, and
#'          sims/xgb_paths_train/sims_xgb_<key>.rds for a training year —
#'          the paths of each model, with the share of each kernel off the
#'          window and the run time
#'
#' Configuration: edit the block below before running, or let
#' run_sims_xgb.R set KEY and N_THREAD.

# Configuration ---------------------------------------------------------------
KEY <- "5227_nb_2021"
N_PATHS <- 10L
# TRUE: simulate only the deer types
DEER_ONLY <- TRUE
STEP_HOURS <- 4
SEED <- 1L
N_THREAD <- 1L

# Load packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(xgboost)
  library(terra)
  library(tidyverse)
})

# helper functions
source("scripts/helper_functions.R")

t_all <- Sys.time()
track <- readRDS(sprintf("data/tracks/data_%s.rds", KEY))
stopifnot(track$set %in% c("train", "test"))
id <- track$id
season <- track$season
stp <- track$stp[[1]]

picks <- readRDS("results/xgb/final_selection.rds")
picks <- picks[picks$season == season & (picks$deer | !DEER_ONLY), ]
models <- lapply(sprintf("results/xgb/final_%s_%s.rds", season, picks$config),
                 xgb_path_model, nthread = N_THREAD)
names(models) <- picks$config
deer_ids <- lapply(models, function(m) {
  if (is.null(m$deer_levels)) {
    return(NA_integer_)
  }
  i <- match(id, m$deer_levels)
  stopifnot(!is.na(i))
  i - 1L
})

# The 4-hour clock from the first fix to the last, through every gap, and
# following the collar where its schedule shifts (xgb_path_clock())
last_of_burst <- !duplicated(stp$burst_, fromLast = TRUE)
fixes <- sort(unique(c(stp$t1_, stp$t2_[last_of_burst])))
times <- xgb_path_clock(fixes, STEP_HOURS)
first <- times[1]
n_steps <- length(times) - 1
reach <- xgb_path_reach(id, season)
cat(sprintf("%s: %d steps from %s, reach %.0f m, models %s\n", KEY, n_steps,
            format(first, "%Y-%m-%d %H:%M"), reach,
            paste(picks$config, collapse = ", ")))

land <- xgb_path_landscape(track, load_landcover(track$year, season),
                           load_ndvi(track$year), load_topo())
v <- track$stp.var.nonp[[1]]
layers <- stats::na.omit(unique(c(
  xgb_path_ndvi_layer(c(times, v$t2_), land$ndvi_times)
)))
statics <- Map(function(m, deer_id) {
  xgb_path_static(m, land, layers, deer_id)
}, models, deer_ids)

gaps <- xgb_path_check(track, land, models, statics, deer_ids)
cat(sprintf("check passed: %d endpoints, largest score difference %.2g\n",
            nrow(v), max(gaps)))

paths <- list()
seconds <- c()
for (nm in names(models)) {
  t0 <- Sys.time()
  set.seed(SEED)
  paths[[nm]] <- xgb_simulate_paths(
    models[[nm]], land, statics[[nm]],
    start = c(stp$x1_[which.min(stp$t1_)], stp$y1_[which.min(stp$t1_)]),
    times = times, reach = reach, n_paths = N_PATHS,
    deer_id = deer_ids[[nm]]
  )
  seconds[nm] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ow <- attr(paths[[nm]], "off_window")
  cat(sprintf(
    paste0("  %-14s %5.0f s | kernel cells off the window: mean %.4f, ",
           "steps over 1%%: %d of %d | stopped paths: %d\n"),
    nm, seconds[nm], mean(ow, na.rm = TRUE), sum(ow > 0.01, na.rm = TRUE),
    sum(!is.na(ow)), length(attr(paths[[nm]], "stopped"))
  ))
}

out_dir <- if (track$set == "test") "sims/xgb_paths" else "sims/xgb_paths_train"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(out_dir, sprintf("sims_xgb_%s.rds", KEY))
saveRDS(list(key = KEY, season = season, picks = picks, reach = reach,
             times = times, paths = paths, seconds = seconds,
             score_gap = gaps),
        out_path)
cat(sprintf("-> %s   total %.1f min\n", out_path,
            as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
