#' @description
#' Repeat the whole importance run several times to see how much of each
#' variable's score is just the procedure's own randomness.
#'
#' The deer and their observed steps are fixed — these 41 deer are the
#' population of interest, so there is no sampling to bootstrap. What changes
#' between replicates is everything the procedure draws:
#'   * a fresh set of N_PTS random points per observed step;
#'   * a fresh bagging and column-sampling seed;
#'   * fresh permutations inside the importance step.
#'
#' Variables whose interval sits clearly above zero carry signal. Variables
#' straddling zero are ones this procedure cannot tell from noise: shuffling
#' them changes an out-of-bag score about as often up as down.
#'
#' Replicates run in parallel, one per worker. Each worker loads its own
#' rasters, since they cannot be shipped between processes.
#'
#' Runs every season-year in RUNS. Variables are set per season, matching
#' run_all_xgb.R: winter (`nb`) leaves NDVI out, because the season runs into
#' the next calendar year while load_ndvi() stacks one year.
#'
#' Inputs: data/xgb/pooled_<season>_<year>.rds, for the list of deer that
#'   passed the steps-per-edf filter (its points are not reused)
#' Outputs:
#'   results/xgb/replicates_xgb_<season>_<year>.rds  every replicate's scores
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
# Season-years to run
RUNS <- tibble::tribble(
  ~season, ~year,
  "fa", 2017L,
  "fa", 2018L,
  "fa", 2019L,
  "fa", 2020L,
  "fa", 2021L,
  "nb", 2017L,
  "nb", 2018L,
  "nb", 2020L,
  "nb", 2021L,
  "pf", 2017L,
  "pf", 2018L
)
# Habitat variables per season; landcover enters as one categorical column
HAB_VARS_BY_SEASON <- list(
  fa = c("HR_center_end", "ndvi_end", "landcover", "forest_edge_end",
         "elevation_end", "northness_end", "eastness_end"),
  pf = c("HR_center_end", "ndvi_end", "landcover", "forest_edge_end",
         "elevation_end", "northness_end", "eastness_end"),
  nb = c("HR_center_end", "landcover", "forest_edge_end",
         "elevation_end", "northness_end", "eastness_end")
)
REPLICATES <- 10L
# Replicate r uses seed SEED_BASE + r, so a run can be reproduced
SEED_BASE <- 1000L
# These must match fit_xgb.R, or the replicates describe a different model
N_PTS <- 100L
N_ROUNDS <- 500L
LEARNING_RATE <- 0.05
MAX_DEPTH <- 2L
ONE_COL_PER_TREE <- TRUE
BAG_FRAC <- 0.632
N_PERM <- 3L
# Workers, and threads per worker. A replicate costs roughly 1.4 minutes per
# 1,000 observed steps, so the big winters are the slow ones.
N_WORKERS <- 5L
N_THREAD <- 2L
MIN_STEPS_PER_EDF <- 10
# Where finished fits leave their marker files
PROGRESS_DIR <- "results/xgb/replicate_progress"
overwrite <- TRUE

# Load packages ---------------------------------------------------------------
library(amt)
library(terra)
library(tidyverse)
library(sf)
library(xgboost)
library(furrr)

# helper functions
source("scripts/helper_functions.R")

dir.create("results/xgb", showWarnings = FALSE, recursive = TRUE)

MOVE_VARS <- c("sl_", "tod_day", "cos_ta")
FAMD_VARS <- sprintf("famd%d_end", 1:5)

tracks <- readRDS("library/SW_filtered_deer.RData") |>
  dplyr::filter(keep == TRUE, excursion == FALSE,
                unstable_hr_center == FALSE) |>
  dplyr::filter(year %in% 2017:2021)
tracks$key <- sprintf("%s_%s_%d", tracks$id, tracks$season, tracks$year)

# The deer that passed the filter, from each season-year's stored pool.
deer_keys <- purrr::map2(RUNS$season, RUNS$year, function(season, year) {
  f <- readRDS(sprintf("data/xgb/pooled_%s_%d.rds", season, year))$filter
  f$key[f$keep]
})
names(deer_keys) <- sprintf("%s_%d", RUNS$season, RUNS$year)

# One replicate of one season-year ---------------------------------------------
one_replicate <- function(season, year, rep_id) {
  source("scripts/helper_functions.R")
  key <- sprintf("%s_%d", season, year)
  hab_vars <- HAB_VARS_BY_SEASON[[season]]
  set.seed(SEED_BASE + rep_id)
  t0 <- Sys.time()

  landcover <- load_landcover(year, season)
  rasters <- list(
    landcover = landcover,
    water = make_water_mask(landcover),
    ndvi = load_ndvi(year),
    landfire = load_landfire(year, season),
    topo = load_topo()
  )
  pooled <- build_xgb_pool(
    keys = deer_keys[[key]], tracks = tracks, rasters = rasters,
    n_pts = N_PTS, hab_vars = hab_vars, famd_vars = FAMD_VARS,
    move_vars = MOVE_VARS
  )
  specs <- make_xgb_specs(
    move_vars = MOVE_VARS, hab_vars = hab_vars, famd_vars = FAMD_VARS,
    categorical = "landcover", learning_rate = LEARNING_RATE,
    max_depth = MAX_DEPTH, one_col_per_tree = ONE_COL_PER_TREE,
    nthread = N_THREAD
  )
  fit <- fit_xgb_boosters(pooled, specs, N_ROUNDS, bag_frac = BAG_FRAC)
  imp <- xgb_perm_importance(fit$boosters, pooled, specs, bag = fit$bag,
                             n_perm = N_PERM)
  imp$season <- season
  imp$year <- year
  imp$rep <- rep_id
  imp$n_deer <- dplyr::n_distinct(pooled$deer)
  imp$n_steps <- max(pooled$stratum)
  imp$loglik <- fit$loglik
  imp$secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Worker output is not shown, so each finished fit leaves a marker file:
  # the run is long, and this is how you see it progressing.
  dir.create(PROGRESS_DIR, showWarnings = FALSE, recursive = TRUE)
  cat(sprintf("%s rep %d: %.0f s, %d steps\n", key, rep_id, imp$secs[1],
              imp$n_steps[1]),
      file = file.path(PROGRESS_DIR, sprintf("%s_rep%02d.txt", key, rep_id)))
  imp
}

# Run -------------------------------------------------------------------------
# Every (season-year, replicate) pair is one task, so the workers stay busy
# even though the season-years differ a lot in size.
tasks <- tidyr::expand_grid(RUNS, rep_id = seq_len(REPLICATES))
start_time <- Sys.time()
cat(sprintf("%d season-years x %d replicates = %d fits on %d workers\n",
            nrow(RUNS), REPLICATES, nrow(tasks), N_WORKERS))
future::plan(future::multisession, workers = N_WORKERS)

reps <- furrr::future_pmap_dfr(
  list(tasks$season, tasks$year, tasks$rep_id),
  one_replicate,
  .options = furrr::furrr_options(
    packages = c("amt", "terra", "tidyverse", "sf", "xgboost"),
    seed = TRUE
  ),
  # No progress bar: these runs are long enough to be logged to a file, and
  # the bar redraws thousands of lines there.
  .progress = FALSE
)

future::plan(future::sequential)
for (k in unique(sprintf("%s_%d", reps$season, reps$year))) {
  saveRDS(reps[sprintf("%s_%d", reps$season, reps$year) == k, ],
          sprintf("results/xgb/replicates_xgb_%s.rds", k))
}

# Summary ---------------------------------------------------------------------
summ <- reps |>
  dplyr::group_by(season, year, variable) |>
  dplyr::summarise(
    median = median(out_of_bag),
    p10 = quantile(out_of_bag, 0.1),
    p90 = quantile(out_of_bag, 0.9),
    share_positive = mean(out_of_bag > 0),
    .groups = "drop"
  ) |>
  dplyr::arrange(season, year, dplyr::desc(median))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("\n%d fits in %.1f min (%.0f s each)\n", nrow(tasks), elapsed,
            median(reps$secs)))
cat("\nout-of-bag importance across replicates:\n")
print(as.data.frame(
  summ |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(., 1)),
                  share_positive = round(share_positive, 2))
), row.names = FALSE)
cat(sprintf("\n-> results/xgb/replicates_xgb_<season>_<year>.rds (%d files)\n",
            dplyr::n_distinct(sprintf("%s_%d", reps$season, reps$year))))
