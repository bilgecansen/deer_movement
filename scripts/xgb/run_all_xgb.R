#' @description
#' Run the pooled importance fit for several season-year combinations in one
#' go: build each pool, fit, score the variables, and write the same files
#' prep_pool_xgb.R and fit_xgb.R write for a single one.
#'
#' Variables are set per season. Winter (`nb`) leaves NDVI out: the season is
#' named for the year it starts in and runs into the next, while load_ndvi()
#' stacks one year, so ndvi_end is NA for 71-98% of winter rows and keeping it
#' would drop those rows from the fit.
#'
#' Season-years run in parallel, one per worker; each worker loads its own
#' rasters, since they cannot be shipped between processes.
#'
#' Outputs, per season-year:
#'   data/xgb/pooled_<season>_<year>.rds
#'   results/xgb/fit_xgb_<season>_<year>.rds
#'   results/xgb/importance_xgb_<season>_<year>.rds
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
FAMD_VARS <- sprintf("famd%d_end", 1:5)
MOVE_VARS <- c("sl_", "tod_day", "cos_ta")
# Keep a deer when observed steps / edf of its GAM null is at least this
MIN_STEPS_PER_EDF <- 10
# These match fit_xgb.R
N_PTS <- 100L
N_ROUNDS <- 500L
LEARNING_RATE <- 0.05
MAX_DEPTH <- 2L
ONE_COL_PER_TREE <- TRUE
BAG_FRAC <- 0.632
N_PERM <- 3L
N_WORKERS <- 4L
N_THREAD <- 3L

# Load packages ---------------------------------------------------------------
library(amt)
library(terra)
library(tidyverse)
library(sf)
library(mgcv)
library(xgboost)
library(furrr)

# helper functions
source("scripts/helper_functions.R")

dir.create("data/xgb", showWarnings = FALSE, recursive = TRUE)
dir.create("results/xgb", showWarnings = FALSE, recursive = TRUE)

tracks_all <- readRDS("library/SW_filtered_deer.RData") |>
  dplyr::filter(keep == TRUE, excursion == FALSE,
                unstable_hr_center == FALSE) |>
  dplyr::filter(year %in% 2017:2021)
tracks_all$key <- sprintf("%s_%s_%d", tracks_all$id, tracks_all$season,
                          tracks_all$year)

# One season-year -------------------------------------------------------------
run_one <- function(season, year) {
  source("scripts/helper_functions.R")
  key <- sprintf("%s_%d", season, year)
  t0 <- Sys.time()
  set.seed(year)

  hab_vars <- HAB_VARS_BY_SEASON[[season]]

  # Deer with enough steps per effective degree of freedom of their GAM null
  null_files <- list.files(
    "results/gam",
    pattern = sprintf("^results_gam_null_.*_%s[.]rds$", key),
    full.names = TRUE
  )
  filter_tab <- purrr::map_dfr(null_files, function(f) {
    gfit <- readRDS(f)$gam
    deer_key <- sub("^results_gam_null_(.*)[.]rds$", "\\1", basename(f))
    if (!inherits(gfit, "gam")) {
      return(tibble::tibble(key = deer_key, n_steps = NA_integer_,
                            edf = NA_real_))
    }
    tibble::tibble(key = deer_key,
                   n_steps = length(unique(gfit$model[[1]][, 2])),
                   edf = sum(gfit$edf))
  }) |>
    dplyr::mutate(
      steps_per_edf = n_steps / edf,
      keep = !is.na(steps_per_edf) & steps_per_edf >= MIN_STEPS_PER_EDF
    )
  keys <- filter_tab$key[filter_tab$keep]

  landcover <- load_landcover(year, season)
  rasters <- list(
    landcover = landcover,
    water = make_water_mask(landcover),
    ndvi = load_ndvi(year),
    landfire = load_landfire(year, season),
    topo = load_topo()
  )
  pooled <- build_xgb_pool(
    keys = keys, tracks = tracks_all, rasters = rasters, n_pts = N_PTS,
    hab_vars = hab_vars, famd_vars = FAMD_VARS, move_vars = MOVE_VARS
  )
  saveRDS(list(pooled = pooled, filter = filter_tab),
          sprintf("data/xgb/pooled_%s.rds", key))

  specs <- make_xgb_specs(
    move_vars = MOVE_VARS, hab_vars = hab_vars, famd_vars = FAMD_VARS,
    categorical = "landcover", learning_rate = LEARNING_RATE,
    max_depth = MAX_DEPTH, one_col_per_tree = ONE_COL_PER_TREE,
    nthread = N_THREAD
  )
  fit <- fit_xgb_boosters(pooled, specs, N_ROUNDS, bag_frac = BAG_FRAC)
  saveRDS(
    list(boosters = fit$boosters, specs = specs, n_trees = fit$n_trees,
         loglik = fit$loglik, bag = fit$bag,
         settings = list(n_rounds = N_ROUNDS, learning_rate = LEARNING_RATE,
                         max_depth = MAX_DEPTH, bag_frac = BAG_FRAC,
                         one_col_per_tree = ONE_COL_PER_TREE,
                         n_perm = N_PERM)),
    sprintf("results/xgb/fit_xgb_%s.rds", key)
  )

  imp <- xgb_perm_importance(fit$boosters, pooled, specs, bag = fit$bag,
                             n_perm = N_PERM)
  imp$season <- season
  imp$year <- year
  imp$n_deer <- dplyr::n_distinct(pooled$deer)
  imp$n_steps <- max(pooled$stratum)
  saveRDS(imp, sprintf("results/xgb/importance_xgb_%s.rds", key))

  tibble::tibble(
    key = key, deer = dplyr::n_distinct(pooled$deer),
    steps = max(pooled$stratum), rows = nrow(pooled),
    loglik = fit$loglik,
    mins = as.numeric(difftime(Sys.time(), t0, units = "mins"))
  )
}

# Run -------------------------------------------------------------------------
start_time <- Sys.time()
cat(sprintf("%d season-years on %d workers\n", nrow(RUNS), N_WORKERS))
future::plan(future::multisession, workers = N_WORKERS)

done <- furrr::future_map2_dfr(
  RUNS$season, RUNS$year, run_one,
  .options = furrr::furrr_options(
    packages = c("amt", "terra", "tidyverse", "sf", "mgcv", "xgboost"),
    seed = TRUE
  ),
  .progress = FALSE
)

future::plan(future::sequential)

cat(sprintf("\nfinished in %.1f min\n",
            difftime(Sys.time(), start_time, units = "mins")))
print(as.data.frame(
  done |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(., 1)))
), row.names = FALSE)
