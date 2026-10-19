#' @description
#' Fit one model type on a season's whole training table, without the noise
#' columns, and save its trees for simulating paths.
#'
#' The two noise columns calibrate importance in model selection
#' (fit_model_xgb.R); a model used to simulate paths has no use for them,
#' so they are left out here. Everything else is as in fit_model_xgb.R: the
#' same columns, blocks, rounds, depths, bagging and seed. The in-sample
#' score is printed next to that of the selection run's whole-data fit
#' (results/xgb/compare_<season>_<type>.rds) for comparison.
#'
#' run_final_xgb.R picks which model types to fit and runs this script for
#' each of them.
#'
#' Inputs: data/xgb/pooled_start_<season>.rds (prep_start_xgb.R),
#'         results/xgb/compare_<season>_<type>.rds (fit_model_xgb.R)
#' Output: results/xgb/final_<season>_<type>.rds — the trees, the block
#'         specs, the deer-id levels and the in-sample score
#'
#' Configuration: edit the block below before running, or let
#' run_final_xgb.R set SEASON, CONFIG and N_THREAD.

# Configuration ---------------------------------------------------------------
SEASON <- "fa"
CONFIG <- "main"
# The same settings as fit_model_xgb.R
N_ROUNDS <- 500L
LEARNING_RATE <- 0.05
MAX_DEPTH <- 2L
MAX_DEPTH_DEER <- 3L
BAG_FRAC <- 0.632
SEED <- 1L
N_THREAD <- max(1L, parallel::detectCores() - 1L)

# Load packages ---------------------------------------------------------------
library(xgboost)
library(tidyverse)

# helper functions
source("scripts/helper_functions.R")

out_path <- sprintf("results/xgb/final_%s_%s.rds", SEASON, CONFIG)
d <- readRDS(sprintf("data/xgb/pooled_start_%s.rds", SEASON))

# The deer types' animal id, coded 0..k-1 for xgboost's categorical splits
with_deer <- grepl("_deer$", CONFIG)
deer_levels <- NULL
if (with_deer) {
  deer_levels <- sort(unique(d$animal))
  d$deer_id <- match(d$animal, deer_levels) - 1L
}

# Which columns go where ------------------------------------------------------
# As in fit_model_xgb.R, less the two noise columns
MOVE_VARS <- c("sl_", "tod_hour", "cos_ta")
HR_VAR <- "HR_center_end"
NOISE_VARS <- c("shadow_gauss", "shadow_start")
FAMD_COLS <- grep("^famd[0-9]+_end$", names(d), value = TRUE)
START_VARS <- intersect(c("wiscland_start", "ndvi_start",
                          "forest_edge_log_start"),
                        names(d))
HAB_VARS <- setdiff(
  names(d),
  c(MOVE_VARS, HR_VAR, FAMD_COLS, START_VARS, NOISE_VARS, "key", "deer",
    "animal", "year", "step_id_", "case_", "stratum", "day_of_season",
    "deer_id")
)
CATEGORICAL <- c("landcover", "wiscland_start", "deer_id")

specs <- make_xgb_specs(
  move_vars = MOVE_VARS,
  hr_var = HR_VAR,
  hab_vars = HAB_VARS,
  config = CONFIG,
  start_vars = START_VARS,
  categorical = CATEGORICAL,
  learning_rate = LEARNING_RATE,
  max_depth = MAX_DEPTH,
  nthread = N_THREAD,
  deer_var = "deer_id",
  deer_depth = MAX_DEPTH_DEER
)
stopifnot(!any(NOISE_VARS %in% unlist(lapply(specs, `[[`, "feats"))))
cat(sprintf("%s / %s (final, no noise columns): %d trees\n", SEASON, CONFIG,
            length(specs) * N_ROUNDS))
cat(sprintf("  %-8s %s\n", vapply(specs, `[[`, "", "name"),
            vapply(specs, function(s) paste(s$feats, collapse = ", "), "")),
    sep = "")

# Fit -------------------------------------------------------------------------
t0 <- Sys.time()
set.seed(SEED)
fit <- fit_xgb_boosters(d, specs, n_rounds = N_ROUNDS, bag_frac = BAG_FRAC,
                        verbose_every = 0)
selection_fit <- readRDS(sprintf("results/xgb/compare_%s_%s.rds", SEASON,
                                 CONFIG))
cat(sprintf("in-sample logLik %.1f (selection run: %.1f)\n",
            fit$loglik, selection_fit$loglik_in))

saveRDS(
  list(season = SEASON, config = CONFIG, specs = specs,
       boosters = fit$boosters, n_trees = fit$n_trees,
       loglik_in = fit$loglik, loglik_in_selection = selection_fit$loglik_in,
       deer_levels = deer_levels, n_steps = sum(d$case_ == 1),
       n_deer_years = dplyr::n_distinct(d$key)),
  out_path
)
cat(sprintf("-> %s   %.1f min\n", out_path,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
