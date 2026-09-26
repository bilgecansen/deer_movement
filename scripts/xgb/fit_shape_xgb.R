#' @description
#' Fit one shape of the pooled model on one season, score it on held-out
#' folds, and record what its trees reached for.
#'
#' Four shapes differ only in whether the start-of-step block is present
#' and whether the ranked trees may combine variables:
#'
#'                modifier block   hab / famd trees
#'   full            yes           may combine
#'   rsf             no            may combine
#'   start           yes           one variable per tree
#'   main            no            one variable per tree
#'
#' A 2x2, so each effect can be read twice — letting the ranked trees
#' combine is (full - start) and (rsf - main), and adding the start block
#' is (full - rsf) and (start - main). On the pooled seasons the two
#' readings agree to within 3% and the effects add, which is what says they
#' are separable.
#'
#' 500 rounds per BOOSTER in every shape, so `full` and `start` carry 2,500
#' trees and the others 2,000. The block is the unit: at equal rounds per
#' block, habitat is fitted the same way in all four, so a difference
#' between them is down to what actually differs. Equalising total trees
#' would instead hand the four-booster shapes 625 rounds each, giving their
#' habitat blocks 25% more trees than the five-booster ones, and push them
#' past the round count that won on held-out score.
#'
#' Three noise columns calibrate the three ranked blocks — a variable earns
#' its place by beating its block's noise column, not by clearing zero.
#' shadow_start is drawn per step, the others per row (see
#' prep_start_xgb.R).
#'
#' Importance is the drop in held-out conditional log-likelihood when a
#' column is shuffled. End-point columns are shuffled within a stratum;
#' start columns between strata, a whole step at a time, because shuffling
#' a stratum-constant column within its own stratum changes nothing and
#' scores exactly zero. That zero is kept in the output as a check.
#'
#' Inputs: data/xgb/pooled_start_<season>.rds (prep_start_xgb.R)
#' Output: results/xgb/compare_<season>_<shape>.rds
#'
#' Configuration: edit the block below before running, or let
#' run_shapes_xgb.R set SEASON and CONFIG.

# Configuration ---------------------------------------------------------------
SEASON <- "fa"
CONFIG <- "rsf"
# Trees per booster. 500 beats 1000 on held-out score in every scheme
# tried; more rounds buy in-sample fit and lose held-out fit.
N_ROUNDS <- 500L
LEARNING_RATE <- 0.05
MAX_DEPTH <- 2L
# Share of steps each round learns from. A regulariser only; the honest
# score comes from the folds.
BAG_FRAC <- 0.632
N_FOLDS <- 5L
N_PERM <- 3L
SEED <- 1L
N_THREAD <- max(1L, parallel::detectCores() - 1L)
overwrite <- TRUE

# Load packages ---------------------------------------------------------------
library(xgboost)
library(tidyverse)

# helper functions
source("scripts/helper_functions.R")

out_path <- sprintf("results/xgb/compare_%s_%s.rds", SEASON, CONFIG)
dir.create("results/xgb", showWarnings = FALSE, recursive = TRUE)
if (!overwrite && file.exists(out_path)) {
  stop(sprintf("%s exists and overwrite is FALSE", out_path))
}

d <- readRDS(sprintf("data/xgb/pooled_start_%s.rds", SEASON))
# Drawn before anything shape-specific, so every shape of the model on a
# season sees the identical noise columns and the identical folds.
set.seed(SEED)
d$shadow_gauss <- stats::rnorm(nrow(d))
set.seed(SEED + 1L)
d$shadow_famd <- stats::rnorm(nrow(d))

# Which columns go where ------------------------------------------------------
# nb carries no ndvi_end or ndvi_start; both are derived from the columns
# present, so nothing here needs to know that.
MOVE_VARS <- c("sl_", "tod_day", "cos_ta")
HR_VAR <- "HR_center_end"
FAMD_VARS <- c(grep("^famd[0-9]+_end$", names(d), value = TRUE),
               "shadow_famd")
START_VARS <- c(intersect(c("wiscland_start", "ndvi_start",
                            "forest_edge_start"), names(d)),
                "shadow_start")
HAB_VARS <- setdiff(
  names(d),
  c(MOVE_VARS, HR_VAR, FAMD_VARS, START_VARS, "key", "deer", "animal",
    "year", "step_id_", "case_", "stratum")
)
CATEGORICAL <- c("landcover", "wiscland_start")

specs <- make_xgb_specs(
  move_vars = MOVE_VARS,
  hr_var = HR_VAR,
  hab_vars = HAB_VARS,
  famd_vars = FAMD_VARS,
  config = CONFIG,
  start_vars = START_VARS,
  categorical = CATEGORICAL,
  learning_rate = LEARNING_RATE,
  max_depth = MAX_DEPTH,
  nthread = N_THREAD
)
uses_start <- CONFIG %in% c("full", "start")
n_steps <- sum(d$case_ == 1)
n_deer_years <- dplyr::n_distinct(d$key)
cat(sprintf("%s / %s: %s | %d deer-years, %s steps, %d trees\n",
            SEASON, CONFIG,
            paste(vapply(specs, `[[`, "", "name"), collapse = " -> "),
            n_deer_years, formatC(n_steps, format = "d", big.mark = ","),
            length(specs) * N_ROUNDS))
cat(sprintf("habitat: %s\nFAMD: %s\nstart: %s\n",
            paste(HAB_VARS, collapse = ", "),
            paste(FAMD_VARS, collapse = ", "),
            if (uses_start) paste(START_VARS, collapse = ", ") else
              "(not used by this shape)"))

# Fit -------------------------------------------------------------------------
# The whole-data fit is what the interaction heat maps read; the folds
# below give the score and the importance.
t0 <- Sys.time()
set.seed(SEED)
fit <- fit_xgb_boosters(d, specs, n_rounds = N_ROUNDS,
                        bag_frac = BAG_FRAC, verbose_every = 0)
cat(sprintf("full fit %.1f min, in-sample logLik %.1f\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            fit$loglik))

ranked <- which(vapply(specs, function(s) {
  s$name %in% c("hab", "famd", "modifier")
}, logical(1)))
struct <- stats::setNames(
  lapply(ranked, function(b) {
    xgb_tree_structure(xgboost::xgb.load.raw(fit$boosters[[b]]))
  }),
  vapply(specs[ranked], `[[`, character(1), "name")
)

# Importance ------------------------------------------------------------------
imp <- xgb_cv_importance(
  d, specs, n_rounds = N_ROUNDS, bag_frac = BAG_FRAC, n_folds = N_FOLDS,
  n_perm = N_PERM, seed = SEED,
  across_vars = if (uses_start) START_VARS else character(0),
  also_within = uses_start, verbose = TRUE
)
imp$season <- SEASON
imp$config <- CONFIG
imp$n_deer <- n_deer_years
imp$n_steps <- n_steps
imp$per_deer <- xgb_scale_value(imp$cv, n_deer_years, n_steps, "deer")

null_ll <- n_steps * log(1 / (nrow(d) / n_steps))
cat(sprintf("\nlogLik  in-sample %.1f  held-out %.1f  null %.1f\n",
            fit$loglik, attr(imp, "ll_cv"), null_ll))
cat(sprintf("=== %s / %s importance ===\n", SEASON, CONFIG))
print(as.data.frame(imp), row.names = FALSE, digits = 4)

saveRDS(
  list(season = SEASON, config = CONFIG, specs = specs, importance = imp,
       structure = struct, loglik_in = fit$loglik,
       ll_cv = attr(imp, "ll_cv"), null_ll = null_ll, n_steps = n_steps,
       n_deer_years = n_deer_years, n_animals = dplyr::n_distinct(d$animal),
       n_trees = length(specs) * N_ROUNDS),
  out_path
)
cat(sprintf("\n-> %s   total %.1f min\n", out_path,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
