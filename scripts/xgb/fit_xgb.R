#' @description
#' Fit the pooled boosted-tree model for one season and year and score how
#' much each variable matters relative to the others.
#'
#' This is exploratory: it ranks candidate variables so you can decide what
#' is worth carrying into the GAM and amt model sets. It is not a filter, and
#' nothing here is held out — the gates in filter_models_gam.R /
#' filter_models_amt.R do the testing, on simulations and the log score.
#'
#' The model is three boosters grown in alternation on one conditional-logit
#' likelihood (see helpers_xgb.R):
#'   * movement — sl_ with time of day, and cos(ta_); always available;
#'   * habitat — HR_center_end, ndvi_end, landcover; single-variable trees,
#'     one column offered per tree;
#'   * FAMD — the five axes, switched off where they are NA (outside forest)
#'     and anchored at a score of 0, so the forest / non-forest level stays
#'     with landcover instead of being absorbed here.
#'
#' Importance is the drop in conditional log-likelihood when a variable is
#' shuffled within strata (FAMD among forest points only) — the same units as
#' delta_logp. Two versions are reported:
#'   * out_of_bag — each step scored only by the trees from rounds that left
#'     it out (BAG_FRAC sets how many). Rank by this one;
#'   * in_sample — every step scored by the whole model, which also counts
#'     the noise the trees fitted.
#' On fa 2021 the two agree on the top variables, while FAMD and aspect fall
#' to zero or below out-of-bag, and elevation drops sharply.
#'
#' Inputs: data/xgb/pooled_<season>_<year>.rds (prep_pool_xgb.R)
#' Outputs:
#'   results/xgb/fit_xgb_<season>_<year>.rds        boosters + settings
#'   results/xgb/importance_xgb_<season>_<year>.rds importance table
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
SEASON <- "fa"
YEAR <- 2021L
# Trees per booster. Moderate on purpose: in-sample importance inflates with
# rounds, and more for continuous variables than categorical ones.
N_ROUNDS <- 500L
LEARNING_RATE <- 0.05
MAX_DEPTH <- 2L
# Offer each habitat / FAMD tree one randomly chosen column, so a dominant
# variable cannot take every tree.
ONE_COL_PER_TREE <- TRUE
# Share of steps each round learns from; the rest are that round's
# out-of-bag steps. 1 disables bagging and leaves only in-sample importance.
BAG_FRAC <- 0.632
# Permutations averaged per variable
N_PERM <- 3L
N_THREAD <- max(1L, parallel::detectCores() - 1L)
overwrite <- TRUE

# Load packages ---------------------------------------------------------------
library(xgboost)
library(tidyverse)

# helper functions
source("scripts/helper_functions.R")

key <- sprintf("%s_%d", SEASON, YEAR)
fit_path <- sprintf("results/xgb/fit_xgb_%s.rds", key)
imp_path <- sprintf("results/xgb/importance_xgb_%s.rds", key)
dir.create("results/xgb", showWarnings = FALSE, recursive = TRUE)

if (!overwrite && file.exists(fit_path)) {
  stop(sprintf("%s exists and overwrite is FALSE", fit_path))
}

pooled <- readRDS(sprintf("data/xgb/pooled_%s.rds", key))$pooled

# Which columns go where ------------------------------------------------------
MOVE_VARS <- c("sl_", "tod_day", "cos_ta")
FAMD_VARS <- grep("^famd[0-9]+_end$", names(pooled), value = TRUE)
HAB_VARS <- setdiff(
  names(pooled),
  c(MOVE_VARS, FAMD_VARS, "key", "deer", "step_id_", "case_", "stratum")
)
cat(sprintf("habitat: %s\n", paste(HAB_VARS, collapse = ", ")))
cat(sprintf("FAMD: %s\n", paste(FAMD_VARS, collapse = ", ")))

specs <- make_xgb_specs(
  move_vars = MOVE_VARS,
  hab_vars = HAB_VARS,
  famd_vars = FAMD_VARS,
  categorical = "landcover",
  learning_rate = LEARNING_RATE,
  max_depth = MAX_DEPTH,
  one_col_per_tree = ONE_COL_PER_TREE,
  nthread = N_THREAD
)

# Fit -------------------------------------------------------------------------
start_time <- Sys.time()
set.seed(YEAR)
cat(sprintf(
  "fitting %d rounds on %d deer, %d steps, %d rows\n", N_ROUNDS,
  dplyr::n_distinct(pooled$deer), max(pooled$stratum), nrow(pooled)
))

fit <- fit_xgb_boosters(pooled, specs, N_ROUNDS, bag_frac = BAG_FRAC,
                        verbose_every = 100)

saveRDS(
  list(
    boosters = fit$boosters,
    specs = specs,
    n_trees = fit$n_trees,
    loglik = fit$loglik,
    bag = fit$bag,
    settings = list(n_rounds = N_ROUNDS, learning_rate = LEARNING_RATE,
                    max_depth = MAX_DEPTH, bag_frac = BAG_FRAC,
                    one_col_per_tree = ONE_COL_PER_TREE, n_perm = N_PERM)
  ),
  fit_path
)

# Importance ------------------------------------------------------------------
imp <- xgb_perm_importance(fit$boosters, pooled, specs, bag = fit$bag,
                           n_perm = N_PERM)
saveRDS(imp, imp_path)

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf(
  "\nin-sample logLik %.1f   elapsed %.1f min\n", fit$loglik, elapsed
))
cat("\nimportance (drop in conditional log-likelihood when shuffled),",
    "ranked by out-of-bag:\n")
print(
  imp |> dplyr::mutate(dplyr::across(where(is.numeric), ~ round(., 1))),
  row.names = FALSE
)
cat(sprintf("\n-> %s\n-> %s\n", fit_path, imp_path))
