#' @description
#' Fit the pooled boosted-tree model for one season, all years stacked, and
#' score how much each candidate variable matters relative to the others.
#'
#' This is exploratory: it ranks candidate variables so you can decide what
#' is worth carrying into the GAM and amt model sets. It is not a filter,
#' and nothing here is a test — the gates in filter_models_gam.R /
#' filter_models_amt.R do the testing, on simulations and the log score.
#'
#' The model is four boosters grown in alternation on one conditional-logit
#' likelihood (see helpers_xgb.R):
#'   * movement — sl_ with time of day, and cos(ta_);
#'   * HR centre — distance to the home-range centre, on its own;
#'   * habitat — the environmental candidates, plus a column of noise;
#'   * FAMD — the five axes, switched off where they are NA (outside
#'     forest) and anchored at a score of 0, so the forest / non-forest
#'     level stays with landcover instead of being absorbed here.
#'
#' Movement and HR centre are the nuisance block and are never ranked. They
#' mirror the GAM null, so an importance here is what a variable adds on top
#' of that null. Trees in the ranked boosters may combine variables;
#' max_depth 2 keeps any root-to-leaf path pairwise.
#'
#' Importance is the drop in held-out conditional log-likelihood over
#' N_FOLDS random folds of whole steps, when a variable is shuffled within
#' strata (FAMD among forest points only). SHADOW_COL is a column of pure
#' noise carried through the whole fit; it lands at zero and marks where
#' "carries nothing" sits, which is the line to read against rather than
#' zero itself.
#'
#' Years are stacked rather than fitted separately because the GAM and amt
#' model sets are fitted per season, not per season-year. Year is constant
#' within a stratum, so it cancels in the softmax — there is no year term to
#' fit. Pooling also cuts the over-fitting sharply: on fa the in-sample and
#' held-out scores differ by about 1% of the distance above the null,
#' against 7% for a single year.
#'
#' Inputs: data/xgb/pooled_<season>_<year>.rds (prep_pool_xgb.R)
#' Output: results/xgb/season_xgb_<season>.rds
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
SEASON <- "fa"
# Trees per booster. 500 beats 1000 on the held-out score in every scheme
# tried; more rounds buy in-sample fit and lose held-out fit.
N_ROUNDS <- 500L
LEARNING_RATE <- 0.05
MAX_DEPTH <- 2L
# Share of steps each round learns from. A regulariser only; the honest
# score comes from the folds.
BAG_FRAC <- 0.632
N_FOLDS <- 5L
# Permutations averaged per variable per fold
N_PERM <- 3L
# Name of the noise column added to the habitat block
SHADOW_COL <- "shadow_gauss"
SEED <- 1L
N_THREAD <- max(1L, parallel::detectCores() - 1L)
overwrite <- TRUE

# Load packages ---------------------------------------------------------------
library(xgboost)
library(tidyverse)

# helper functions
source("scripts/helper_functions.R")

out_path <- sprintf("results/xgb/season_xgb_%s.rds", SEASON)
dir.create("results/xgb", showWarnings = FALSE, recursive = TRUE)
if (!overwrite && file.exists(out_path)) {
  stop(sprintf("%s exists and overwrite is FALSE", out_path))
}

pooled <- xgb_season_pool(SEASON)
set.seed(SEED)
pooled[[SHADOW_COL]] <- stats::rnorm(nrow(pooled))

n_steps <- sum(pooled$case_ == 1)
# A deer tracked in two years counts once per year, so these are deer-years
# rather than animals; both are reported because the per-deer scale is per
# deer-year.
per_year <- pooled |>
  filter(case_ == 1) |>
  group_by(year) |>
  summarise(deer = n_distinct(animal), steps = n(), .groups = "drop")
n_deer_years <- sum(per_year$deer)
n_animals <- n_distinct(pooled$animal)
cat(sprintf("season %s: %d deer-years, %d animals, %s observed steps\n",
            SEASON, n_deer_years, n_animals,
            formatC(n_steps, format = "d", big.mark = ",")))
print(as.data.frame(per_year), row.names = FALSE)

# Which columns go where ------------------------------------------------------
# nb carries no ndvi_end column: it is dropped when the per-year files are
# built, because the winter season is named for the year it starts in and
# runs into the next, so that year's NDVI stack does not cover it. HAB_VARS
# is derived from the columns present, so nothing here needs to know that.
MOVE_VARS <- c("sl_", "tod_day", "cos_ta")
HR_VAR <- "HR_center_end"
FAMD_VARS <- grep("^famd[0-9]+_end$", names(pooled), value = TRUE)
HAB_VARS <- setdiff(
  names(pooled),
  c(MOVE_VARS, HR_VAR, FAMD_VARS, "key", "deer", "animal", "year",
    "step_id_", "case_", "stratum")
)
cat(sprintf("\nhabitat (%d): %s\nFAMD (%d): %s\n",
            length(HAB_VARS), paste(HAB_VARS, collapse = ", "),
            length(FAMD_VARS), paste(FAMD_VARS, collapse = ", ")))

specs <- make_xgb_specs(
  move_vars = MOVE_VARS,
  hr_var = HR_VAR,
  hab_vars = HAB_VARS,
  famd_vars = FAMD_VARS,
  categorical = "landcover",
  learning_rate = LEARNING_RATE,
  max_depth = MAX_DEPTH,
  nthread = N_THREAD
)

# Fit -------------------------------------------------------------------------
t0 <- Sys.time()
set.seed(SEED)
fit <- fit_xgb_boosters(pooled, specs, n_rounds = N_ROUNDS,
                        bag_frac = BAG_FRAC, verbose_every = 100)
cat(sprintf("\nfit in %.1f min, in-sample logLik %.1f\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")),
            fit$loglik))

# Importance ------------------------------------------------------------------
imp <- xgb_cv_importance(pooled, specs, n_rounds = N_ROUNDS,
                         bag_frac = BAG_FRAC, n_folds = N_FOLDS,
                         n_perm = N_PERM, seed = SEED, verbose = TRUE)
imp$season <- SEASON
imp$n_deer <- n_deer_years
imp$n_steps <- n_steps
imp$per_deer <- xgb_scale_value(imp$cv, n_deer_years, n_steps, "deer")
imp$per_100 <- xgb_scale_value(imp$cv, n_deer_years, n_steps, "steps100")

null_ll <- n_steps * log(1 / (nrow(pooled) / n_steps))
cat(sprintf("\nlogLik  in-sample %.1f  held-out %.1f  null %.1f\n",
            fit$loglik, attr(imp, "ll_cv"), null_ll))
cat(sprintf("above null: in-sample %.0f, held-out %.0f (gap %.1f%%)\n",
            fit$loglik - null_ll, attr(imp, "ll_cv") - null_ll,
            100 * (fit$loglik - attr(imp, "ll_cv")) /
              (fit$loglik - null_ll)))
cat(sprintf("\n=== importance, pooled %s ===\n", SEASON))
print(as.data.frame(imp), row.names = FALSE, digits = 4)

# What the trees reached for, free once the fit is here.
structure_of <- function(b) {
  xgb_tree_structure(xgboost::xgb.load.raw(fit$boosters[[b]]))
}
ranked <- which(vapply(specs, function(s) s$name %in% c("hab", "famd"),
                       logical(1)))
struct <- stats::setNames(lapply(ranked, structure_of),
                          vapply(specs[ranked], `[[`, character(1), "name"))

saveRDS(
  list(fit = fit, specs = specs, importance = imp, structure = struct,
       season = SEASON, n_deer_years = n_deer_years, n_animals = n_animals,
       n_steps = n_steps, per_year = per_year, ll_cv = attr(imp, "ll_cv"),
       null_ll = null_ll),
  out_path
)
cat(sprintf("\n-> %s   total %.1f min\n", out_path,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
