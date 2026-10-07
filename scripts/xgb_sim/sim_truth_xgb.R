#' @description
#' Build the TRUE step-selection models the simulation draws from, and
#' calibrate their strengths to the real data.
#'
#' Every endpoint gets a true score:
#'
#'   score = movement and home-range part  +  habitat part
#'
#' The movement and home-range part is the null model type fitted to the
#' real season, so the simulated deer move like the real ones and the
#' probabilities within a step are as concentrated as they really are. The
#' habitat part is written down here, from functions of the real columns,
#' so what the framework should find is known exactly.
#'
#' Four scenarios:
#'   A  no habitat effect at all
#'   B  main effects only: strong (elevation, forest edge), moderate
#'      (NDVI, landcover), just above the threshold (northness), just below
#'      it (FAMD3), and none (eastness)
#'   C  B plus two interactions: elevation x forest edge, and a stronger
#'      elevation x NDVI
#'   D  B with elevation's smooth hump replaced by a sharp band, which a
#'      tree can only follow by cutting elevation again and again
#'
#' Each variable's TRUE importance is computed the way the framework
#' computes importance, but with the true model: shuffle the variable
#' within each step and measure the expected drop in log-likelihood, per
#' deer-year. The strengths are tuned until those hit TARGET.
#'
#' Input:  data/xgb/pooled_start_<season>.rds
#' Output: sims/xgb/<season>/truth.rds
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
SEASON <- "pf"
N_ROUNDS <- 500L
SEED <- 1L
# True importance to aim for, per deer-year, in scenario B
TARGET <- c(elevation = 15, forest_edge = 15, ndvi = 6, landcover = 6,
            northness = 4, famd3 = 2)
# Each interaction's true worth in scenario C, per deer-year
TARGET_INTERACTION <- c(elev_x_edge = 6, elev_x_ndvi = 15)
N_SHUFFLE <- 3L
N_CALIBRATE <- 8L

# Load packages ---------------------------------------------------------------
library(xgboost)
library(tidyverse)

# helper functions
source("scripts/helper_functions.R")

out_dir <- file.path("sims/xgb", SEASON)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

d <- readRDS(sprintf("data/xgb/pooled_start_%s.rds", SEASON))
sizes <- xgb_strata_sizes(d)
sor <- rep(seq_along(sizes), sizes)
n_deer <- dplyr::n_distinct(d$key)

# The movement and home-range part: the null type fitted to the real data
specs_null <- make_xgb_specs(
  move_vars = c("sl_", "tod_day", "cos_ta"), hr_var = "HR_center_end",
  hab_vars = character(0), config = "null", nthread = 6L
)
set.seed(SEED)
fit_null <- fit_xgb_boosters(d, specs_null, n_rounds = N_ROUNDS)
eta_null <- Reduce(`+`, xgb_contributions(
  lapply(fit_null$boosters, xgboost::xgb.load.raw), d, specs_null
))

# Standardised columns, and the habitat functions written in them, shared
# with the summary through helpers_xgb_sim.R so the two cannot drift apart
Z <- xgb_sim_columns(d)
famd_on <- Z$famd_on
terms_of <- xgb_sim_terms
habitat_score <- function(Z, scenario, coef) {
  tt <- terms_of(Z, scenario)
  Reduce(`+`, lapply(names(coef), function(k) coef[[k]] * tt[[k]]))
}

# Expected log-likelihood, per deer-year, of a model with score `eta` when
# the steps are drawn from the truth p_true
expected_ll <- function(p_true, eta) {
  sum(p_true * log(xgb_softmax(eta, sizes))) / n_deer
}

# A within-step shuffle of one column, as the framework does it. The FAMD
# axis is shuffled only among the endpoints that have it.
shuffle_within <- function(x, rows = seq_along(x)) {
  sh <- rows[order(sor[rows], stats::runif(length(rows)))]
  x[rows] <- x[sh]
  x
}
true_importance <- function(var, scenario, coef, p_true, ll_true) {
  mean(replicate(N_SHUFFLE, {
    Zs <- Z
    Zs[[var]] <- if (var == "famd3") {
      shuffle_within(Z[[var]], which(famd_on))
    } else {
      shuffle_within(Z[[var]])
    }
    ll_true - expected_ll(p_true, eta_null + habitat_score(Zs, scenario,
                                                           coef))
  }))
}
all_importance <- function(scenario, coef) {
  eta <- eta_null + habitat_score(Z, scenario, coef)
  p_true <- xgb_softmax(eta, sizes)
  ll_true <- expected_ll(p_true, eta)
  set.seed(SEED)
  vapply(c(names(TARGET), "eastness"), true_importance, numeric(1),
         scenario = scenario, coef = coef, p_true = p_true,
         ll_true = ll_true)
}

# Calibrate B: each strength is rescaled by the square root of how far its
# importance is from the target, which is roughly how importance scales
calibrate <- function(scenario, coef, targets) {
  for (it in seq_len(N_CALIBRATE)) {
    imp <- all_importance(scenario, coef)[names(targets)]
    ratio <- pmin(pmax(sqrt(targets / pmax(imp, 1e-6)), 0.5), 2)
    coef[names(targets)] <- coef[names(targets)] * ratio
    cat(sprintf("  %s iteration %d: %s\n", scenario, it,
                paste(sprintf("%s %.2f", names(imp), imp),
                      collapse = ", ")))
  }
  coef
}

cat("\ncalibrating B\n")
coef_B <- calibrate("B", setNames(rep(1, length(TARGET)), names(TARGET)),
                    TARGET)
cat("calibrating D (the band replaces the hump)\n")
coef_D <- calibrate("D", coef_B, TARGET["elevation"])

# C: B's strengths plus the interactions, each tuned so that removing it
# alone, with the other kept, costs its TARGET_INTERACTION per deer-year
interaction_worth <- function(term, coef) {
  eta <- eta_null + habitat_score(Z, "C", coef)
  p_true <- xgb_softmax(eta, sizes)
  expected_ll(p_true, eta) -
    expected_ll(p_true, eta - coef[[term]] * terms_of(Z, "C")[[term]])
}
all_worth <- function(coef) {
  vapply(names(TARGET_INTERACTION), interaction_worth, numeric(1),
         coef = coef)
}
coef_C <- c(coef_B, stats::setNames(rep(0.1, length(TARGET_INTERACTION)),
                                    names(TARGET_INTERACTION)))
cat("calibrating C (the interactions)\n")
for (it in seq_len(N_CALIBRATE)) {
  worth <- all_worth(coef_C)
  ratio <- pmin(pmax(sqrt(TARGET_INTERACTION / worth), 0.5), 2)
  coef_C[names(worth)] <- coef_C[names(worth)] * ratio
  cat(sprintf("  C iteration %d: %s\n", it,
              paste(sprintf("%s %.2f", names(worth), worth),
                    collapse = ", ")))
}

scenarios <- list(
  A = list(coef = c(), score = rep(0, nrow(d))),
  B = list(coef = coef_B, score = habitat_score(Z, "B", coef_B)),
  C = list(coef = coef_C, score = habitat_score(Z, "C", coef_C)),
  D = list(coef = coef_D, score = habitat_score(Z, "D", coef_D))
)

# Report: strengths, true importance, and how much habitat adds over the
# null in total, per step, next to the real fits
p_null <- xgb_softmax(eta_null, sizes)
report <- purrr::imap_dfr(scenarios, function(s, nm) {
  eta <- eta_null + s$score
  p_true <- xgb_softmax(eta, sizes)
  imp <- if (nm == "A") {
    stats::setNames(rep(0, length(TARGET) + 1), c(names(TARGET), "eastness"))
  } else {
    all_importance(nm, s$coef)
  }
  maxp <- tapply(p_true, sor, max)
  worth <- if (nm == "C") all_worth(s$coef) else TARGET_INTERACTION * 0
  tibble::tibble(
    scenario = nm,
    habitat_over_null_per_step = sum(p_true * (log(p_true) - log(p_null))) /
      length(sizes),
    median_top_p = stats::median(maxp),
    !!!as.list(worth),
    !!!as.list(imp)
  )
})
cat("\n=== strengths ===\n")
no_interactions <- TARGET_INTERACTION * 0
print(rbind(B = c(coef_B, no_interactions), C = coef_C,
            D = c(coef_D, no_interactions)), digits = 3)
cat("\n=== true values, per deer-year (importance) and per step ===\n")
print(as.data.frame(report), row.names = FALSE, digits = 3)

saveRDS(
  list(season = SEASON, scenarios = scenarios, eta_null = eta_null,
       report = report, target = TARGET,
       target_interaction = TARGET_INTERACTION, n_deer = n_deer),
  file.path(out_dir, "truth.rds")
)
cat(sprintf("\n-> %s\n", file.path(out_dir, "truth.rds")))
