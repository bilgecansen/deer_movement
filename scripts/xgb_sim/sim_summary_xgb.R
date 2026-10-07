#' @description
#' Compare what the framework found on simulated data with the truth it was
#' simulated from. Works on whatever replicates have finished, so it can be
#' run part-way through run_sim_xgb.R.
#'
#' Five comparisons, per scenario:
#'
#'   selection   which model type the rule picks (simplest within THRESHOLD
#'               per deer of the best). The right answer: null in A, main
#'               in B and D, an interaction type in C.
#'   signal      how much of the true model's log-likelihood the fitted
#'               models reach, and what share of the true habitat signal,
#'               over the null, they capture.
#'   variables   each variable's estimated importance against its true
#'               importance, and how often it clears THRESHOLD against
#'               whether it truly does. Read from the selected model, as the
#'               real workflow does, and from main on its own.
#'   pairs       the interaction log score of every pair in the rsf model:
#'               each of the two true pairs in C (its score, how often it
#'               clears THRESHOLD, its rank among all pairs), and pairs that
#'               do not exist anywhere, which is where the approximate H
#'               could invent them.
#'   curves      the fitted effect of elevation and of forest edge in the
#'               main model against the true function, both centred: their
#'               correlation, and the error relative to the true curve's
#'               spread.
#'
#' Inputs:  sims/xgb/<season>/truth.rds,
#'          sims/xgb/<season>/<scenario>/rep_<k>/,
#'          data/xgb/pooled_start_<season>.rds
#' Output:  sims/xgb/<season>/summary.rds, and the tables printed
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
SEASON <- "pf"
THRESHOLD <- 3
COMPLEXITY <- c("null", "main", "rsf", "rsf_hr", "start", "full")
EXPECTED <- c(A = "null", B = "main", C = "rsf", D = "main")
# Rows averaged over, and grid points, for the fitted curves
CURVE_ROWS <- 3000L
CURVE_GRID <- 25L
SEED <- 1L

# Load packages ---------------------------------------------------------------
library(xgboost)
library(tidyverse)

# helper functions
source("scripts/helper_functions.R")

season_dir <- file.path("sims/xgb", SEASON)
truth <- readRDS(file.path(season_dir, "truth.rds"))
d <- readRDS(sprintf("data/xgb/pooled_start_%s.rds", SEASON))
sizes <- xgb_strata_sizes(d)
n_deer <- truth$n_deer
p_null <- xgb_softmax(truth$eta_null, sizes)

# Every finished fit
files <- list.files(season_dir, pattern = "^compare_.*[.]rds$",
                    recursive = TRUE, full.names = TRUE)
if (!length(files)) {
  stop("No simulation results yet; run run_sim_xgb.R")
}
fits <- tibble::tibble(
  file = files,
  scenario = basename(dirname(dirname(files))),
  rep = basename(dirname(files)),
  config = sub(sprintf("^compare_%s_(.*)[.]rds$", SEASON), "\\1",
               basename(files))
)
fits$ll_cv <- vapply(fits$file, function(f) readRDS(f)$ll_cv, numeric(1))
# Only replicates with every type fit take part
n_types <- dplyr::n_distinct(fits$config)
fits <- fits |>
  dplyr::group_by(scenario, rep) |>
  dplyr::filter(dplyr::n() == n_types) |>
  dplyr::ungroup()
cat(sprintf("complete replicates: %s\n", paste(
  sprintf("%s %d", names(table(unique(fits[c("scenario", "rep")])$scenario)),
          table(unique(fits[c("scenario", "rep")])$scenario)),
  collapse = ", "
)))

# 1. Selection -----------------------------------------------------------------
selection <- fits |>
  dplyr::group_by(scenario, rep) |>
  dplyr::mutate(per_deer = (ll_cv - max(ll_cv)) / n_deer) |>
  dplyr::filter(per_deer > -THRESHOLD) |>
  dplyr::slice_min(match(config, COMPLEXITY), n = 1) |>
  dplyr::ungroup() |>
  dplyr::select(scenario, rep, selected = config)
cat("\n=== 1. selection: how often each type is picked ===\n")
print(as.data.frame(
  selection |>
    dplyr::count(scenario, selected) |>
    tidyr::pivot_wider(names_from = selected, values_from = n,
                       values_fill = 0) |>
    dplyr::mutate(expected = EXPECTED[scenario])
), row.names = FALSE)

# 2. Signal captured -----------------------------------------------------------
signal <- fits |>
  dplyr::group_by(scenario, rep) |>
  dplyr::summarise(best = max(ll_cv), null_fit = ll_cv[config == "null"],
                   .groups = "drop") |>
  dplyr::mutate(
    choice = file.path(season_dir, scenario, rep, "choice.rds"),
    ll_true = purrr::map2_dbl(choice, scenario, function(f, s) {
      used <- readRDS(f)
      p <- xgb_softmax(truth$eta_null + truth$scenarios[[s]]$score, sizes)
      sum(log(p[used]))
    }),
    ll_true_null = purrr::map_dbl(choice, function(f) {
      sum(log(p_null[readRDS(f)]))
    }),
    gap_to_truth = (ll_true - best) / n_deer,
    habitat_true = (ll_true - ll_true_null) / n_deer,
    habitat_found = (best - null_fit) / n_deer
  )
cat("\n=== 2. signal, per deer-year (mean over replicates) ===\n")
print(as.data.frame(
  signal |>
    dplyr::group_by(scenario) |>
    dplyr::summarise(dplyr::across(c(gap_to_truth, habitat_true,
                                     habitat_found), mean),
                     .groups = "drop")
), row.names = FALSE, digits = 3)

# 3. Variables -----------------------------------------------------------------
TRUE_NAME <- c(elevation_end = "elevation", forest_edge_end = "forest_edge",
               ndvi_end = "ndvi", landcover = "landcover",
               northness_end = "northness", eastness_end = "eastness",
               famd3_end = "famd3")
true_imp <- truth$report |>
  dplyr::select(scenario, dplyr::all_of(unname(TRUE_NAME))) |>
  tidyr::pivot_longer(-scenario, names_to = "truth_name",
                      values_to = "true")
importance_of <- function(f) {
  r <- readRDS(f)
  r$importance |>
    dplyr::filter(scheme == "within") |>
    dplyr::transmute(variable, est = cv / r$n_deer_years)
}
est <- fits |>
  dplyr::filter(config == "main") |>
  dplyr::mutate(imp = purrr::map(file, importance_of)) |>
  tidyr::unnest(imp) |>
  dplyr::select(scenario, rep, variable, main = est) |>
  dplyr::left_join(
    selection |>
      dplyr::left_join(fits, by = c("scenario", "rep",
                                    "selected" = "config")) |>
      dplyr::mutate(imp = purrr::map(file, importance_of)) |>
      tidyr::unnest(imp) |>
      dplyr::select(scenario, rep, variable, selected = est),
    by = c("scenario", "rep", "variable")
  ) |>
  # a replicate that selects the null keeps nothing
  dplyr::mutate(selected = dplyr::coalesce(selected, 0),
                truth_name = TRUE_NAME[variable]) |>
  dplyr::left_join(true_imp, by = c("scenario", "truth_name")) |>
  # columns with no true effect: the other FAMD axes and the noise
  dplyr::mutate(true = dplyr::coalesce(true, 0))
variables <- est |>
  dplyr::group_by(scenario, variable) |>
  dplyr::summarise(
    true = dplyr::first(true),
    main_mean = mean(main),
    main_sd = stats::sd(main),
    clears_main = mean(main >= THRESHOLD),
    clears_selected = mean(selected >= THRESHOLD),
    .groups = "drop"
  ) |>
  dplyr::mutate(truly_clears = true >= THRESHOLD) |>
  dplyr::arrange(scenario, dplyr::desc(true), dplyr::desc(main_mean))
cat("\n=== 3. variables, per deer-year: share of replicates clearing",
    THRESHOLD, "===\n")
print(as.data.frame(variables), row.names = FALSE, digits = 3)

# 4. Pairs ---------------------------------------------------------------------
pairs_of <- function(f) {
  r <- readRDS(f)
  r$split_gain$pairs |>
    dplyr::filter(parent != child) |>
    dplyr::mutate(a = pmin(parent, child), b = pmax(parent, child)) |>
    dplyr::group_by(booster, a, b) |>
    dplyr::summarise(score = sum(heldout) / r$n_deer_years,
                     .groups = "drop")
}
# The pairs that carry a true interaction in C, named as in the truth, with
# their variables in the alphabetical order pairs_of() gives them
TRUE_PAIRS <- c(elev_x_edge = "elevation_end x forest_edge_end",
                elev_x_ndvi = "elevation_end x ndvi_end")
pair_tab <- fits |>
  dplyr::filter(config == "rsf") |>
  dplyr::mutate(pp = purrr::map(file, pairs_of)) |>
  tidyr::unnest(pp) |>
  dplyr::group_by(scenario, rep) |>
  dplyr::mutate(rank = rank(-score)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    true_name = ifelse(scenario == "C",
                       names(TRUE_PAIRS)[match(paste(a, "x", b),
                                               TRUE_PAIRS)],
                       NA_character_),
    is_true = !is.na(true_name),
    with_elevation = a == "elevation_end" | b == "elevation_end"
  )
true_pairs <- pair_tab |>
  dplyr::filter(is_true) |>
  dplyr::select(scenario, rep, true_name, score, rank)
pairs <- pair_tab |>
  dplyr::group_by(scenario, rep) |>
  dplyr::summarise(
    best_false_pair = max(score[!is_true]),
    false_pairs_over = sum(score[!is_true] >= THRESHOLD),
    best_false_with_elevation = max(score[!is_true & with_elevation]),
    .groups = "drop"
  )
cat("\n=== 4a. true pairs in the rsf model, per deer-year (mean over",
    "replicates; top = share where it is the highest pair) ===\n")
print(as.data.frame(
  true_pairs |>
    dplyr::group_by(scenario, true_name) |>
    dplyr::summarise(
      score_mean = mean(score),
      score_min = min(score),
      score_max = max(score),
      clears = mean(score >= THRESHOLD),
      top = mean(rank == 1),
      .groups = "drop"
    )
), row.names = FALSE, digits = 3)
cat("\n=== 4b. false pairs in the rsf model, per deer-year (mean over",
    "replicates; over = share with a false pair clearing", THRESHOLD,
    ") ===\n")
print(as.data.frame(
  pairs |>
    dplyr::group_by(scenario) |>
    dplyr::summarise(
      best_false_pair = mean(best_false_pair),
      false_over = mean(false_pairs_over > 0),
      best_false_with_elevation = mean(best_false_with_elevation),
      .groups = "drop"
    )
), row.names = FALSE, digits = 3)

# 5. Curves --------------------------------------------------------------------
# The main model's habitat trees use one variable each, so its fitted effect
# of a variable is read off by moving that variable alone over a grid,
# averaging over a fixed sample of rows.
set.seed(SEED)
sample_rows <- sample.int(nrow(d), CURVE_ROWS)
# The habitat trees also read the noise column fit_model_xgb.R adds; any
# draw from the same distribution serves for averaging over rows
d$shadow_gauss <- stats::rnorm(nrow(d))
curve_vars <- c(elevation = "elevation_end",
                forest_edge = "forest_edge_end")
curve_fit <- function(f, scenario) {
  r <- readRDS(f)
  nm <- vapply(r$specs, `[[`, "", "name")
  spec <- r$specs[[which(nm == "hab")]]
  hab <- xgboost::xgb.load.raw(r$boosters[[which(nm == "hab")]])
  X <- as.matrix(d[sample_rows, spec$feats])
  coef <- truth$scenarios[[scenario]]$coef
  purrr::imap_dfr(curve_vars, function(col, tn) {
    x <- d[[col]]
    grid <- stats::quantile(x, seq(0.025, 0.975, length.out = CURVE_GRID),
                            names = FALSE)
    fitted <- vapply(grid, function(g) {
      Xg <- X
      Xg[, col] <- g
      mean(xgb_raw(hab, xgb_matrix(Xg, spec$categorical), 1,
                   xgboost::xgb.get.num.boosted.rounds(hab)))
    }, numeric(1))
    zg <- (grid - mean(x)) / stats::sd(x)
    # xgb_sim_terms() computes every term but only this one is read, so
    # the other columns are placeholders of the grid's length
    zlist <- list(elevation = zg, forest_edge = zg, ndvi = zg,
                  landcover = integer(length(zg)), northness = zg,
                  famd3 = zg)
    true <- coef[[tn]] * xgb_sim_terms(zlist, scenario)[[tn]]
    fitted <- fitted - mean(fitted)
    true <- true - mean(true)
    tibble::tibble(
      term = tn,
      correlation = stats::cor(fitted, true),
      error_vs_spread = sqrt(mean((fitted - true)^2)) / stats::sd(true)
    )
  })
}
# A has no true curve, so until another scenario finishes there is nothing
# to compare
curve_fits <- fits |>
  dplyr::filter(config == "main", scenario != "A")
curves <- tibble::tibble(scenario = character(), term = character(),
                         correlation = numeric(),
                         error_vs_spread = numeric())
if (nrow(curve_fits)) {
  curves <- curve_fits |>
    dplyr::mutate(cc = purrr::map2(file, scenario, curve_fit)) |>
    tidyr::unnest(cc) |>
    dplyr::group_by(scenario, term) |>
    dplyr::summarise(correlation = mean(correlation),
                     error_vs_spread = mean(error_vs_spread),
                     .groups = "drop")
}
cat("\n=== 5. fitted curves in the main model against the truth",
    "(mean over replicates) ===\n")
print(as.data.frame(curves), row.names = FALSE, digits = 3)

saveRDS(list(selection = selection, signal = signal, variables = variables,
             estimates = est, pairs = pairs, true_pairs = true_pairs,
             pair_table = pair_tab, curves = curves),
        file.path(season_dir, "summary.rds"))
cat(sprintf("\n-> %s\n", file.path(season_dir, "summary.rds")))
