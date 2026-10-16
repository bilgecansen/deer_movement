#' @description
#' The log-score filter: the one-step log score of every deer-year's real
#' path under each of its season's selected deer models, against the deer
#' null. Like the HR overlap and svf filters it takes a path and a model, so
#' the training years and the test year are scored the same way.
#'
#' The models are the final ones (fit_final_xgb.R), fitted on all the
#' training years: the season's two deer models (final_selection.rds, the
#' simplest and the best) and the deer null (null_deer: movement and the
#' home-range centre, with the animal id). At every real step,
#' xgb_onestep_logscore() builds the simulator's kernel at the real start,
#' with the reach the paths were simulated with (REACH, see
#' sim_paths_xgb.R), and takes the log probability of the cell where the
#' deer ended. Before scoring, xgb_path_check() confirms the kernel reads
#' the deer's own steps as the models were fitted on them.
#'
#' A step counts only if every model could score it: one that ends beyond
#' the reach, in its own start cell, or in a cell with no weight under any
#' model is left out of every model's total. A model passes for a deer-year
#' when its total beats the null's by at least DELTA_MIN.
#'
#' Deer-years are scored in parallel, one at a time per worker, the longest
#' first, and each is saved as it finishes, so an interrupted run resumes
#' where it stopped.
#'
#' Inputs:  data/tracks/data_<key>.rds, results/xgb/final_selection.rds,
#'          results/xgb/final_<season>_<type>.rds (the two picks and
#'          null_deer), the rasters the wrangle read
#' Output:  filters/xgb/logscore_<tag>/logscore_xgb_<key>.rds, one per
#'          deer-year (every step's score under every model, and the
#'          totals), and filters/xgb/logscore_xgb_<tag>.rds, the totals of
#'          all of them; <tag> is xgb_scores_tag(SET, REACH)
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
# The deer-years to score: "test" or "train"
SET <- "train"
# The reach rule, as the paths were simulated: "longest" or "gamma99"
REACH <- "gamma99"
NULL_TYPE <- "null_deer"
# A model must beat the null by this many log units over a deer-year
DELTA_MIN <- 3
# Deer-years to score: NULL for every one of the set not yet scored
KEYS <- NULL
N_PAR <- parallel::detectCores() - 1L

# Load packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(furrr)
})

# helper functions
source("scripts/helper_functions.R")

tag <- xgb_scores_tag(SET, REACH)
out_dir <- sprintf("filters/xgb/logscore_%s", tag)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (is.null(KEYS)) {
  deer <- canonical_deer(readRDS("library/SW_filtered_deer.RData"))
  deer <- deer[deer$set == SET, ]
  KEYS <- sprintf("%s_%s_%d", deer$id, deer$season, as.integer(deer$year))
  # The longest deer-years first
  n_steps <- vapply(deer$stp, nrow, integer(1))
  KEYS <- KEYS[order(-n_steps)]
  KEYS <- KEYS[!file.exists(file.path(out_dir,
                                      sprintf("logscore_xgb_%s.rds", KEYS)))]
}
cat(sprintf("scoring %d %s deer-years on %d workers (%s)\n", length(KEYS),
            tag, N_PAR, format(Sys.time(), "%H:%M")))

score_one <- function(key) {
  source("scripts/helper_functions.R")
  track <- readRDS(sprintf("data/tracks/data_%s.rds", key))
  stopifnot(track$set == SET)
  season <- track$season
  stp <- track$stp[[1]]
  picks <- readRDS("results/xgb/final_selection.rds")
  picks <- picks[picks$season == season & picks$deer, ]
  configs <- c(picks$config, NULL_TYPE)
  models <- lapply(sprintf("results/xgb/final_%s_%s.rds", season, configs),
                   xgb_path_model, nthread = 1L)
  names(models) <- configs
  deer_ids <- lapply(models, function(m) {
    i <- match(track$id, m$deer_levels)
    stopifnot(!is.na(i))
    i - 1L
  })

  reach_longest <- xgb_path_reach(track$id, season)
  reach <- if (REACH == "longest") {
    reach_longest
  } else {
    min(xgb_path_reach_gamma(track), reach_longest)
  }
  land <- xgb_path_landscape(track, load_landcover(track$year, season),
                             load_ndvi(track$year), load_topo())
  v <- track$stp.var.nonp[[1]]
  layers <- stats::na.omit(unique(xgb_path_ndvi_layer(
    c(stp$t1_, stp$t2_, v$t2_), land$ndvi_times
  )))
  statics <- Map(function(m, id) xgb_path_static(m, land, layers, id),
                 models, deer_ids)
  gaps <- xgb_path_check(track, land, models, statics, deer_ids)

  steps <- purrr::imap_dfr(models, function(m, cf) {
    xgb_onestep_logscore(m, land, statics[[cf]], stp, reach,
                         deer_ids[[cf]]) |>
      dplyr::mutate(config = cf)
  })
  # The steps every model could score
  common <- steps |>
    dplyr::group_by(step) |>
    dplyr::summarise(all_ok = all(status == "ok"), .groups = "drop") |>
    dplyr::filter(all_ok) |>
    dplyr::pull(step)
  # Reach and the start cell are the same for every model, so the null's
  # count of those steps is every model's; a cell with no weight can differ
  # by model
  null_status <- table(factor(steps$status[steps$config == NULL_TYPE],
                              c("ok", "beyond_reach", "own_cell",
                                "no_weight")))
  n_no_weight <- dplyr::n_distinct(steps$step[steps$status == "no_weight"])
  totals <- steps |>
    dplyr::filter(step %in% common) |>
    dplyr::group_by(config) |>
    dplyr::summarise(logp = sum(logp), .groups = "drop")
  logp_null <- totals$logp[totals$config == NULL_TYPE]
  summary <- totals |>
    dplyr::filter(config != NULL_TYPE) |>
    dplyr::left_join(picks[, c("config", "best", "simplest")],
                     by = "config") |>
    dplyr::mutate(
      set = SET, reach_rule = REACH, reach = reach, key = key,
      season = season, n_steps = sum(!is.na(stp$ta_)),
      n_scored = length(common),
      n_beyond_reach = null_status[["beyond_reach"]],
      n_own_cell = null_status[["own_cell"]],
      n_no_weight = n_no_weight,
      logp_null = logp_null, delta_logp = logp - logp_null,
      pass = delta_logp >= DELTA_MIN, check_gap = max(gaps)
    ) |>
    dplyr::relocate(set, reach_rule, reach, key, season, config, best,
                    simplest)
  saveRDS(list(key = key, season = season, reach = reach, steps = steps,
               summary = summary),
          file.path(out_dir, sprintf("logscore_xgb_%s.rds", key)))
  "ok"
}

if (length(KEYS)) {
  future::plan(future::multisession, workers = min(N_PAR, length(KEYS)))
  status <- furrr::future_map_chr(
    KEYS,
    function(key) {
      tryCatch(score_one(key),
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

files <- list.files(out_dir, pattern = "^logscore_xgb_.*[.]rds$",
                    full.names = TRUE)
totals <- purrr::map_dfr(files, function(f) readRDS(f)$summary)
out_path <- sprintf("filters/xgb/logscore_xgb_%s.rds", tag)
saveRDS(totals, out_path)
cat(sprintf("-> %s: %d deer-years, %d rows; pass: delta_logp >= %g (%s)\n",
            out_path, dplyr::n_distinct(totals$key), nrow(totals), DELTA_MIN,
            format(Sys.time(), "%H:%M")))
