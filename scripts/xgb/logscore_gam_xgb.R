#' @description
#' The log-score filters against per-deer-year GAMs: the one-step log score
#' of every deer-year's real path under its own GAM, set against the scores
#' of its season's two selected deer models. GAM picks which:
#'   "null_gam"        the GAM null (fit_null_gam_xgb.R)
#'   "gilbertson_gam"  the Gilbertson et al. (2025) habitat model
#'                     (fit_gilbertson_gam_xgb.R)
#'
#' The GAM scoring a deer-year is the one fitted to the training deer-year
#' closest in time (xgb_closest_train_key()): the deer-year itself in
#' training, and for a test deer-year a training year of the same animal and
#' season, with that year's gamma and von Mises. The kernel reads the scored
#' deer-year's own landscape, home-range centre included.
#'
#' The xgb models' scores are the ones logscore_xgb.R saved for each real
#' step. The GAM is scored here on the same cells, with the same reach and
#' the same statuses (xgb_gam_onestep_logscore()); only its cell weight
#' differs, because it was fitted on points drawn from a gamma and von Mises
#' rather than evenly in distance. Before scoring, xgb_gam_check() confirms
#' the kernel reads every endpoint column the GAM reads as the wrangle
#' stored it for the scored deer-year's own random points.
#'
#' A step counts only if the two models and the GAM could all score it. A
#' model passes for a deer-year when its total beats the GAM's by at least
#' DELTA_MIN.
#'
#' Deer-years are scored in parallel, one at a time per worker, and each is
#' saved as it finishes, so an interrupted run resumes where it stopped.
#'
#' Inputs:  data/tracks/data_<key>.rds, results/xgb/<GAM>/,
#'          results/xgb/final_selection.rds,
#'          filters/xgb/logscore_<tag>/logscore_xgb_<key>.rds
#'          (logscore_xgb.R), the rasters the wrangle read
#' Output:  filters/xgb/logscore_<GAM>_<tag>/logscore_<GAM>_<key>.rds, one
#'          per deer-year (the GAM's score at every step, and the totals),
#'          and filters/xgb/logscore_<GAM>_<tag>.rds, the totals of all of
#'          them; <tag> is xgb_scores_tag(SET, REACH)
#'
#' Configuration: edit the block below before running. The GAM and the set
#' can also be given on the command line:
#'   Rscript scripts/xgb/logscore_gam_xgb.R null_gam test

# Configuration ---------------------------------------------------------------
# The GAM to score against: "null_gam" or "gilbertson_gam"
GAM <- "gilbertson_gam"
# The deer-years to score: "test" or "train"
SET <- "train"
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1) {
  GAM <- .args[1]
}
if (length(.args) >= 2) {
  SET <- .args[2]
}
stopifnot(SET %in% c("test", "train"))
# The reach rule, as the paths were simulated: "longest" or "gamma99"
REACH <- "gamma99"
# A model must beat the GAM by this many log units over a deer-year
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

stopifnot(GAM %in% c("null_gam", "gilbertson_gam"))
tag <- xgb_scores_tag(SET, REACH)
xgb_dir <- sprintf("filters/xgb/logscore_%s", tag)
out_dir <- sprintf("filters/xgb/logscore_%s_%s", GAM, tag)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (is.null(KEYS)) {
  deer <- canonical_deer(readRDS("library/SW_filtered_deer.RData"))
  deer <- deer[deer$set == SET, ]
  KEYS <- sprintf("%s_%s_%d", deer$id, deer$season, as.integer(deer$year))
  # The longest deer-years first
  n_steps <- vapply(deer$stp, nrow, integer(1))
  gam_keys <- mapply(xgb_closest_train_key, deer$id, deer$season,
                     as.integer(deer$year))
  KEYS <- KEYS[order(-n_steps)]
  gam_keys <- gam_keys[order(-n_steps)]
  # A deer-year whose GAM did not fit has nothing to score
  fitted <- file.exists(sprintf("results/xgb/%s/%s_%s.rds", GAM, GAM,
                                gam_keys))
  if (any(!fitted)) {
    cat(sprintf("%d deer-years have no fitted %s and are not scored\n",
                sum(!fitted), GAM))
  }
  KEYS <- KEYS[fitted]
  KEYS <- KEYS[!file.exists(file.path(
    out_dir, sprintf("logscore_%s_%s.rds", GAM, KEYS)
  ))]
}
cat(sprintf("scoring %d %s deer-years against %s on %d workers (%s)\n",
            length(KEYS), tag, GAM, N_PAR, format(Sys.time(), "%H:%M")))

score_one <- function(key) {
  source("scripts/helper_functions.R")
  track <- readRDS(sprintf("data/tracks/data_%s.rds", key))
  stopifnot(track$set == SET)
  season <- track$season
  picks <- readRDS("results/xgb/final_selection.rds")
  picks <- picks[picks$season == season, ]
  xgb <- readRDS(file.path(xgb_dir, sprintf("logscore_xgb_%s.rds", key)))

  reach <- xgb_path_reach_rule(track, REACH)
  stopifnot(reach == xgb$reach)
  land <- xgb_path_landscape(track, load_landcover(track$year, season),
                             load_ndvi(track$year), load_topo())
  gam_key <- xgb_closest_train_key(track$id, season, as.integer(track$year))
  m <- readRDS(sprintf("results/xgb/%s/%s_%s.rds", GAM, GAM, gam_key))
  end_vars <- xgb_gam_end_vars(m$fit)
  # Every one of the scored deer-year's own random points
  gap <- xgb_gam_check(xgb_gam_data(track, 25L, 1L, end_vars), land,
                       end_vars)

  gam_steps <- xgb_gam_onestep_logscore(m$fit, land, track$stp[[1]], reach,
                                        m$sl_distr, m$ta_distr) |>
    dplyr::mutate(config = GAM)
  steps <- dplyr::bind_rows(
    dplyr::filter(xgb$steps, config %in% picks$config),
    gam_steps
  )
  stopifnot(setequal(steps$step[steps$config == GAM],
                     steps$step[steps$config == picks$config[1]]))
  # The steps every model and the GAM could score
  common <- steps |>
    dplyr::group_by(step) |>
    dplyr::summarise(all_ok = all(status == "ok"), .groups = "drop") |>
    dplyr::filter(all_ok) |>
    dplyr::pull(step)
  gam_status <- table(factor(gam_steps$status,
                             c("ok", "beyond_reach", "own_cell",
                               "no_weight")))
  totals <- steps |>
    dplyr::filter(step %in% common) |>
    dplyr::group_by(config) |>
    dplyr::summarise(logp = sum(logp), .groups = "drop")
  logp_null <- totals$logp[totals$config == GAM]
  summary <- totals |>
    dplyr::filter(config != GAM) |>
    dplyr::left_join(picks[, c("config", "best", "simplest")],
                     by = "config") |>
    dplyr::mutate(
      set = SET, reach_rule = REACH, reach = reach, key = key,
      season = season, null = GAM, gam_key = gam_key,
      n_steps = sum(!is.na(track$stp[[1]]$ta_)),
      n_scored = length(common),
      n_beyond_reach = gam_status[["beyond_reach"]],
      n_own_cell = gam_status[["own_cell"]],
      n_no_weight = dplyr::n_distinct(steps$step[steps$status ==
                                                   "no_weight"]),
      logp_null = logp_null, delta_logp = logp - logp_null,
      pass = delta_logp >= DELTA_MIN, check_gap = gap
    ) |>
    dplyr::relocate(set, reach_rule, reach, key, season, config, best,
                    simplest, null)
  saveRDS(list(key = key, season = season, reach = reach,
               steps = gam_steps, summary = summary),
          file.path(out_dir, sprintf("logscore_%s_%s.rds", GAM, key)))
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

files <- list.files(out_dir, pattern = sprintf("^logscore_%s_.*[.]rds$", GAM),
                    full.names = TRUE)
totals <- purrr::map_dfr(files, function(f) readRDS(f)$summary)
out_path <- sprintf("filters/xgb/logscore_%s_%s.rds", GAM, tag)
saveRDS(totals, out_path)
cat(sprintf("-> %s: %d deer-years, %d rows; pass: delta_logp >= %g (%s)\n",
            out_path, dplyr::n_distinct(totals$key), nrow(totals), DELTA_MIN,
            format(Sys.time(), "%H:%M")))
