#' @description
#' Build the pooled step tables the hierarchical GAM reads: one per season,
#' every deer-year of the canonical list's training set (canonical_deer(),
#' set == "train") stacked. The test year stays out.
#'
#' The random points are the wrangle's amt design, read from data/tracks/
#' (stp.var): drawn from the deer-year's own fitted gamma (step length) and
#' von Mises (turning angle), 25 per step. N_PTS of them are kept per step,
#' at random. The 25 are independent draws from the same distributions, so
#' keeping N_PTS of them is the same as having drawn N_PTS. They are picked
#' before any step is dropped, so which ones are kept does not depend on the
#' covariates.
#'
#' A step with any listed variable missing, on the observed end or on a kept
#' point, is dropped whole. A habitat column present on fewer than
#' MIN_PRESENT of a season's rows is left out for that season: ndvi_end in
#' nb, as in prep_pool_xgb.R.
#'
#' Columns written: key, animal, year, step_id_, case_, stratum, the
#' movement columns (sl_, log_sl, ta_, cos_ta, tod_), the habitat columns
#' kept, wiscland_end, and one 0/1 indicator per landcover class other than
#' forest, the reference (<class>_end, in LANDCOVER_LEVELS order).
#'
#' Output: data/gam/pooled_<season>.rds, a list of
#'   * pooled   — the step table, sorted by stratum with the observed step
#'                first
#'   * filter   — steps in and steps kept for every deer-year
#'   * hab_vars — the habitat columns kept for the season
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
SEASONS <- "fa"
# Random points kept per step, of the 25 the wrangle drew
N_PTS <- 10L
HAB_VARS <- c(
  "HR_center_end",
  "ndvi_end",
  "forest_edge_end",
  "elevation_end",
  "northness_end",
  "eastness_end"
)
# A habitat column is kept for a season only if present on at least this
# share of its rows
MIN_PRESENT <- 0.95
SEED <- 20261016L

# Load packages ---------------------------------------------------------------
library(tidyverse)

# helper functions
source("scripts/helper_functions.R")

dir.create("data/gam", showWarnings = FALSE, recursive = TRUE)

train <- canonical_deer(readRDS("library/SW_filtered_deer.RData"))
train <- train[train$set == "train", ]
train$key <- sprintf("%s_%s_%d", train$id, train$season,
                     as.integer(train$year))
lc_cols <- paste0(LANDCOVER_LEVELS[-1], "_end")
move_cols <- c("sl_", "log_sl", "ta_", "cos_ta", "tod_")

# One deer-year's steps, each with N_PTS of its random points
read_steps <- function(key, id) {
  v <- as.data.frame(
    readRDS(sprintf("data/tracks/data_%s.rds", key))$stp.var[[1]]
  )
  kept <- v |>
    dplyr::filter(!case_) |>
    dplyr::slice_sample(n = N_PTS, by = step_id_)
  v <- dplyr::bind_rows(dplyr::filter(v, case_), kept)
  lc <- factor(v$wiscland_end, levels = LANDCOVER_LEVELS)
  out <- tibble::tibble(
    key = key,
    animal = id,
    year = as.integer(sub(".*_", "", key)),
    step_id_ = v$step_id_,
    case_ = as.integer(v$case_),
    sl_ = v$sl_,
    log_sl = log(v$sl_),
    ta_ = v$ta_,
    cos_ta = cos(v$ta_),
    tod_ = v$tod_
  )
  for (h in HAB_VARS) {
    out[[h]] <- v[[h]]
  }
  out$wiscland_end <- lc
  for (i in seq_along(lc_cols)) {
    ind <- as.integer(lc == LANDCOVER_LEVELS[i + 1])
    # The wrangle's own indicator columns agree wherever the class is known
    known <- !is.na(lc)
    stopifnot(all(v[[lc_cols[i]]][known] == ind[known]))
    out[[lc_cols[i]]] <- ind
  }
  out
}

for (season in SEASONS) {
  cat(sprintf("\n########## %s ##########\n", season))
  set.seed(SEED)
  sub <- train[train$season == season, ]
  sub <- sub[order(sub$key), ]
  raw <- purrr::map2_dfr(sub$key, sub$id, read_steps)

  present <- vapply(HAB_VARS, function(h) mean(!is.na(raw[[h]])), numeric(1))
  hab_vars <- HAB_VARS[present >= MIN_PRESENT]
  cat("habitat columns, share present:\n")
  print(round(present, 3))
  cat(sprintf("keeping: %s\n", paste(hab_vars, collapse = ", ")))
  raw <- raw[, setdiff(names(raw), setdiff(HAB_VARS, hab_vars))]

  pooled <- raw |>
    dplyr::group_by(key, step_id_) |>
    dplyr::filter(
      !any(dplyr::if_any(dplyr::all_of(c(move_cols, hab_vars,
                                         "wiscland_end")), is.na)),
      sum(case_) == 1,
      dplyr::n() >= 2
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(key, step_id_, dplyr::desc(case_)) |>
    dplyr::mutate(
      stratum = as.integer(factor(paste(key, step_id_),
                                  unique(paste(key, step_id_))))
    ) |>
    dplyr::relocate(stratum, .after = case_)

  filter_tab <- raw |>
    dplyr::group_by(key) |>
    dplyr::summarise(steps_in = sum(case_), .groups = "drop") |>
    dplyr::left_join(
      pooled |>
        dplyr::group_by(key) |>
        dplyr::summarise(steps_kept = sum(case_), .groups = "drop"),
      by = "key"
    ) |>
    dplyr::mutate(steps_kept = dplyr::coalesce(steps_kept, 0L))

  rows_per_step <- table(table(pooled$stratum))
  out_path <- sprintf("data/gam/pooled_%s.rds", season)
  saveRDS(list(pooled = pooled, filter = filter_tab, hab_vars = hab_vars),
          out_path)
  cat(sprintf(
    paste0("%s: %d deer, %d deer-years | %d of %d steps kept | %d rows ",
           "| rows per step: %s -> %s\n"),
    season, dplyr::n_distinct(pooled$animal), nrow(filter_tab),
    sum(filter_tab$steps_kept), sum(filter_tab$steps_in), nrow(pooled),
    paste(sprintf("%s (%d)", names(rows_per_step), rows_per_step),
          collapse = ", "),
    out_path
  ))
}
