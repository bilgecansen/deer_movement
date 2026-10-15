#' @description
#' Build the pooled step tables the xgboost importance run reads: one per
#' season and training year, every deer-year of the canonical list's training
#' set (canonical_deer(), set == "train") in it. The test year stays out.
#'
#' The random points are the wrangle's nonp design, read from data/tracks/
#' (stp.var.nonp): 100 per step, uniform in distance out to the deer's
#' longest step with a turning angle, every one counting the same. The xgb
#' path draws no points of its own.
#'
#' A habitat column present on fewer than MIN_PRESENT of a season's rows is
#' left out for that season, for every year alike, so the years still stack.
#' ndvi_end in nb is the case in point: the winter season runs into a year
#' whose NDVI stack it does not load.
#'
#' Columns written: the movement block (sl_, tod_hour, cos_ta), the habitat
#' variables in HAB_VARS, landcover as integer codes for xgboost's
#' categorical splits, and the FAMD axes. FAMD stays NA off forest; nothing
#' is filled in, and no model type uses it.
#'
#' Output: data/xgb/pooled_<season>_<year>.rds, a list of
#'   * pooled — the step table, sorted by stratum with the observed step first
#'   * filter — steps in and steps kept for every deer-year
#' A season's files for other years, left from an earlier build, are moved
#' to data/xgb/archive/ so xgb_season_pool() cannot stack them in.
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
SEASONS <- c("br", "fa", "nb", "pf")
# Habitat variables; landcover enters as one categorical column
HAB_VARS <- c(
  "HR_center_end",
  "ndvi_end",
  "landcover",
  "forest_edge_end",
  "elevation_end",
  "northness_end",
  "eastness_end"
)
# FAMD axes; NA off forest, so they are fit through the switch
FAMD_VARS <- sprintf("famd%d_end", 1:5)
# A habitat column is kept for a season only if present on at least this
# share of its rows
MIN_PRESENT <- 0.95

# Load packages ---------------------------------------------------------------
library(tidyverse)

# helper functions
source("scripts/helper_functions.R")

dir.create("data/xgb", showWarnings = FALSE, recursive = TRUE)
dir.create("data/xgb/archive", showWarnings = FALSE, recursive = TRUE)

train <- canonical_deer(readRDS("library/SW_filtered_deer.RData"))
train <- train[train$set == "train", ]
train$key <- sprintf("%s_%s_%d", train$id, train$season,
                     as.integer(train$year))
missing <- train$key[!file.exists(sprintf("data/tracks/data_%s.rds",
                                          train$key))]
if (length(missing)) {
  stop(sprintf("No wrangled track for %d training deer-years, e.g. %s",
               length(missing), paste(utils::head(missing, 3),
                                      collapse = ", ")))
}

for (season in SEASONS) {
  cat(sprintf("\n########## %s ##########\n", season))
  keys <- train$key[train$season == season]
  years <- sort(unique(as.integer(train$year[train$season == season])))

  # Which habitat columns the season carries: rows with a value, over all
  # of the season's rows
  counts <- sapply(keys, function(k) {
    v <- readRDS(sprintf("data/tracks/data_%s.rds", k))$stp.var.nonp[[1]]
    lc <- as.integer(factor(v$wiscland_end, levels = LANDCOVER_LEVELS))
    c(vapply(HAB_VARS, function(h) {
      x <- if (h == "landcover") lc else v[[h]]
      sum(!is.na(x))
    }, numeric(1)), rows = nrow(v))
  })
  present <- rowSums(counts[HAB_VARS, , drop = FALSE]) / sum(counts["rows", ])
  hab_vars <- HAB_VARS[present >= MIN_PRESENT]
  cat("habitat columns, share present:\n")
  print(round(present, 3))
  cat(sprintf("keeping: %s\n", paste(hab_vars, collapse = ", ")))

  for (yr in years) {
    yr_keys <- keys[as.integer(train$year[match(keys, train$key)]) == yr]
    pooled <- build_xgb_pool(yr_keys, hab_vars, FAMD_VARS)

    on <- !is.na(pooled[[FAMD_VARS[1]]])
    stopifnot(all(vapply(FAMD_VARS,
                         function(v) identical(!is.na(pooled[[v]]), on),
                         logical(1))))

    filter_tab <- tibble::tibble(
      key = yr_keys,
      steps_in = vapply(yr_keys, function(k) {
        v <- readRDS(sprintf("data/tracks/data_%s.rds", k))$stp.var.nonp[[1]]
        sum(v$case_)
      }, numeric(1)),
      steps_kept = vapply(yr_keys, function(k) {
        sum(pooled$case_[pooled$key == k])
      }, numeric(1))
    )
    out_path <- sprintf("data/xgb/pooled_%s_%d.rds", season, yr)
    saveRDS(list(pooled = pooled, filter = filter_tab), out_path)
    cat(sprintf(
      "%s %d: %d deer-years | %d of %d steps kept | %d rows -> %s\n",
      season, yr, length(yr_keys), sum(filter_tab$steps_kept),
      sum(filter_tab$steps_in), nrow(pooled), out_path
    ))
  }

  # Files for this season's other years would be stacked in by
  # xgb_season_pool(); they go to the archive
  old <- list.files("data/xgb", sprintf("^pooled_%s_[0-9]{4}[.]rds$", season),
                    full.names = TRUE)
  stale <- old[!as.integer(sub(".*_([0-9]{4})[.]rds$", "\\1", old)) %in% years]
  if (length(stale)) {
    file.rename(stale, file.path("data/xgb/archive", basename(stale)))
    cat(sprintf("archived: %s\n", paste(basename(stale), collapse = ", ")))
  }
}
