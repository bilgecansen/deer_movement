#' @description
#' Stack each season's years into one table and attach the start-of-step
#' covariates the model-shape comparison needs.
#'
#' Start covariates are properties of the step's origin, identical for all
#' 101 candidates and independent of which random points were drawn, so
#' they join onto the pooled tables from the stored extraction
#' (stp.var.nonp) rather than needing anything regenerated. That is minutes
#' rather than the hours a full re-extraction would cost.
#'
#' Two things are added beyond the joined columns:
#'
#'   shadow_start  one N(0,1) draw per STEP, repeated across its stratum.
#'                 It has to be drawn per step, not per row: a step-level
#'                 column is constant within a stratum, so a split on it
#'                 partitions whole strata, whose gradients sum to zero, and
#'                 its gain at the root is exactly zero. It can only be used
#'                 below a split on something that varies within a stratum.
#'                 A per-row shadow has no such restriction and would win
#'                 root splits the real start columns cannot contest.
#'   animal        the stable id, from the key prefix; `deer` is only an
#'                 index within each year's file.
#'
#' A start column that is mostly missing for a season is dropped for that
#' season. ndvi_start is about 80% missing in nb for the same reason
#' ndvi_end is absent there: the winter season is named for the year it
#' starts in and runs into the next, and that year's NDVI stack does not
#' cover it.
#'
#' Inputs:  data/xgb/pooled_<season>_<year>.rds (prep_pool_xgb.R)
#'          data/tracks/data_<key>.rds
#' Output:  data/xgb/pooled_start_<season>.rds
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
SEASONS <- c("fa", "nb", "pf")
# Numeric start columns to carry across; the cover type is always taken
NUMERIC_START <- c("ndvi_start", "forest_edge_start")
# A start column is kept for a season only if present on at least this
# share of that season's rows
MIN_PRESENT <- 0.95
SEED <- 1L

# Load packages ---------------------------------------------------------------
library(tidyverse)

# helper functions
source("scripts/helper_functions.R")

start_of_key <- function(k) {
  f <- file.path("data/tracks", sprintf("data_%s.rds", k))
  if (!file.exists(f)) {
    return(NULL)
  }
  v <- readRDS(f)$stp.var.nonp[[1]]
  v <- v[v$case_ == TRUE, ]
  out <- tibble::tibble(
    key = k,
    step_id_ = v$step_id_,
    wiscland_start_chr = as.character(v$wiscland_start)
  )
  for (nm in NUMERIC_START) {
    out[[nm]] <- v[[nm]]
  }
  out
}

for (season in SEASONS) {
  cat(sprintf("\n########## %s ##########\n", season))
  pooled <- xgb_season_pool(season)
  starts <- purrr::map_dfr(unique(pooled$key), start_of_key)
  d <- dplyr::left_join(pooled, starts, by = c("key", "step_id_"))
  d$wiscland_start <- as.integer(factor(d$wiscland_start_chr,
                                        levels = LANDCOVER_LEVELS)) - 1L
  d$wiscland_start_chr <- NULL

  candidates <- c("wiscland_start", NUMERIC_START)
  present <- vapply(d[, candidates], function(x) mean(!is.na(x)),
                    numeric(1))
  keep <- names(present)[present >= MIN_PRESENT]
  cat("start columns, share present:\n")
  print(round(present, 3))
  cat(sprintf("keeping: %s\n", paste(keep, collapse = ", ")))
  d <- d[, setdiff(names(d), setdiff(candidates, keep))]

  # A start column is constant within a stratum, so a missing one is
  # missing for the whole step. Those steps go, rather than leaving a
  # stratum with a column the trees would read as a category of its own.
  has_na <- Reduce(`|`, lapply(d[, keep, drop = FALSE], is.na))
  bad <- unique(d$stratum[has_na])
  if (length(bad)) {
    cat(sprintf("dropping %d of %d strata with a missing start value\n",
                length(bad), dplyr::n_distinct(d$stratum)))
    d <- d[!d$stratum %in% bad, ]
  }

  set.seed(SEED)
  sid <- as.integer(factor(d$stratum, levels = unique(d$stratum)))
  d$shadow_start <- stats::rnorm(max(sid))[sid]

  # Every added column must be constant within a stratum, or the argument
  # about where they can be split on stops holding.
  const_within <- function(x, stratum) {
    all(tapply(x, stratum, function(v) length(unique(v)) == 1))
  }
  for (nm in c(keep, "shadow_start")) {
    stopifnot(const_within(d[[nm]], d$stratum), !any(is.na(d[[nm]])))
  }

  cat(sprintf("final: %d rows, %d strata, %d deer-years, %d animals\n",
              nrow(d), dplyr::n_distinct(d$stratum),
              dplyr::n_distinct(d$key), dplyr::n_distinct(d$animal)))
  out <- sprintf("data/xgb/pooled_start_%s.rds", season)
  saveRDS(d, out)
  cat(sprintf("-> %s\n", out))
}
