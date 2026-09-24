#' @description
#' Run fit_season_xgb.R for every season, one after another.
#'
#' Each season is one pooled fit plus N_FOLDS more for the held-out scoring,
#' so the seasons run in sequence with all the threads rather than in
#' parallel with a share each. On the current cohort that is about an hour
#' for the three: nb is the largest (~54 min), fa next (~42 min), pf the
#' smallest (~11 min).
#'
#' pf 2020 and 2021 are excluded upstream — prep_pool_xgb.R never writes
#' them — so nothing here has to know about it.
#'
#' Inputs: data/xgb/pooled_<season>_<year>.rds (prep_pool_xgb.R)
#' Output: results/xgb/season_xgb_<season>.rds, one per season
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
SEASONS <- c("fa", "nb", "pf")

# Load packages ---------------------------------------------------------------
library(tidyverse)

script <- "scripts/xgb/fit_season_xgb.R"
lines <- readLines(script)
season_line <- grep('^SEASON <- ', lines)
stopifnot(length(season_line) == 1)

for (s in SEASONS) {
  cat(sprintf("\n########## season %s ##########\n", s))
  # fit_season_xgb.R is configured by the block at its top, so the season is
  # set by rewriting that one line into a temporary copy. The script on disk
  # is left alone.
  tmp <- tempfile(fileext = ".R")
  patched <- lines
  patched[season_line] <- sprintf('SEASON <- "%s"', s)
  writeLines(patched, tmp)
  status <- system2("Rscript", tmp)
  if (status != 0) {
    stop(sprintf("season %s failed with status %d", s, status))
  }
  unlink(tmp)
}

cat("\nall seasons done\n")
