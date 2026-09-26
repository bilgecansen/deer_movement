#' @description
#' Run fit_shape_xgb.R for every shape and season, a few at a time.
#'
#' Twelve jobs — four shapes x three seasons — each one whole-data fit plus
#' N_FOLDS more for the held-out score. The jobs are independent, so they
#' run in parallel with a share of the threads each rather than in sequence
#' with all of them; xgboost scales sublinearly with threads, so that is
#' the faster way round. On the current cohort the whole set takes about
#' two and a half hours at N_PAR 3.
#'
#' fa and nb are interleaved so at most two of the large-season jobs hold
#' memory at once.
#'
#' Inputs: data/xgb/pooled_start_<season>.rds (prep_start_xgb.R)
#' Output: results/xgb/compare_<season>_<shape>.rds, one per cell
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
SEASONS <- c("fa", "nb", "pf")
CONFIGS <- c("full", "rsf", "start", "main")
# Jobs at a time, and threads each. N_PAR * N_THREAD should leave the
# machine a core.
N_PAR <- 3L
N_THREAD <- 3L
POLL_SECONDS <- 20

# Load packages ---------------------------------------------------------------
library(tidyverse)

script <- "scripts/xgb/fit_shape_xgb.R"
lines <- readLines(script)
line_of <- function(pattern) {
  i <- grep(pattern, lines)
  stopifnot(length(i) == 1)
  i
}
i_season <- line_of("^SEASON <- ")
i_config <- line_of("^CONFIG <- ")
i_thread <- line_of("^N_THREAD <- ")

dir.create("logs/xgb", showWarnings = FALSE, recursive = TRUE)

# Interleave the seasons so the two largest do not all land together.
jobs <- tidyr::expand_grid(config = CONFIGS, season = SEASONS) |>
  dplyr::arrange(match(config, CONFIGS), match(season, SEASONS))

running <- function() {
  length(system2("pgrep", c("-f", shQuote("fit_shape_xgb")),
                 stdout = TRUE, stderr = FALSE))
}

for (i in seq_len(nrow(jobs))) {
  while (running() >= N_PAR) {
    Sys.sleep(POLL_SECONDS)
  }
  s <- jobs$season[i]
  cf <- jobs$config[i]
  # fit_shape_xgb.R is configured by the block at its top, so the season,
  # the shape and the thread count are set by rewriting those three lines
  # into a temporary copy. The script on disk is left alone.
  tmp <- tempfile(pattern = sprintf("fit_shape_xgb_%s_%s_", s, cf),
                  fileext = ".R")
  patched <- lines
  patched[i_season] <- sprintf('SEASON <- "%s"', s)
  patched[i_config] <- sprintf('CONFIG <- "%s"', cf)
  patched[i_thread] <- sprintf("N_THREAD <- %dL", N_THREAD)
  writeLines(patched, tmp)
  log <- sprintf("logs/xgb/compare_%s_%s.log", s, cf)
  cat(sprintf("start %s %s  (%s)\n", s, cf, format(Sys.time(), "%H:%M")))
  system2("Rscript", tmp, stdout = log, stderr = log, wait = FALSE)
  Sys.sleep(POLL_SECONDS)
}

while (running() > 0) {
  Sys.sleep(POLL_SECONDS)
}
cat(sprintf("\nall shapes done (%s)\n", format(Sys.time(), "%H:%M")))
