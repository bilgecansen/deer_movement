#' @description
#' Run fit_model_xgb.R for every model type and season, a few at a time.
#'
#' Eighteen jobs — six model types x three seasons — each one whole-data fit
#' plus N_FOLDS more for the held-out score, with every fold replayed on its
#' held-out steps for the split gain. The jobs are independent, so they run
#' in parallel with a share of the threads each rather than in sequence
#' with all of them; xgboost scales sublinearly with threads, so that is
#' the faster way round.
#'
#' fa and nb are interleaved so at most two of the large-season jobs hold
#' memory at once.
#'
#' SKIP_EXISTING fits only the cells with no result yet, for adding a model
#' type. It is safe only while the season tables are unchanged: every type
#' is compared on the same folds of the same table, and a result fit to an
#' older table would no longer be comparable. After rebuilding a table,
#' refit everything.
#'
#' Inputs: data/xgb/pooled_start_<season>.rds (prep_start_xgb.R)
#' Output: results/xgb/compare_<season>_<type>.rds, one per cell
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
SEASONS <- c("fa", "nb", "pf")
CONFIGS <- c("full", "rsf", "start", "main", "rsf_hr", "null")
# TRUE: fit only the cells with no result file yet (see above)
SKIP_EXISTING <- FALSE
# Jobs at a time, and threads each. N_PAR * N_THREAD should leave the
# machine a core.
N_PAR <- 3L
N_THREAD <- 3L
POLL_SECONDS <- 20

# Load packages ---------------------------------------------------------------
library(tidyverse)

script <- "scripts/xgb/fit_model_xgb.R"
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
if (SKIP_EXISTING) {
  done <- file.exists(sprintf("results/xgb/compare_%s_%s.rds", jobs$season,
                              jobs$config))
  cat(sprintf("skipping %d cells that already have a result\n", sum(done)))
  jobs <- jobs[!done, ]
}

running <- function() {
  length(system2("pgrep", c("-f", shQuote("fit_model_xgb")),
                 stdout = TRUE, stderr = FALSE))
}

for (i in seq_len(nrow(jobs))) {
  while (running() >= N_PAR) {
    Sys.sleep(POLL_SECONDS)
  }
  s <- jobs$season[i]
  cf <- jobs$config[i]
  # fit_model_xgb.R is configured by the block at its top, so the season,
  # the model type and the thread count are set by rewriting those three
  # lines into a temporary copy. The script on disk is left alone.
  tmp <- tempfile(pattern = sprintf("fit_model_xgb_%s_%s_", s, cf),
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
cat(sprintf("\nall models done (%s)\n", format(Sys.time(), "%H:%M")))
