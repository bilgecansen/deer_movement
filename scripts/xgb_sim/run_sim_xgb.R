#' @description
#' Run the framework on simulated data. For every scenario and replicate,
#' draw which endpoint of each step was used from the true model
#' (sim_truth_xgb.R), then fit the model types with fit_model_xgb.R exactly
#' as on the real data: same steps, columns, folds, rounds and settings.
#' Only the choices differ, and they are known to follow the truth.
#'
#' Jobs run a few at a time, replicate by replicate, so a partial run
#' already holds complete replicates of every scenario. A job whose result
#' exists is skipped, so an interrupted run resumes where it stopped.
#'
#' fit_model_xgb.R is configured by rewriting lines of a temporary copy, as
#' run_models_xgb.R does: the season, type and threads, where the table
#' comes from (the real table with the simulated choices), where the result
#' goes, and one addition — the whole-data boosters are kept, so the fitted
#' curves can be compared with the true ones afterwards.
#'
#' Inputs: sims/xgb/<season>/truth.rds (sim_truth_xgb.R),
#'         data/xgb/pooled_start_<season>.rds
#' Output: sims/xgb/<season>/<scenario>/rep_<k>/choice.rds and
#'         compare_<season>_<type>.rds
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
SEASON <- "pf"
SCENARIOS <- c("A", "B", "C", "D")
N_REP <- 10L
CONFIGS <- c("null", "main", "rsf")
SEED_BASE <- 1000L
# Most of a fit's time is single-threaded R, so many one-thread jobs use the
# machine better than a few multi-thread ones. Each job needs about 3.5 GB.
N_PAR <- 7L
N_THREAD <- 1L
POLL_SECONDS <- 20
OUT_ROOT <- "sims/xgb"
# NA runs the production rounds; a small number is for testing this script
TEST_ROUNDS <- NA

# Load packages ---------------------------------------------------------------
library(tidyverse)

# helper functions
source("scripts/helper_functions.R")

season_dir <- file.path(OUT_ROOT, SEASON)
pool_path <- sprintf("data/xgb/pooled_start_%s.rds", SEASON)
truth <- readRDS(file.path("sims/xgb", SEASON, "truth.rds"))

# Draw the choices -------------------------------------------------------------
sizes <- xgb_strata_sizes(readRDS(pool_path))
for (s in seq_along(SCENARIOS)) {
  eta <- truth$eta_null + truth$scenarios[[SCENARIOS[s]]]$score
  for (k in seq_len(N_REP)) {
    rep_dir <- file.path(season_dir, SCENARIOS[s], sprintf("rep_%02d", k))
    dir.create(rep_dir, showWarnings = FALSE, recursive = TRUE)
    choice_path <- file.path(rep_dir, "choice.rds")
    if (!file.exists(choice_path)) {
      saveRDS(xgb_sim_choose(eta, sizes, SEED_BASE + 100L * s + k),
              choice_path)
    }
  }
}

# Fit --------------------------------------------------------------------------
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
i_rounds <- line_of("^N_ROUNDS <- ")
i_out <- line_of("^out_path <- ")
i_read <- line_of("^d <- readRDS\\(")
i_save <- line_of("n_trees = length\\(specs\\) \\* N_ROUNDS\\)")

dir.create("logs/xgb_sim", showWarnings = FALSE, recursive = TRUE)
jobs <- tidyr::expand_grid(rep = seq_len(N_REP), scenario = SCENARIOS,
                           config = CONFIGS)

# Number of fits running whose temporary script name contains `pattern`.
# pgrep exits with status 1 when nothing matches, which is not an error here.
running <- function(pattern = "fit_sim_xgb_") {
  length(suppressWarnings(system2("pgrep", c("-f", shQuote(pattern)),
                                  stdout = TRUE, stderr = FALSE)))
}

for (i in seq_len(nrow(jobs))) {
  rep_dir <- file.path(season_dir, jobs$scenario[i],
                       sprintf("rep_%02d", jobs$rep[i]))
  out <- file.path(rep_dir, sprintf("compare_%s_%s.rds", SEASON,
                                    jobs$config[i]))
  job_name <- sprintf("fit_sim_xgb_%s_%02d_%s_", jobs$scenario[i],
                      jobs$rep[i], jobs$config[i])
  # Skip a fit that is finished, or still running from an earlier start of
  # this script
  if (file.exists(out) || running(job_name) > 0) {
    next
  }
  while (running() >= N_PAR) {
    Sys.sleep(POLL_SECONDS)
  }
  patched <- lines
  patched[i_season] <- sprintf('SEASON <- "%s"', SEASON)
  patched[i_config] <- sprintf('CONFIG <- "%s"', jobs$config[i])
  patched[i_thread] <- sprintf("N_THREAD <- %dL", N_THREAD)
  if (!is.na(TEST_ROUNDS)) {
    patched[i_rounds] <- sprintf("N_ROUNDS <- %dL", TEST_ROUNDS)
  }
  patched[i_out] <- sprintf('out_path <- "%s"', out)
  patched[i_read] <- sprintf('d <- xgb_sim_pool("%s", "%s")', pool_path,
                             file.path(rep_dir, "choice.rds"))
  patched[i_save] <- sub(
    "n_trees = length\\(specs\\) \\* N_ROUNDS\\)",
    "n_trees = length(specs) * N_ROUNDS, boosters = fit$boosters)",
    patched[i_save]
  )
  tmp <- tempfile(pattern = job_name, fileext = ".R")
  writeLines(patched, tmp)
  log <- sprintf("logs/xgb_sim/%s_%s_rep%02d_%s.log", SEASON,
                 jobs$scenario[i], jobs$rep[i], jobs$config[i])
  cat(sprintf("start %s rep %02d %s  (%s)\n", jobs$scenario[i], jobs$rep[i],
              jobs$config[i], format(Sys.time(), "%H:%M")))
  system2("Rscript", tmp, stdout = log, stderr = log, wait = FALSE)
  Sys.sleep(POLL_SECONDS)
}

while (running() > 0) {
  Sys.sleep(POLL_SECONDS)
}
cat(sprintf("\nall simulation fits done (%s)\n", format(Sys.time(), "%H:%M")))
