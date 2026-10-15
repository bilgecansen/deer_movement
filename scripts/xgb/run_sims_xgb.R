#' @description
#' Run sim_paths_xgb.R for every deer-year of one set, several at a time:
#' the test year, or the training years. Simulating the training years
#' shows how well the models reproduce the paths they were fitted on, which
#' separates a model that fits badly from one that fits well but transfers
#' badly to a new year.
#'
#' One job per deer-year, each simulating every final model of its season
#' (results/xgb/final_selection.rds). A job's cost grows with its number of
#' steps and with the area within reach of a step, so the costliest jobs
#' start first and the run does not end waiting on one of them. The jobs are
#' mostly single-threaded R, so they run one thread each, as many at a time
#' as leaves the machine a core.
#'
#' SKIP_EXISTING runs only the deer-years with no output yet, to resume an
#' interrupted run.
#'
#' Inputs:  as sim_paths_xgb.R
#' Output:  sims/xgb_paths/sims_xgb_<key>.rds, one per test deer-year, or
#'          sims/xgb_paths_train/sims_xgb_<key>.rds, one per training
#'          deer-year, and logs/xgb/sims_<key>.log
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
# The deer-years to simulate: "test" or "train"
SET <- "test"
SKIP_EXISTING <- TRUE
N_PAR <- parallel::detectCores() - 1L
N_THREAD <- 1L
POLL_SECONDS <- 20

# Load packages ---------------------------------------------------------------
library(tidyverse)

# helper functions
source("scripts/helper_functions.R")

deer <- canonical_deer(readRDS("library/SW_filtered_deer.RData"))
longest <- deer |>
  filter(set == "train") |>
  mutate(R = purrr::map_dbl(stp, function(s) max(s$sl_[!is.na(s$ta_)]))) |>
  group_by(id, season) |>
  summarise(reach = max(R), .groups = "drop")
jobs <- deer |>
  filter(set == SET) |>
  mutate(key = sprintf("%s_%s_%d", id, season, as.integer(year)),
         hours = purrr::map_dbl(stp, function(s) {
           as.numeric(difftime(max(s$t2_), min(s$t1_), units = "hours"))
         })) |>
  select(key, id, season, hours) |>
  left_join(longest, by = c("id", "season")) |>
  # Steps times the cells within reach, the kernel's work
  arrange(desc(hours * reach^2))
stopifnot(!anyNA(jobs$reach))
out_dir <- if (SET == "test") "sims/xgb_paths" else "sims/xgb_paths_train"
if (SKIP_EXISTING) {
  done <- file.exists(file.path(out_dir, sprintf("sims_xgb_%s.rds", jobs$key)))
  cat(sprintf("skipping %d deer-years that already have paths\n", sum(done)))
  jobs <- jobs[!done, ]
}

script <- "scripts/xgb/sim_paths_xgb.R"
lines <- readLines(script)
line_of <- function(pattern) {
  i <- grep(pattern, lines)
  stopifnot(length(i) == 1)
  i
}
i_key <- line_of("^KEY <- ")
i_thread <- line_of("^N_THREAD <- ")

dir.create("logs/xgb", showWarnings = FALSE, recursive = TRUE)

running <- function() {
  length(system2("pgrep", c("-f", shQuote("sim_paths_xgb_")),
                 stdout = TRUE, stderr = FALSE))
}

for (i in seq_len(nrow(jobs))) {
  while (running() >= N_PAR) {
    Sys.sleep(POLL_SECONDS)
  }
  k <- jobs$key[i]
  # The deer-year and thread count are set by rewriting those two lines
  # into a temporary copy; the script on disk is left alone.
  tmp <- tempfile(pattern = sprintf("sim_paths_xgb_%s_", k), fileext = ".R")
  patched <- lines
  patched[i_key] <- sprintf('KEY <- "%s"', k)
  patched[i_thread] <- sprintf("N_THREAD <- %dL", N_THREAD)
  writeLines(patched, tmp)
  log <- sprintf("logs/xgb/sims_%s.log", k)
  cat(sprintf("start %s  (%d of %d, %s)\n", k, i, nrow(jobs),
              format(Sys.time(), "%H:%M")))
  system2("Rscript", tmp, stdout = log, stderr = log, wait = FALSE)
  Sys.sleep(2)
}

while (running() > 0) {
  Sys.sleep(POLL_SECONDS)
}
cat(sprintf("\nall paths done (%s)\n", format(Sys.time(), "%H:%M")))
