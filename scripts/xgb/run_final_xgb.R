#' @description
#' Pick the model types to simulate paths from and fit each of them with
#' fit_final_xgb.R, several at a time.
#'
#' For every season, within each family (the types without deer, and their
#' deer twins), two model types:
#'   * best      the highest held-out score
#'   * simplest  the simplest type within THRESHOLD log units per deer of
#'               the best, the type plot_models_xgb.R selects
#' Simulating both shows what the simpler model gives up in prediction. A
#' type that is both is fitted once.
#'
#' The fits are mostly single-threaded R, so they run one thread each, as
#' many at a time as leaves the machine a core.
#'
#' Inputs: results/xgb/compare_<season>_<type>.rds (fit_model_xgb.R)
#' Output: results/xgb/final_<season>_<type>.rds, one per pick, and
#'         results/xgb/final_selection.rds, the picks and their roles
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
# Log units per deer, and the order of simplicity, as in plot_models_xgb.R
THRESHOLD <- 3
COMPLEXITY <- c("null", "main", "rsf", "rsf_hr", "start", "full", "full_hr")
N_PAR <- parallel::detectCores() - 1L
N_THREAD <- 1L
POLL_SECONDS <- 20

# Load packages ---------------------------------------------------------------
library(tidyverse)

files <- list.files("results/xgb", pattern = "^compare_.*[.]rds$",
                    full.names = TRUE)
scores <- purrr::map_dfr(files, function(f) {
  x <- readRDS(f)
  tibble(season = x$season, config = x$config, held_out = x$ll_cv,
         n_deer = x$n_deer_years, n_steps = x$n_steps)
})
scores <- scores |>
  mutate(deer = grepl("_deer$", config), type = sub("_deer$", "", config))
# Every family must carry every type, or a missing one could change a pick
stopifnot(all(count(scores, season, deer)$n == length(COMPLEXITY)),
          all(scores$type %in% COMPLEXITY))

picks <- scores |>
  group_by(season, deer) |>
  mutate(per_deer = (held_out - max(held_out)) / n_deer,
         best = held_out == max(held_out),
         simplest = type == COMPLEXITY[min(match(type[per_deer > -THRESHOLD],
                                                 COMPLEXITY))]) |>
  ungroup() |>
  filter(best | simplest) |>
  # Largest seasons first, so the longest fits do not start last
  arrange(desc(n_steps), deer, match(type, COMPLEXITY))
stopifnot(all(count(picks, season, deer, best)$n == 1),
          all(count(picks, season, deer, simplest)$n == 1))

cat("=== model types to fit ===\n")
print(as.data.frame(select(picks, season, config, per_deer, best, simplest)),
      row.names = FALSE, digits = 4)
dir.create("results/xgb", showWarnings = FALSE, recursive = TRUE)
saveRDS(select(picks, season, config, deer, type, per_deer, best, simplest),
        "results/xgb/final_selection.rds")

script <- "scripts/xgb/fit_final_xgb.R"
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

running <- function() {
  length(system2("pgrep", c("-f", shQuote("fit_final_xgb")),
                 stdout = TRUE, stderr = FALSE))
}

for (i in seq_len(nrow(picks))) {
  while (running() >= N_PAR) {
    Sys.sleep(POLL_SECONDS)
  }
  s <- picks$season[i]
  cf <- picks$config[i]
  # The season, model type and thread count are set by rewriting those
  # three lines into a temporary copy; the script on disk is left alone.
  tmp <- tempfile(pattern = sprintf("fit_final_xgb_%s_%s_", s, cf),
                  fileext = ".R")
  patched <- lines
  patched[i_season] <- sprintf('SEASON <- "%s"', s)
  patched[i_config] <- sprintf('CONFIG <- "%s"', cf)
  patched[i_thread] <- sprintf("N_THREAD <- %dL", N_THREAD)
  writeLines(patched, tmp)
  log <- sprintf("logs/xgb/final_%s_%s.log", s, cf)
  cat(sprintf("start %s %s  (%s)\n", s, cf, format(Sys.time(), "%H:%M")))
  system2("Rscript", tmp, stdout = log, stderr = log, wait = FALSE)
  Sys.sleep(POLL_SECONDS)
}

while (running() > 0) {
  Sys.sleep(POLL_SECONDS)
}
cat(sprintf("\nall final fits done (%s)\n", format(Sys.time(), "%H:%M")))
