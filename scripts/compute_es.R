#' @description
#' Compute one-step energy scores for every deer with simulations and every
#' model that was simulated. Writes a single combined data frame to
#' filters/es.rds.
#'
#' By default, deer that already appear in the existing filters/es.rds are
#' skipped (resumable reruns). Pass --overwrite to reprocess everything.
#'
#' Usage:
#'   Rscript scripts/compute_es.R              # resumable (default)
#'   Rscript scripts/compute_es.R --overwrite  # reprocess every deer

# Parse CLI args ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args

# Load packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(scoringRules)
})

source("scripts/helper_functions.R")

# Setup ------------------------------------------------------------------------
dir.create("filters", showWarnings = FALSE)
out_path <- "filters/es.rds"

existing <- if (!overwrite && file.exists(out_path)) {
  readRDS(out_path)
} else {
  tibble::tibble(
    id = character(),
    season = character(),
    year = integer(),
    model = integer(),
    energy_score = numeric()
  )
}

# Discover deer with simulations
sim_files <- list.files("sims", pattern = "^sims_.*\\.rds$", full.names = TRUE)
keys      <- gsub("^sims_(.*)\\.rds$", "\\1", basename(sim_files))

cat(sprintf("Found %d deer with simulations\n", length(sim_files)))

if (overwrite) {
  cat("[--overwrite] Recomputing ES for every deer.\n")
}

# Loop ------------------------------------------------------------------------
start_time <- Sys.time()
n_done    <- 0L
n_skipped <- 0L
n_failed  <- 0L

acc <- list()

for (i in seq_along(sim_files)) {
  key <- keys[i]
  parts  <- strsplit(key, "_")[[1]]
  id     <- parts[1]
  season <- parts[2]
  year   <- as.integer(parts[3])

  # Skip if this deer's rows already exist in `existing` (and not overwriting)
  if (!overwrite &&
    nrow(existing) > 0 &&
    any(
      existing$id == id &
        existing$season == season &
        existing$year == year
    )) {
    cat(sprintf("[skip-exists] %s\n", key))
    n_skipped <- n_skipped + 1L
    next
  }

  cat(sprintf("[run]         %s\n", key))

  one <- tryCatch(
    {
      results_sim <- readRDS(sim_files[i])
      one_deer    <- readRDS(sprintf("data/tracks/data_%s.rds", key))
      obs         <- one_deer$stp[[1]]

      per_model <- purrr::map_dfr(names(results_sim), function(m_chr) {
        sim_m <- results_sim[[m_chr]]
        if (length(sim_m) == 1 && is.na(sim_m)) {
          return(tibble::tibble(
            model = as.integer(m_chr),
            energy_score = NA_real_
          ))
        }
        tibble::tibble(
          model        = as.integer(m_chr),
          energy_score = calc_energy_score(obs = obs, sim = sim_m)
        )
      })

      per_model |>
        dplyr::mutate(
          id     = id,
          season = season,
          year   = year,
          .before = 1
        )
    },
    error = function(e) {
      cat(sprintf("[fail]        %s: %s\n", key, conditionMessage(e)))
      NULL
    }
  )

  if (is.null(one)) {
    n_failed <- n_failed + 1L
  } else {
    acc[[length(acc) + 1L]] <- one
    n_done <- n_done + 1L
  }
}

new_results <- dplyr::bind_rows(acc)

# When overwriting, drop rows for any deer we just recomputed before binding
if (overwrite && nrow(new_results) > 0) {
  combined <- dplyr::anti_join(
    existing,
    new_results |> dplyr::distinct(id, season, year),
    by = c("id", "season", "year")
  ) |>
    dplyr::bind_rows(new_results)
} else {
  combined <- dplyr::bind_rows(existing, new_results)
}

saveRDS(combined, out_path)

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf(
  "Done: %d   Skipped: %d   Failed: %d   Elapsed: %.1f min\n",
  n_done,
  n_skipped,
  n_failed,
  elapsed
))
cat(sprintf(
  "Wrote %s with %d total rows (%d new).\n",
  out_path,
  nrow(combined),
  nrow(new_results)
))
