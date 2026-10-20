#' @description
#' Range-residency screen after Gilbertson et al. (2025): ctmm's initial
#' estimate of the positional autocorrelation timescale, tau, for every
#' deer-season-year. tau is roughly the time the deer takes to cross its
#' range; a deer-year whose initial tau is over TAU_MAX_DAYS
#' (helpers_constants.R) is flagged as not range resident (a range shift, a
#' long trip away or a drift), and canonical_deer() leaves it out.
#'
#' The initial estimate is ctmm.guess(interactive = FALSE), the starting
#' values ctmm fits from, on the same fixes create_hr_rasters.R fits the
#' home range to: every step start plus the last fix of each burst,
#' training and test steps alike.
#'
#' Screened: every deer-year on the canonical list without the screen
#' (canonical_deer(residency = FALSE)), plus the deer-seasons removed by
#' hand (EXCLUDED_DEER_SEASONS) in the training and test years, with
#' on_list FALSE, for comparison.
#'
#' Input:  library/SW_filtered_deer.RData
#' Output: data/residency_tau.rds, one row per deer-year: key, id, season,
#'         year, set, on_list, n_fixes, span_days, tau_days, flagged
#'
#' Usage: Rscript scripts/residency_tau.R

# Load packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(ctmm)
  library(tidyverse)
  library(sf)
})

# helper functions
source("scripts/helper_functions.R")

# Deer-years ------------------------------------------------------------------
tracks <- readRDS("library/SW_filtered_deer.RData")
on_list <- canonical_deer(tracks, residency = FALSE)
on_list$on_list <- TRUE
year <- as.integer(tracks$year)
by_hand <- tracks[paste(tracks$id, tracks$season, sep = "_") %in%
                    EXCLUDED_DEER_SEASONS &
                    year %in% c(TRAIN_YEARS, TEST_YEAR), ]
by_hand$set <- ifelse(as.integer(by_hand$year) %in% TRAIN_YEARS, "train",
                      "test")
by_hand$on_list <- FALSE
deer <- dplyr::bind_rows(on_list, by_hand)
rm(tracks)

# Initial tau -----------------------------------------------------------------
# The fixes of one deer-year as a ctmm telemetry object, the way
# create_hr_rasters.R builds them
deer_telemetry <- function(stp) {
  locs_start <- dplyr::select(stp, x = x1_, y = y1_, timestamp = t1_)
  locs_end <- stp |>
    dplyr::group_by(burst_) |>
    dplyr::slice_tail(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(x = x2_, y = y2_, timestamp = t2_)
  dplyr::bind_rows(locs_start, locs_end) |>
    dplyr::arrange(timestamp) |>
    dplyr::distinct(timestamp, .keep_all = TRUE) |>
    sf::st_as_sf(coords = c("x", "y"), crs = 6610) |>
    sf::st_transform(4326) |>
    dplyr::mutate(longitude = sf::st_coordinates(geometry)[, 1],
                  latitude = sf::st_coordinates(geometry)[, 2]) |>
    sf::st_drop_geometry() |>
    as.data.frame() |>
    ctmm::as.telemetry()
}

# ctmm.guess returns tau in seconds, position first; a guess without a
# position timescale (an IID model) has tau 0
initial_tau_days <- function(tele) {
  guess <- ctmm::ctmm.guess(tele, interactive = FALSE)
  tau <- guess$tau
  if (length(tau) == 0) {
    return(0)
  }
  unname(tau[1]) / 86400
}

out <- purrr::map_dfr(seq_len(nrow(deer)), function(i) {
  tele <- suppressMessages(suppressWarnings(deer_telemetry(deer$stp[[i]])))
  tibble::tibble(
    key = sprintf("%s_%s_%d", deer$id[i], deer$season[i],
                  as.integer(deer$year[i])),
    id = deer$id[i],
    season = deer$season[i],
    year = as.integer(deer$year[i]),
    set = deer$set[i],
    on_list = deer$on_list[i],
    n_fixes = nrow(tele),
    span_days = as.numeric(diff(range(tele$t))) / 86400,
    tau_days = initial_tau_days(tele)
  )
})
out$flagged <- out$tau_days > TAU_MAX_DAYS

saveRDS(out, "data/residency_tau.rds")

# Summary ---------------------------------------------------------------------
cat("=== initial tau (days), deer-years on the canonical list ===\n")
print(as.data.frame(out |>
                      dplyr::filter(on_list) |>
                      dplyr::group_by(season, set) |>
                      dplyr::summarise(
                        deer_years = dplyr::n(),
                        median = stats::median(tau_days),
                        q90 = stats::quantile(tau_days, 0.9),
                        max = max(tau_days),
                        over_1 = sum(tau_days > 1),
                        flagged = sum(flagged),
                        .groups = "drop"
                      )),
      digits = 3, row.names = FALSE)
cat(sprintf("\n=== flagged: tau over %g days ===\n", TAU_MAX_DAYS))
print(as.data.frame(out |>
                      dplyr::filter(flagged) |>
                      dplyr::arrange(dplyr::desc(tau_days)) |>
                      dplyr::select(key, set, on_list, n_fixes, span_days,
                                    tau_days)),
      digits = 3, row.names = FALSE)
cat("\n=== deer-seasons removed by hand ===\n")
print(as.data.frame(out[!out$on_list, c("key", "set", "tau_days", "flagged")]),
      digits = 3, row.names = FALSE)
cat("-> data/residency_tau.rds\n")
