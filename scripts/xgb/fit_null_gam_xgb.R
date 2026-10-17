#' @description
#' Fit the GAM null for every deer-year of one set, each on its own: the
#' null the xgb log-score filter scores its models against
#' (logscore_gam_xgb.R).
#'
#' The null is the GAM path's (fit_GAM.R): movement plus distance to the
#' home-range centre,
#'   log(sl_) + cos(ta_) + s(tod_, bs = "cc", by = sl_) + s(HR_center_end)
#' fitted as a stratified Cox model in mgcv (cox.ph, Klappstein et al. 2024)
#' by REML, without the double penalty (select = FALSE), as there. The
#' random points are the wrangle's amt design (stp.var), drawn from the
#' deer-year's own gamma and von Mises, N_PTS of the 25 kept per step at
#' random (xgb_gam_data()). The movement terms therefore correct those
#' distributions. The time-of-day basis is capped below the deer-year's
#' distinct start times, as in fit_GAM.R.
#'
#' Each deer-year's points are drawn from its own seed, SEED plus its place
#' in the set's sorted keys, so a rerun or a partial run draws the same ones.
#'
#' Inputs:  data/tracks/data_<key>.rds
#' Output:  results/xgb/null_gam/null_gam_<key>.rds, one per deer-year: the
#'          fit (fit_gam_mod()), the rows it was fitted on, the tentative
#'          distributions, k_tod, N_PTS and the seed
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
# The deer-years to fit: "test" or "train"
SET <- "train"
# Random points kept per step, of the 25 the wrangle drew
N_PTS <- 10L
K_TOD <- 10L
SEED <- 20261016L
# Deer-years to fit: NULL for every one of the set not yet fitted
KEYS <- NULL
N_PAR <- parallel::detectCores() - 1L

# Load packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(furrr)
})

# helper functions
source("scripts/helper_functions.R")

out_dir <- "results/xgb/null_gam"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
deer <- canonical_deer(readRDS("library/SW_filtered_deer.RData"))
deer <- deer[deer$set == SET, ]
all_keys <- sort(sprintf("%s_%s_%d", deer$id, deer$season,
                         as.integer(deer$year)))
if (is.null(KEYS)) {
  KEYS <- all_keys[!file.exists(file.path(out_dir,
                                          sprintf("null_gam_%s.rds",
                                                  all_keys)))]
}
cat(sprintf("fitting %d %s deer-years on %d workers (%s)\n", length(KEYS),
            SET, N_PAR, format(Sys.time(), "%H:%M")))

fit_one <- function(key) {
  source("scripts/helper_functions.R")
  track <- readRDS(sprintf("data/tracks/data_%s.rds", key))
  stopifnot(track$set == SET)
  seed <- SEED + match(key, all_keys)
  d <- xgb_gam_data(track, N_PTS, seed)
  k_tod <- max(3L, min(K_TOD, length(unique(d$tod_)) - 1L))
  formula <- sprintf(
    paste0("log(sl_) + cos(ta_) + s(tod_, bs = 'cc', k = %d, by = sl_) + ",
           "s(HR_center_end)"),
    k_tod
  )
  res <- fit_gam_mod(d, formula, select = FALSE)
  if (is.character(res$gam)) {
    stop("the fit failed")
  }
  saveRDS(
    list(key = key, fit = res$gam, smooth_diag = res$smooth_diag,
         formula = formula, k_tod = k_tod, n_pts = N_PTS, seed = seed,
         data = d, sl_distr = attr(track$stp.var[[1]], "sl_"),
         ta_distr = attr(track$stp.var[[1]], "ta_")),
    file.path(out_dir, sprintf("null_gam_%s.rds", key))
  )
  "ok"
}

if (length(KEYS)) {
  future::plan(future::multisession, workers = min(N_PAR, length(KEYS)))
  status <- furrr::future_map_chr(
    KEYS,
    function(key) {
      tryCatch(fit_one(key),
               error = function(e) paste("failed:", conditionMessage(e)))
    },
    .options = furrr::furrr_options(seed = TRUE, stdout = FALSE,
                                    chunk_size = 1L)
  )
  future::plan(future::sequential)
  failed <- status != "ok"
  if (any(failed)) {
    cat(sprintf("%d deer-years failed:\n", sum(failed)))
    cat(sprintf("  %s: %s\n", KEYS[failed], status[failed]), sep = "")
  }
}

files <- list.files(out_dir, pattern = "^null_gam_.*[.]rds$",
                    full.names = TRUE)
# A cox.ph fit records its REML convergence in outer.info, not converged
conv <- vapply(files, function(f) {
  identical(readRDS(f)$fit$outer.info$conv, "full convergence")
}, logical(1))
cat(sprintf("-> %s: %d deer-years fitted, %d not fully converged (%s)\n",
            out_dir, length(files), sum(!conv),
            format(Sys.time(), "%H:%M")))
if (any(!conv)) {
  cat(sprintf("  %s\n", basename(files[!conv])), sep = "")
}
