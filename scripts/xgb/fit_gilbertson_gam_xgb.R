#' @description
#' Fit the Gilbertson et al. (2025) habitat model for every deer-year of
#' one set, each on its own: the second model the xgb log-score filter
#' scores its models against (logscore_gam_xgb.R), after the GAM null.
#'
#' Their final step-selection model's habitat terms, with the GAM null's
#' movement and its distance to the home-range centre, on our own layers:
#'   log(sl_) + cos(ta_) + s(tod_, bs = "cc", by = sl_) +
#'   eastness_end + I(eastness_end^2) + northness_end +
#'   I(log(abs(forest_edge_end))) + ndvi_c + s(wiscland_end, bs = "re") +
#'   s(wiscland_end, ndvi_c, bs = "re") + s(HR_center_end)
#' Landcover is our ten classes, and forest edge their log distance (m) to
#' the nearest edge, unsigned: the log of the absolute value of our signed
#' forest_edge_end (+ inside forest, - outside), never under 15 m, as the
#' xgb models read it. Their neighbour and snow-depth terms have no
#' counterpart here.
#'
#' Their landcover by NDVI interaction is fitted hierarchically within the
#' deer-year: every class has its own effect and its own NDVI slope, drawn
#' around a slope shared by all classes (ndvi_c) and a common effect, with
#' how far they spread estimated by REML from the deer-year's own classes.
#' A class the deer-year used rarely is pulled towards the shared values,
#' one it used often keeps its own, and one it never had sits at them. NDVI
#' enters as ndvi_c, the deer-year's NDVI less its mean over the rows
#' fitted (`ndvi_centre`, stored on the fit), so a class's effect is its
#' effect at the deer-year's typical NDVI. nb has no NDVI, so it keeps only
#' the class effects.
#'
#' It is fitted as the GAM null is (fit_null_gam_xgb.R): a stratified Cox
#' model in mgcv by REML, on the same N_PTS of each step's 25 amt random
#' points (same seed), so the movement terms correct the deer-year's gamma
#' and von Mises, with the null's cap on the time-of-day basis.
#'
#' Inputs:  data/tracks/data_<key>.rds
#' Output:  results/xgb/gilbertson_gam/gilbertson_gam_<key>.rds, one per
#'          deer-year: the fit (fit_gam_mod(), with `ndvi_centre`), the rows
#'          it was fitted on, the tentative distributions, the per-class
#'          table (counts, and each class's effect and NDVI slope relative to
#'          forest), convergence and warnings, N_PTS and the seed
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
# The deer-years to fit: "test" or "train"
SET <- "train"
# Random points kept per step, of the 25 the wrangle drew; with SEED, the
# same points as the GAM null
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

out_dir <- "results/xgb/gilbertson_gam"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
deer <- canonical_deer(readRDS("library/SW_filtered_deer.RData"))
deer <- deer[deer$set == SET, ]
all_keys <- sort(sprintf("%s_%s_%d", deer$id, deer$season,
                         as.integer(deer$year)))
if (is.null(KEYS)) {
  KEYS <- all_keys[!file.exists(file.path(
    out_dir, sprintf("gilbertson_gam_%s.rds", all_keys)
  ))]
}
cat(sprintf("fitting %d %s deer-years on %d workers (%s)\n", length(KEYS),
            SET, N_PAR, format(Sys.time(), "%H:%M")))

fit_one <- function(key) {
  source("scripts/helper_functions.R")
  track <- readRDS(sprintf("data/tracks/data_%s.rds", key))
  stopifnot(track$set == SET)
  seed <- SEED + match(key, all_keys)
  with_ndvi <- track$season != "nb"
  end_vars <- c("HR_center_end", if (with_ndvi) "ndvi_end",
                "forest_edge_end", "northness_end", "eastness_end",
                "wiscland_end")
  d <- xgb_gam_data(track, N_PTS, seed, end_vars)
  k_tod <- max(3L, min(K_TOD, length(unique(d$tod_)) - 1L))
  ndvi_centre <- NA_real_
  if (with_ndvi) {
    ndvi_centre <- signif(mean(d$ndvi_end), 6)
    d$ndvi_c <- d$ndvi_end - ndvi_centre
    lc_terms <- paste("ndvi_c + s(wiscland_end, bs = 're') +",
                      "s(wiscland_end, ndvi_c, bs = 're')")
  } else {
    lc_terms <- "s(wiscland_end, bs = 're')"
  }
  formula <- paste(
    sprintf("log(sl_) + cos(ta_) + s(tod_, bs = 'cc', k = %d, by = sl_) +",
            k_tod),
    "eastness_end + I(eastness_end^2) + northness_end +",
    "I(log(abs(forest_edge_end))) +",
    lc_terms, "+ s(HR_center_end)"
  )
  warns <- character(0)
  res <- withCallingHandlers(
    fit_gam_mod(d, formula, select = FALSE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (is.character(res$gam)) {
    stop("the fit failed")
  }
  fit <- res$gam
  # The kernel helpers read the centre from the fit (xgb_gam_end_columns())
  fit$ndvi_centre <- ndvi_centre

  # Each class's effect and NDVI slope relative to forest
  b <- stats::coef(fit)
  effect <- b[grep("^s\\(wiscland_end\\)\\.", names(b))]
  names(effect) <- LANDCOVER_LEVELS
  slope <- if (with_ndvi) {
    v <- b[grep("^s\\(wiscland_end,ndvi_c\\)\\.", names(b))]
    names(v) <- LANDCOVER_LEVELS
    v
  } else {
    stats::setNames(rep(NA_real_, length(LANDCOVER_LEVELS)), LANDCOVER_LEVELS)
  }
  lc <- LANDCOVER_LEVELS[-1]
  classes <- tibble::tibble(
    class = lc,
    used = as.integer(table(d$wiscland_end[d$case_ == 1])[lc]),
    available = as.integer(table(d$wiscland_end[d$case_ == 0])[lc]),
    effect = unname(effect[lc] - effect[["forest"]]),
    ndvi_slope = unname(slope[lc] - slope[["forest"]])
  ) |>
    dplyr::mutate(status = dplyr::case_when(
      used + available == 0 ~ "absent",
      used == 0 ~ "never_used",
      TRUE ~ "used"
    ))
  saveRDS(
    list(key = key, fit = fit, smooth_diag = res$smooth_diag,
         formula = formula, end_vars = end_vars, k_tod = k_tod,
         n_pts = N_PTS, seed = seed, ndvi_centre = ndvi_centre, data = d,
         classes = classes, conv = fit$outer.info$conv,
         warnings = unique(warns),
         sl_distr = attr(track$stp.var[[1]], "sl_"),
         ta_distr = attr(track$stp.var[[1]], "ta_")),
    file.path(out_dir, sprintf("gilbertson_gam_%s.rds", key))
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

files <- list.files(out_dir, pattern = "^gilbertson_gam_.*[.]rds$",
                    full.names = TRUE)
fits <- lapply(files, function(f) {
  x <- readRDS(f)
  x[c("key", "conv", "classes", "warnings")]
})
cat(sprintf("-> %s: %d deer-years fitted (%s)\n", out_dir, length(files),
            format(Sys.time(), "%H:%M")))
cat("REML convergence:\n")
print(table(vapply(fits, `[[`, character(1), "conv")))
classes <- purrr::map_dfr(fits, function(x) {
  dplyr::mutate(x$classes, key = x$key)
})
cat("largest class effect and NDVI slope relative to forest, per deer-year:\n")
print(summary(classes |>
                dplyr::group_by(key) |>
                dplyr::summarise(effect = max(abs(effect)),
                                 slope = suppressWarnings(
                                   max(abs(ndvi_slope), na.rm = TRUE)
                                 ),
                                 .groups = "drop") |>
                dplyr::select(effect, slope)))
