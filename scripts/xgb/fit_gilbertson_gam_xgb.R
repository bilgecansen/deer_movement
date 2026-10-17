#' @description
#' Fit the Gilbertson et al. (2025) habitat model for every deer-year of
#' one set, each on its own: the second model the xgb log-score filter
#' scores its models against (logscore_gam_xgb.R), after the GAM null.
#'
#' Their final step-selection model, with distance to the home-range centre
#' added as in the GAM null and our own layers:
#'   log(sl_) + cos(ta_) + eastness_end + I(eastness_end^2) +
#'   northness_end + forest_edge_end + wiscland_end * ndvi_end +
#'   s(HR_center_end)
#' Landcover is our ten classes, forest the reference, and forest_edge_end
#' our signed distance (m; + inside forest, - outside), both as the xgb
#' models read them. Their neighbour and snow-depth terms have no
#' counterpart here. nb has no NDVI, so it drops ndvi_end and the landcover
#' by NDVI interaction.
#'
#' It is fitted as the GAM null is (fit_null_gam_xgb.R): a stratified Cox
#' model in mgcv by REML, on the same N_PTS of each step's 25 amt random
#' points (same seed), so the movement terms correct the deer-year's gamma
#' and von Mises.
#'
#' Landcover is kept at all ten levels, and its coefficients carry a weak
#' fixed ridge (mgcv paraPen): a normal prior with mean 0 and sd PRIOR_SD
#' on each class's effect relative to forest at the deer-year's mean NDVI,
#' and on each class's NDVI slope per sd of the deer-year's NDVI. NDVI is
#' centred at that mean inside the model, I(ndvi_end - mean), which leaves
#' the unpenalized model unchanged but puts the class effects where the
#' data are. Without the ridge a class with available points but no used
#' step has its coefficient run off towards minus infinity, and REML fails
#' or stops early on most deer-years. With it, such a class stops at a
#' finite value and a class the deer-year never has stays at forest's. The
#' per-class counts and coefficients (`classes`), REML's convergence and
#' any warnings are recorded per deer-year.
#'
#' Inputs:  data/tracks/data_<key>.rds
#' Output:  results/xgb/gilbertson_gam/gilbertson_gam_<key>.rds, one per
#'          deer-year: the fit (fit_gam_mod()), the rows it was fitted on,
#'          the tentative distributions, the per-class table, convergence
#'          and warnings, N_PTS and the seed
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
# The deer-years to fit: "test" or "train"
SET <- "train"
# Random points kept per step, of the 25 the wrangle drew; with SEED, the
# same points as the GAM null
N_PTS <- 10L
SEED <- 20261016L
# Prior sd of the landcover terms (log scale), see the description
PRIOR_SD <- 5
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
  # The ridge, one fixed penalty per landcover term
  n_lc <- length(LANDCOVER_LEVELS) - 1L
  pen <- list(wiscland_end = list(diag(n_lc), sp = 1 / PRIOR_SD^2))
  ndvi_centre <- NA_real_
  if (with_ndvi) {
    ndvi_centre <- signif(mean(d$ndvi_end), 6)
    ndvi_term <- sprintf("I(ndvi_end - %s)", ndvi_centre)
    lc_terms <- paste("wiscland_end *", ndvi_term)
    pen[[paste0("wiscland_end:", ndvi_term)]] <- list(
      diag(n_lc), sp = stats::sd(d$ndvi_end)^2 / PRIOR_SD^2
    )
  } else {
    lc_terms <- "wiscland_end"
  }
  formula <- paste(
    "log(sl_) + cos(ta_) + eastness_end + I(eastness_end^2) +",
    "northness_end + forest_edge_end +", lc_terms, "+ s(HR_center_end)"
  )
  warns <- character(0)
  res <- withCallingHandlers(
    fit_gam_mod(d, formula, select = FALSE, paraPen = pen),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (is.character(res$gam)) {
    stop("the fit failed")
  }
  b <- stats::coef(res$gam)
  se <- sqrt(diag(res$gam$Vp))
  lc <- LANDCOVER_LEVELS[-1]
  classes <- tibble::tibble(
    class = lc,
    used = as.integer(table(d$wiscland_end[d$case_ == 1])[lc]),
    available = as.integer(table(d$wiscland_end[d$case_ == 0])[lc]),
    coef = unname(b[paste0("wiscland_end", lc)]),
    se = unname(se[paste0("wiscland_end", lc)])
  ) |>
    dplyr::mutate(status = dplyr::case_when(
      used + available == 0 ~ "absent",
      used == 0 ~ "never_used",
      TRUE ~ "used"
    ))
  saveRDS(
    list(key = key, fit = res$gam, smooth_diag = res$smooth_diag,
         formula = formula, end_vars = end_vars, n_pts = N_PTS,
         seed = seed, prior_sd = PRIOR_SD, ndvi_centre = ndvi_centre,
         data = d, classes = classes,
         conv = res$gam$outer.info$conv, warnings = unique(warns),
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
cat("landcover classes per deer-year: absent (no rows), never used ",
    "(available, no used step), used\n", sep = "")
classes <- purrr::map_dfr(fits, function(x) {
  dplyr::mutate(x$classes, key = x$key)
})
print(as.data.frame(table(classes$class, classes$status)) |>
        tidyr::pivot_wider(names_from = Var2, values_from = Freq) |>
        dplyr::rename(class = Var1), row.names = FALSE)
cat(sprintf("deer-years with at least one never-used class: %d\n",
            dplyr::n_distinct(classes$key[classes$status == "never_used"])))
