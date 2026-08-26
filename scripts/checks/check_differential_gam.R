#' @description
#' Differential test on REAL deer: the amt and GAM paths are two independent
#' implementations of the same estimator, so where they are given the same data
#' and the same model structure they must agree.
#'
#' The synthetic version of this (check_synthetic_gam.R, Part 2) runs on a
#' generated world. This one runs on real tracks and real rasters, so it also
#' exercises the actual covariate values, the real cropped grids, and the fitted
#' tentative distributions.
#'
#' Two comparisons per deer, both with an identical LINEAR structure so the two
#' fits are literally the same likelihood (conditional logistic == stratified
#' Cox PH with one event per stratum):
#'   (a) coefficients from amt::fit_issf vs mgcv::gam(family = cox.ph)
#'   (b) the per-step probabilities gate 3 is built on -- each model's
#'       redistribution kernel normalised over the disc and read at the observed
#'       endpoint, which is what onestep_logscore*() does.
#'
#' Both kernels are given the SAME covariate function and the SAME observed
#' incoming heading, so the only thing that differs is the kernel implementation
#' itself -- amt's (its own Jacobian, normalisation, rasterisation) against
#' ours. A shared simple cov_fun is used rather than the two production ones,
#' precisely so that a disagreement can only come from the kernel maths.
#'
#' Usage: Rscript scripts/checks/check_differential_gam.R [n_deer] [n_steps]

suppressPackageStartupMessages({
  library(mgcv)
  library(amt)
  library(terra)
  library(tidyverse)
  library(sf)
  library(survival)
})
source("scripts/helper_functions.R")
source("scripts/checks/check_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
N_DEER <- if (length(args) >= 1) as.integer(args[1]) else 5L
N_STEPS <- if (length(args) >= 2) as.integer(args[2]) else 30L

# Spread the deer over seasons and track lengths rather than taking the first N.
all_keys <- gsub("^data_(.*)\\.rds$", "\\1",
                 basename(list.files("data/tracks", "^data_.*\\.rds$")))
all_keys <- all_keys[
  file.exists(sprintf("results/gam/results_gam_%s.rds", all_keys))
]
seasons <- vapply(strsplit(all_keys, "_"), `[`, character(1), 2)
set.seed(20)
keys <- unlist(lapply(split(all_keys, seasons), function(k) {
  sample(k, min(length(k), ceiling(N_DEER / dplyr::n_distinct(seasons)) + 1))
}))
keys <- head(unique(keys), N_DEER)

# Linear RHS. HR_center_end is present for every deer and carries real signal,
# so the comparison is not run on a degenerate near-zero coefficient.
RHS <- "HR_center_end + sl_ + log(sl_) + cos(ta_)"
TERMS <- c("HR_center_end", "sl_", "log(sl_)", "cos(ta_)")

# Shared covariate function: plain extraction, identical for both kernels.
cov_fun <- function(xy, map) amt::extract_covariates(xy, map, where = "both")

summary_rows <- list()

for (key in keys) {
  p <- strsplit(key, "_")[[1]]
  id <- p[1]
  season <- p[2]
  year <- as.integer(p[3])
  check_section(sprintf("deer %s", key))

  deer <- readRDS(sprintf("data/tracks/data_%s.rds", key))
  d <- deer$stp.var[[1]]
  sl_distr <- attr(d, "sl_")
  ta_distr <- attr(d, "ta_")
  stp <- deer$stp[[1]]

  d <- d[!is.na(d$HR_center_end) & !is.na(d$sl_) & !is.na(d$ta_) & d$sl_ > 0, ]
  attr(d, "sl_") <- sl_distr
  attr(d, "ta_") <- ta_distr
  attr(d, "crs") <- 6610
  class(d) <- unique(c("random_steps", "steps_xyt", "steps_xy", class(d)))

  # ---- (a) coefficients ----------------------------------------------------
  iss <- tryCatch(
    amt::fit_issf(
      d,
      stats::as.formula(paste("case_ ~", RHS, "+ strata(step_id_)")),
      model = TRUE
    ),
    error = function(e) NULL)
  gd <- prepare_gam_data(d)
  gm <- tryCatch(
    mgcv::gam(stats::as.formula(paste("cbind(times, stratum) ~", RHS)),
              data = gd, family = mgcv::cox.ph(), weights = obs,
              method = "REML"),
    error = function(e) NULL)

  if (is.null(iss) || is.null(gm)) {
    check(sprintf("%s: both models fit", key), NA, "a fit failed")
    next
  }
  ci <- coef(iss$model)[TERMS]
  cg <- coef(gm)[TERMS]
  max_cd <- max(abs(as.numeric(cg - ci)))
  cat(sprintf("  coefficients (amt | gam): %s\n",
      paste(sprintf("%s %.5f|%.5f", TERMS, ci, cg), collapse = "  ")))
  check(sprintf("%s: coefficients agree", key), max_cd < 1e-4,
        sprintf("max |diff| = %.2e", max_cd))

  # ---- (b) per-step kernel probabilities -----------------------------------
  env <- load_landcover(year, season)
  ce <- sf::st_buffer(sf::st_as_sf(stp, coords = c("x1_", "y1_"), crs = 6610),
                      CROP_BUFFER_M)
  envc <- terra::crop(env, ce)
  envc$HR_center <- load_hr_center_raster(id, season, year, envc)

  # Real incoming heading, as the corrected scoring path now uses.
  s <- stp |> dplyr::group_by(burst_) |>
    dplyr::mutate(ph = dplyr::lag(atan2(y2_ - y1_, x2_ - x1_))) |>
    dplyr::ungroup() |> dplyr::filter(!is.na(ph))
  set.seed(31)
  idx <- sort(sample(nrow(s), min(N_STEPS, nrow(s))))

  dens <- function(k, at) {
    e <- terra::extract(k, at)
    as.numeric(e[1, ncol(e)])
  }
  lp <- purrr::map_dfr(idx, function(i) {
    st <- amt::make_start(c(s$x1_[i], s$y1_[i]), ta_ = s$ph[i], time = s$t1_[i],
                          dt = lubridate::hours(4), crs = 6610)
    ki <- tryCatch(amt::redistribution_kernel(
      x = iss, map = envc, fun = cov_fun, start = st,
      landscape = "discrete", as.rast = TRUE)$redistribution.kernel,
      error = function(e) NULL)
    kg <- tryCatch(redistribution_kernel_gam(
      x = gm, map = envc, start = st, fun = cov_fun,
      sl_distr = sl_distr, ta_distr = ta_distr,
      compensate.movement = TRUE, as.rast = TRUE)$redistribution.kernel,
      error = function(e) NULL)
    if (is.null(ki) || is.null(kg)) return(NULL)
    at <- cbind(s$x2_[i], s$y2_[i])
    tibble(amt = log(dens(ki, at)), gam = log(dens(kg, at)))
  }) |> dplyr::filter(is.finite(amt), is.finite(gam))

  if (nrow(lp) < 5) {
    check(sprintf("%s: enough steps scored by both", key), NA,
          sprintf("only %d", nrow(lp)))
    next
  }
  dd <- lp$gam - lp$amt
  cat(sprintf("  per-step log p over %d steps: mean amt %.4f, mean gam %.4f\n",
              nrow(lp), mean(lp$amt), mean(lp$gam)))
  check(sprintf("%s: per-step probabilities agree", key), max(abs(dd)) < 1e-6,
        sprintf("max |diff| = %.2e, cor = %.8f",
                max(abs(dd)), cor(lp$amt, lp$gam)))

  summary_rows[[key]] <- tibble(
    key = key, season = season, n_strata = max(gd$stratum),
    max_coef_diff = max_cd, n_steps = nrow(lp),
    max_logp_diff = max(abs(dd))
  )
}

cat("\n=== summary ===\n")
print(as.data.frame(dplyr::bind_rows(summary_rows) |>
  dplyr::mutate(max_coef_diff = signif(max_coef_diff, 3),
                max_logp_diff = signif(max_logp_diff, 3))), row.names = FALSE)

quit(status = if (check_summary() > 0) 1L else 0L)
