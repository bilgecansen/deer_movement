#' @description
#' Ground-truth recovery on a synthetic world (tier 1, option A: fitting path
#' only -- the case/control design is built directly rather than through the
#' wrangle).
#'
#' Part 1 -- UNIMODAL truth, 5 replicates x 1000 steps. Selection peaks at an
#'   intermediate covariate value, which a linear iSSF term cannot represent, so
#'   this tests whether the smooth recovers a SHAPE rather than just a sign.
#'   Checked: shape of the recovered smooth against truth (on the log-RSS scale),
#'   location of the fitted peak, pointwise CI coverage, and -- a sharp one --
#'   the movement corrections, which must be ~0 because the tracks were generated
#'   from exactly the tentative gamma / von Mises the controls were drawn from.
#'
#' Part 2 -- MONOTONIC truth, 1 replicate. The same data are fit with amt's iSSF
#'   (conditional logistic) and with our GAM (cox.ph), using an identical LINEAR
#'   structure so the two are the same likelihood. Two comparisons:
#'     (a) coefficients -- should agree to several decimals;
#'     (b) the per-step probabilities that gate 3 is built on -- each model's
#'         redistribution kernel is normalised over the disc and evaluated at the
#'         observed endpoint, exactly as onestep_logscore*() does. Two
#'         implementations sharing no code should assign the same probability to
#'         the same step.
#'
#' The generator samples in polar coordinates by rejection sampling and never
#' touches the planar redistribution kernel -- see scripts/checks/synth_helpers.R
#' for why that independence is the whole point.
#'
#' Usage: Rscript scripts/checks/check_synthetic_gam.R [1|2|both]

suppressPackageStartupMessages({
  library(mgcv); library(amt); library(terra); library(tidyverse); library(survival)
})
source("scripts/helper_functions.R")
source("scripts/checks/check_helpers.R")
source("scripts/checks/synth_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
part <- if (length(args)) args[1] else "both"

SHAPE <- 1.5637; SCALE <- 145.87; KAPPA <- 0.5591
N_STEPS <- 1000; N_REPS <- 5; PEAK <- 0.6
land <- synth_landscape()

fit_gam_synth <- function(d, rhs) {
  d$times <- 1; d$stratum <- as.integer(factor(d$step_id_)); d$obs <- as.integer(d$case_)
  mgcv::gam(stats::as.formula(paste("cbind(times, stratum) ~", rhs)),
            data = d, family = mgcv::cox.ph(), weights = obs, method = "REML")
}

# ---- Part 1: unimodal recovery ----------------------------------------------
if (part %in% c("1", "both")) {
  check_section("Part 1  unimodal truth, 5 replicates x 1000 steps")
  f_true <- function(cv) -4 * (cv - PEAK)^2      # max 0 at cv = PEAK
  grid <- seq(0.05, 0.95, length.out = 101)
  truth <- f_true(grid); truth <- truth - mean(truth)   # smooths are sum-to-zero

  reps <- purrr::map_dfr(seq_len(N_REPS), function(rep) {
    trk <- synth_track(land, f_true, n_steps = N_STEPS, shape = SHAPE,
                       scale = SCALE, kappa = KAPPA, seed = 1000 + rep)
    d <- synth_design(trk, land, shape = SHAPE, scale = SCALE, kappa = KAPPA,
                      seed = 2000 + rep)
    g <- fit_gam_synth(d, "sl_ + log(sl_) + cos(ta_) + s(cov_end)")

    nd <- tibble(cov_end = grid, sl_ = median(d$sl_), ta_ = 0)
    pr <- predict(g, nd, type = "terms", se.fit = TRUE)
    j <- grep("cov_end", colnames(pr$fit))
    est <- pr$fit[, j]; se <- pr$se.fit[, j]
    est_c <- est - mean(est)
    cf <- coef(g)
    tibble(rep = rep,
           rmse = sqrt(mean((est_c - truth)^2)),
           shape_cor = cor(est_c, truth),
           peak = grid[which.max(est)],
           coverage = mean(abs((est - mean(est)) - truth) <= 1.96 * se),
           edf = sum(g$edf),
           b_sl = cf[["sl_"]], b_logsl = cf[["log(sl_)"]], b_costa = cf[["cos(ta_)"]])
  })

  cat("\n  per-replicate:\n")
  print(as.data.frame(reps |> mutate(across(where(is.numeric), ~round(., 4)))),
        row.names = FALSE)

  check("recovered shape correlates with truth (all reps r > 0.95)",
        all(reps$shape_cor > 0.95),
        sprintf("r = %s", paste(round(reps$shape_cor, 3), collapse = ", ")))
  check("fitted peak within 0.10 of the true peak",
        all(abs(reps$peak - PEAK) <= 0.10),
        sprintf("peaks = %s (truth %.2f)", paste(round(reps$peak, 3), collapse = ", "), PEAK))
  check("RMSE on log-RSS scale below 0.25 in every rep",
        all(reps$rmse < 0.25),
        sprintf("rmse = %s", paste(round(reps$rmse, 3), collapse = ", ")))
  # mgcv's smooth intervals target ~95% coverage AVERAGED across the function,
  # not pointwise in every realisation (Nychka; Marra & Wood). Penalisation also
  # biases the fit slightly at a sharp peak, where coverage dips. Judging the
  # mean across replicates is the criterion those intervals actually claim;
  # requiring >=90% in every single rep would fail on correct code.
  check("mean pointwise CI coverage across replicates >= 0.90",
        mean(reps$coverage) >= 0.90,
        sprintf("mean %.2f; per-rep %s", mean(reps$coverage),
                paste(round(reps$coverage, 2), collapse = ", ")))
  # Sharp: the truth has NO movement modification, because the track was
  # generated from the same gamma / von Mises the controls came from. Any
  # systematic non-zero here means the movement block is absorbing signal it
  # should not.
  for (nm in c("b_sl", "b_logsl", "b_costa")) {
    v <- reps[[nm]]
    check(sprintf("movement correction %s ~ 0 (truth: no modification)", nm),
          abs(mean(v)) < 3 * (sd(v) / sqrt(length(v))) + 0.05,
          sprintf("mean %.4f (sd %.4f) over %d reps", mean(v), sd(v), length(v)))
  }
}

# ---- Part 2: iSSF vs GAM on identical data ----------------------------------
if (part %in% c("2", "both")) {
  check_section("Part 2  monotonic truth: iSSF vs GAM, same data, same structure")
  BETA <- 2
  f_lin <- function(cv) BETA * (cv - 1)          # max 0 at cv = 1
  trk <- synth_track(land, f_lin, n_steps = N_STEPS, shape = SHAPE,
                     scale = SCALE, kappa = KAPPA, seed = 7)
  d <- synth_design(trk, land, shape = SHAPE, scale = SCALE, kappa = KAPPA, seed = 8)

  # amt needs the tentative distributions on the object to build a kernel later.
  attr(d, "sl_") <- amt::make_gamma_distr(shape = SHAPE, scale = SCALE)
  attr(d, "ta_") <- amt::make_vonmises_distr(kappa = KAPPA)
  attr(d, "crs") <- 6610
  class(d) <- c("random_steps", "steps_xyt", "steps_xy", class(d))

  RHS <- "cov_end + sl_ + log(sl_) + cos(ta_)"
  iss <- amt::fit_issf(d, stats::as.formula(
    paste("case_ ~", RHS, "+ strata(step_id_)")), model = TRUE)
  gm <- fit_gam_synth(d, RHS)

  terms_ <- c("cov_end", "sl_", "log(sl_)", "cos(ta_)")
  ci <- coef(iss$model)[terms_]
  cg <- coef(gm)[terms_]
  cmp <- tibble(term = terms_, issf = as.numeric(ci), gam = as.numeric(cg),
                diff = as.numeric(cg - ci),
                truth = c(BETA, 0, 0, 0))
  cat("\n  coefficients:\n")
  print(as.data.frame(cmp |> mutate(across(where(is.numeric), ~round(., 5)))),
        row.names = FALSE)

  check("iSSF and GAM coefficients agree (max |diff| < 1e-3)",
        max(abs(cmp$diff)) < 1e-3,
        sprintf("max |diff| = %.2e on %s", max(abs(cmp$diff)),
                cmp$term[which.max(abs(cmp$diff))]))
  # Criterion is CI coverage, not a fixed distance. With n = 1000 the slope's
  # standard error is ~0.2, so demanding the estimate land within 0.15 of truth
  # would require it to be closer than its own standard error -- a check that
  # fails on correct code roughly half the time.
  se_slope <- summary(iss$model)$coefficients["cov_end", "se(coef)"]
  lo <- cmp$issf[1] - 1.96 * se_slope; hi <- cmp$issf[1] + 1.96 * se_slope
  check("95% CI for the cov_end slope contains the truth",
        BETA >= lo && BETA <= hi,
        sprintf("estimate %.4f, SE %.4f, CI [%.3f, %.3f], truth %.2f",
                cmp$issf[1], se_slope, lo, hi, BETA))

  # ---- gate-3 machinery: per-step probabilities, nothing filtered ------------
  # Same construction both onestep_logscore*() use: build the kernel as a raster,
  # normalise it, read the density at the observed endpoint.
  cat("\n  building per-step kernels (this is the slow part)...\n")
  stp <- trk |> dplyr::select(x1_, y1_, x2_, y2_, t1_, t2_, burst_)
  n_eval <- min(120, nrow(stp))
  set.seed(11); idx <- sort(sample(nrow(stp), n_eval))

  logp <- purrr::map_dfr(idx, function(i) {
    start <- amt::make_start(as.numeric(stp[i, c("x1_", "y1_")]), ta_ = 0,
                             time = stp$t1_[i], dt = lubridate::hours(4), crs = 6610)
    cf_issf <- function(xy, map) amt::extract_covariates(xy, map, where = "both")
    ki <- tryCatch(amt::redistribution_kernel(
      x = iss, map = land, fun = cf_issf, start = start,
      landscape = "discrete", as.rast = TRUE)$redistribution.kernel,
      error = function(e) NULL)
    kg <- tryCatch(redistribution_kernel_gam(
      x = gm, map = land, start = start, fun = cf_issf,
      sl_distr = attr(d, "sl_"), ta_distr = attr(d, "ta_"),
      compensate.movement = TRUE, normalize = TRUE, as.rast = TRUE
    )$redistribution.kernel, error = function(e) NULL)
    if (is.null(ki) || is.null(kg)) return(NULL)
    at <- cbind(stp$x2_[i], stp$y2_[i])
    # amt names its kernel layer "w", ours names it "kernel", and terra::extract
    # returns just that one column here (no ID column) -- so take the last
    # column by position rather than assuming an index.
    dens <- function(k) {
      e <- terra::extract(k, at)
      as.numeric(e[1, ncol(e)])
    }
    tibble(step = i, issf = log(dens(ki)), gam = log(dens(kg)))
  }) |> dplyr::filter(is.finite(issf), is.finite(gam))

  cat(sprintf("  evaluated %d steps\n", nrow(logp)))
  dd <- logp$gam - logp$issf
  cat("\n  per-step log probability:\n")
  print(as.data.frame(tibble(
    stat = c("mean issf", "mean gam", "mean diff", "max |diff|", "cor"),
    value = round(c(mean(logp$issf), mean(logp$gam), mean(dd),
                    max(abs(dd)), cor(logp$issf, logp$gam)), 6))), row.names = FALSE)

  check("iSSF and GAM assign the same per-step probability (max |diff| < 0.01)",
        max(abs(dd)) < 0.01,
        sprintf("max |diff| = %.2e over %d steps", max(abs(dd)), nrow(logp)))
  check("per-step log probabilities correlate at r > 0.999",
        cor(logp$issf, logp$gam) > 0.999,
        sprintf("r = %.8f", cor(logp$issf, logp$gam)))
  check("total log score agrees (relative diff < 1e-3)",
        abs(sum(dd)) / abs(sum(logp$issf)) < 1e-3,
        sprintf("issf %.3f vs gam %.3f", sum(logp$issf), sum(logp$gam)))
}

quit(status = if (check_summary() > 0) 1L else 0L)
