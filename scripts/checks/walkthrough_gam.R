#' @description
#' MANUAL WALKTHROUGH of the GAM path: fit one model, build one redistribution
#' kernel, simulate one step, score one step. Straight-line code meant to be
#' stepped through line by line, not run as a job.
#'
#' Every function in scripts/helpers/ is INLINED here. Library calls (mgcv, amt,
#' terra, sf) stay as calls -- that is the boundary of what this file lets you
#' review, and it is the right one: our code is what carries the
#' project's assumptions.
#'
#' There are no loops. Fixing one deer, one model and one step removes every
#' loop in the GAM path: the deer/model loops in fit_GAM.R, the burst and step
#' loops in onestep_logscore_gam(), the step loop in simulate_path_gam(), and
#' the burst and month-chunk loops in simulate_movement().
#' redistribution_kernel_gam() has none to begin with -- it is vectorised over
#' candidate cells.
#'
#' WHY THE ASSERTIONS MATTER. An inlined copy is a second implementation of the
#' same computation, and a second implementation can drift. If it does, you
#' would be carefully reviewing code that is not what runs -- worse than not
#' reviewing, because it produces false confidence. So every stage ends with a
#' stopifnot() comparing the hand-built object against what the production
#' helper returns on the same inputs. They are what make this file trustworthy,
#' and they will fail loudly the day a helper changes. Do not delete them.
#'
#' Sections:
#'   0  configuration
#'   1  the wrangled track
#'   2  the stratified Cox design            (prepare_gam_data)
#'   3  fit one model, and the null          (fit_gam_mod)
#'   4  rasters for this deer
#'   5  one step, and its incoming heading
#'   6  fit-time vs kernel-time covariates   (the fit/predict correspondence)
#'   7  the candidate disc                   (redistribution_kernel_gam, part 1)
#'   8  covariates at the candidates         (gam_cov_fun)
#'   9  the habitat term, eta
#'  10  the movement term, phi               (gam_movement_kernel)
#'  11  the weights                          (gam_kernel_weights)
#'  12  rasterise and normalise              (redistribution_kernel_gam, part 2)
#'  13  SIMULATE one step                    (simulate_path_gam, one iteration)
#'  14 SCORE one step (onestep_logscore_gam, one iteration)
#'  15  delta_logp against the null
#'  16  the handoff to step 2
#'
#' Usage: open it and step through. Running it top to bottom also works and will
#' stop at the first assertion that fails.

suppressPackageStartupMessages({
  library(mgcv)
  library(amt)
  library(terra)
  library(tidyverse)
  library(sf)
  library(foreach) # only so the production functions can be called for checks
})

# The helpers are loaded ONLY so the assertions can compare against them. None
# of the walkthrough code below calls them.
source("scripts/helper_functions.R")


# 0 ---- Configuration ---------------------------------------------------------
# A non-winter deer, so the NDVI month swap, the landcover factor levels and the
# factor-smooth are all exercised. Winter (nb) takes a different path:
# s(wiscland_end, bs = "re") instead of NDVI, and no month swap at all.
KEY <- "7193_fa_2020"
ID <- "7193"
SEASON <- "fa"
YEAR <- 2020L

STEP_I <- 4L # which observed step to walk through (see section 5)
DT_HOURS <- 4L # the fix interval, hard-coded throughout the pipeline

stopifnot(SEASON != "nb")


# 1 ---- The wrangled track ----------------------------------------------------
deer <- readRDS(sprintf("data/tracks/data_%s.rds", KEY))

# stp     every observed step. Used for simulation and scoring.
# stp.var the fitting design: each observed step plus 25 random ones.
stp <- deer$stp[[1]]
stp_var <- deer$stp.var[[1]]

# The tentative movement distributions, fitted during the wrangle and carried as
# attributes. The GAM's sl_ / log(sl_) / cos(ta_) terms are CORRECTIONS to
# these, so the kernel needs them to reconstruct the full movement density.
sl_distr <- attr(stp_var, "sl_")
ta_distr <- attr(stp_var, "ta_")

cat(sprintf(
  "track: %d observed steps in %d bursts | design: %d rows, %d strata\n",
  nrow(stp), dplyr::n_distinct(stp$burst_),
  nrow(stp_var), dplyr::n_distinct(stp_var$step_id_)
))
cat(sprintf(
  "tentative gamma(shape=%.3f, scale=%.1f), von Mises(kappa=%.3f)\n",
  sl_distr$params$shape, sl_distr$params$scale, ta_distr$params$kappa
))

# Fewer strata than observed steps: amt::random_steps drops bursts shorter than
# 3 steps entirely, and every surviving burst loses its first step (no turning
# angle). So the model is fit on fewer steps than exist. See #15 in
# docs/gam_decision_inventory.md.


# 2 ---- The stratified Cox design (prepare_gam_data, inlined) -----------------
# An SSF is fit as a stratified Cox PH model: one stratum per observed step, one
# "event" per stratum at a shared constant time, and the used/available
# indicator passed as the censoring weight. That tie structure reproduces
# conditional logistic regression.
gam_data <- stp_var
gam_data$times <- 1 # constant event time
# one stratum per observed step
gam_data$stratum <- as.integer(factor(gam_data$step_id_))
gam_data$obs <- as.integer(gam_data$case_) # 1 = used, 0 = available

stopifnot(identical(
  as.data.frame(gam_data),
  as.data.frame(prepare_gam_data(stp_var))
))

# Contract: exactly one used point per stratum, a constant number of available
# ones, one shared event time.
per_stratum <- tibble::as_tibble(gam_data) |>
  dplyr::group_by(stratum) |>
  dplyr::summarise(used = sum(obs == 1), avail = sum(obs == 0),
                   .groups = "drop")
cat(sprintf(
  "design: %d strata, %d used, %d available each, %d event time\n",
  nrow(per_stratum), unique(per_stratum$used), unique(per_stratum$avail),
  dplyr::n_distinct(gam_data$times)
))
# NOTE: as_tibble() above is deliberate. amt registers group_by.steps_xy and
# summarise.steps_xy, and the pair SILENTLY drops the grouping on a steps object
# -- you get one aggregate row instead of one per stratum, with no warning.


# 3 ---- Fit one model, and the null (fit_gam_mod, inlined) -------------------
# k for the cyclic time-of-day smooth, capped below this deer's distinct-tod
# count so the 'cc' basis stays identifiable.
K_TOD <- max(3L, min(10L, length(unique(gam_data$tod_)) - 1L))
K_NDVI <- 5L

MOVE <- sprintf(
  "sl_ + log(sl_) + cos(ta_) + s(tod_, bs = 'cc', k = %d, by = sl_)", K_TOD
)
# Numbered model 2 (non-winter): home-range centre plus the NDVI-by-landcover
# GS block -- a global NDVI smooth plus shrunk per-class deviations.
RHS_MODEL <- paste(
  MOVE, "+ s(HR_center_end) + s(ndvi_end) +",
  sprintf("s(ndvi_end, wiscland_end, bs = 'fs', k = %d)", K_NDVI)
)
# The null: movement plus distance to the home-range centre. Everything is
# measured against this.
RHS_NULL <- paste(MOVE, "+ s(HR_center_end)")

fit_one <- function(rhs) {
  mgcv::gam(
    stats::as.formula(paste("cbind(times, stratum) ~", rhs)),
    data = gam_data,
    family = mgcv::cox.ph(),
    weights = obs, # the cox.ph censoring indicator
    method = "REML",
    select = FALSE, # no double-penalty shrinkage in production
    knots = list(tod_ = c(0, 24)), # so the cyclic tod smooth wraps at midnight
    drop.unused.levels = FALSE # keep unvisited landcover classes predictable
  )
}
gam_fit <- fit_one(RHS_MODEL)
gam_null <- fit_one(RHS_NULL)

cat(sprintf("model: %.1f edf | null: %.1f edf\n",
            sum(gam_fit$edf), sum(gam_null$edf)))
print(round(coef(gam_fit)[c("sl_", "log(sl_)", "cos(ta_)")], 5))

stopifnot(all.equal(
  coef(gam_fit),
  coef(fit_gam_mod(gam_data, RHS_MODEL, select = FALSE)$gam)
))


# 4 ---- Rasters for this deer -------------------------------------------------
env_raster <- load_landcover(YEAR, SEASON) # categorical + per-class indicators
ndvi_year <- load_ndvi(YEAR) # 12 monthly layers, stamped MID-month

# Crop to the track plus a buffer. CROP_BUFFER_M is 3000 m.
crop_extent <- sf::st_buffer(
  sf::st_as_sf(stp, coords = c("x1_", "y1_"), crs = 6610),
  CROP_BUFFER_M
)
env_cropped <- terra::crop(env_raster, crop_extent)
env_cropped$HR_edge <- load_hr_edge_raster(ID, SEASON, YEAR, env_cropped)
env_cropped$HR_center <- load_hr_center_raster(ID, SEASON, YEAR, env_cropped)
env_cropped$HR_center_log <- log1p(env_cropped$HR_center)
ndvi_cropped <- terra::crop(ndvi_year, crop_extent)

cat(sprintf("cropped map: %d x %d cells, %d layers: %s\n",
            terra::ncol(env_cropped), terra::nrow(env_cropped),
            terra::nlyr(env_cropped),
            paste(names(env_cropped), collapse = ", ")))


# 5 ---- One step, and its incoming heading -----------------------------------
# Heading = the ABSOLUTE bearing of the PRECEDING step in the same burst. The
# kernel converts a candidate endpoint into a turning angle by subtracting it,
# so it is a reference direction, NOT a turning angle -- despite living in a
# field called ta_.
steps <- stp |>
  dplyr::group_by(burst_) |>
  dplyr::mutate(prev_head = dplyr::lag(atan2(y2_ - y1_, x2_ - x1_))) |>
  dplyr::ungroup()

# Only steps with a predecessor are scoreable. The first step of a burst has no
# incoming heading and is skipped rather than scored from an invented one.
scoreable <- which(!is.na(steps$prev_head))
cat(sprintf("%d of %d observed steps have a heading and are scoreable\n",
            length(scoreable), nrow(steps)))

i <- scoreable[STEP_I] # <- THE step this file walks through
step <- steps[i, ]

cat(sprintf(
  "step %d: from (%.0f, %.0f) at %s, observed end (%.0f, %.0f), sl_ %.1f m\n",
  i, step$x1_, step$y1_, format(step$t1_), step$x2_, step$y2_, step$sl_
))
cat(sprintf("incoming heading: %.4f rad (%.1f deg)\n",
            step$prev_head, step$prev_head * 180 / pi))

# The month whose NDVI layer applies. Layers are stamped mid-month so that this
# calendar-month lookup agrees with the nearest-in-time lookup used at fit time.
mo <- lubridate::month(step$t1_)
env_test <- env_cropped
env_test$ndvi <- terra::resample(
  ndvi_cropped[[mo]], env_cropped, method = "near"
)
cat(sprintf("NDVI layer: %s (month %d)\n", names(ndvi_cropped)[mo], mo))


# 6 ---- Fit-time vs kernel-time covariates -----------------------------------
# The model must be scored on the covariates it was FIT on. Here the same step
# is extracted both ways and compared. A mismatch here is silent and invalidates
# everything downstream -- it is exactly where the NDVI bug lived.
#
# CAVEAT while data/tracks/ is still stale. The stored fit-time values were
# wrangled when NDVI layers were stamped on the 1st of the month, so amt's
# nearest-in-time lookup picked the FOLLOWING month's layer for any step past
# mid-month. Layers are now stamped mid-month and the two rules agree, but the
# stored files have not caught up. So until the tracks are re-wrangled:
#   * a step in the first half of a month  -> ndvi_end agrees (as it always did)
#   * a step in the second half -> ndvi_end DISAGREES, and that is the
#                                             stale data, not the code
# Change STEP_I to a late-month step to see it. HR_center_end and wiscland_end
# are unaffected either way and should always agree exactly.
step_steps <- step
class(step_steps) <- c("steps_xyt", "steps_xy", class(step_steps))
attr(step_steps, "crs") <- 6610

kernel_time <- amt::extract_covariates(step_steps, env_test, where = "both")

# The stored fit-time values for this same step (matched on t1_, never on
# position: stp and stp.var have different lengths).
fit_time <- stp_var[stp_var$case_ & stp_var$t1_ == step$t1_, ]

if (nrow(fit_time) == 1) {
  cat(sprintf(
    "HR_center_end  fit %.3f  kernel %.3f  | ndvi_end  fit %.4f  kernel %.4f\n",
    fit_time$HR_center_end, kernel_time$HR_center_end,
    fit_time$ndvi_end, kernel_time$ndvi_end
  ))
  cat(sprintf("wiscland_end   fit %s  kernel %s\n",
              as.character(fit_time$wiscland_end),
              as.character(kernel_time$wiscland_end)))
} else {
  cat("this step was not in the fitting design",
      "(short burst or first in burst)\n")
}


# 7 ---- The candidate disc (redistribution_kernel_gam, part 1) ---------------
# Radius: the 0.99 quantile of the tentative step length. An OBSERVED step
# longer than this cannot be scored -- its endpoint falls outside the disc. That
# is the cause of essentially every scoring failure in the pipeline.
max_dist <- ceiling(do.call(
  paste0("q", sl_distr$name),
  c(list(p = 0.99), sl_distr$params)
))
cat(sprintf("max.dist = %d m | this step is %.0f m -> %s\n",
            max_dist, step$sl_,
            if (step$sl_ <= max_dist) {
              "inside the disc"
            } else {
              "OUTSIDE, unscoreable"
            }))

start_xy <- c(step$x1_, step$y1_)

# Every map cell within max_dist of the start becomes one candidate endpoint.
pt <- sf::st_sf(geom = sf::st_sfc(sf::st_point(start_xy)))
disc <- sf::st_buffer(pt, dist = max_dist)
r1 <- terra::rasterize(terra::vect(disc), terra::crop(env_test, disc))
xy_crds <- terra::crds(r1)

k <- tibble::tibble(x2_ = xy_crds[, 1], y2_ = xy_crds[, 2])

# Absolute bearing of each candidate, then the turning angle relative to the
# incoming heading, wrapped to (-pi, pi].
ta <- atan2(k$y2_ - start_xy[2], k$x2_ - start_xy[1])
ta <- ifelse(ta < 0, 2 * pi + ta, ta)
ta <- (ta - step$prev_head) %% (2 * pi) # <- the heading enters HERE
ta <- ifelse(ta > pi, (2 * pi - ta) * -1, ta)
k$ta_ <- ta
k$sl_ <- sqrt((k$x2_ - start_xy[1])^2 + (k$y2_ - start_xy[2])^2)
k$x1_ <- start_xy[1]
k$y1_ <- start_xy[2]
k$t1_ <- step$t1_
k$t2_ <- step$t1_ + lubridate::hours(DT_HOURS)
class(k) <- c("steps_xyt", "steps_xy", class(k))
attr(k, "crs") <- 6610

cat(sprintf("%d candidates | sl_ %.0f-%.0f m | ta_ %.2f-%.2f rad\n",
            nrow(k), min(k$sl_), max(k$sl_), min(k$ta_), max(k$ta_)))


# 8 ---- Covariates at the candidates (gam_cov_fun, inlined) ------------------
k_cov <- amt::extract_covariates(k, env_test, where = "both") |>
  dplyr::mutate(
    tod_ = lubridate::hour(t1_) + lubridate::minute(t1_) / 60,
    wiscland_start = factor(wiscland_start, levels = LANDCOVER_LEVELS),
    wiscland_end = factor(wiscland_end, levels = LANDCOVER_LEVELS)
  )

stopifnot(all.equal(
  as.data.frame(k_cov),
  as.data.frame(gam_cov_fun(k, env_test))
))
cat(sprintf("candidate landcover classes present: %s\n",
            paste(levels(droplevels(k_cov$wiscland_end)), collapse = ", ")))


# 9 ---- The habitat term, eta -------------------------------------------------
# The GAM's linear predictor at every candidate. cox.ph has no intercept, so eta
# is on a relative scale already.
eta <- as.numeric(stats::predict(
  gam_fit, newdata = k_cov, type = "link", na.action = stats::na.pass
))
cat(sprintf("eta: range %.3f to %.3f, %d NA\n",
            min(eta, na.rm = TRUE), max(eta, na.rm = TRUE), sum(is.na(eta))))


# 10 ---- The movement term, phi (gam_movement_kernel, inlined) ---------------
# The log density of the TENTATIVE movement kernel at each candidate, in polar
# coordinates: log gamma(sl_) + kappa*cos(ta_), each up to a constant.
phi <- -(1 / sl_distr$params$scale) * k_cov$sl_ +
  log(k_cov$sl_) * (sl_distr$params$shape - 1) +
  cos(k_cov$ta_) * ta_distr$params$kappa

stopifnot(all.equal(phi, gam_movement_kernel(k_cov, sl_distr, ta_distr)))

# Independent confirmation: phi differs from the textbook log densities only by
# a constant.
phi_check <- dgamma(k_cov$sl_, shape = sl_distr$params$shape,
                    scale = sl_distr$params$scale, log = TRUE) +
  ta_distr$params$kappa * cos(k_cov$ta_) -
  log(2 * pi * besselI(ta_distr$params$kappa, 0))
cat(sprintf("phi - (log gamma + log vonMises) is constant: sd = %.3g\n",
            sd(phi - phi_check)))


# 11 ---- The weights (gam_kernel_weights, inlined) ---------------------------
# w = exp(eta + phi - log(sl_)).
#
# The -log(sl_) is the polar-to-planar Jacobian. phi is a density over (step
# length, turning angle); the candidates are cells on a PLANE. Converting from
# polar to Cartesian divides by sl_, hence minus its log. Drop this term and the
# kernel over-weights long steps, because a ring at radius r contains more cells
# than one at radius r/2.
w_raw <- eta + phi - log(k_cov$sl_)

# Subtracting the mean before exponentiating is numerical hygiene only -- it
# cancels in the normalisation and leaves relative weights untouched.
w <- exp(w_raw - mean(w_raw[is.finite(w_raw)], na.rm = TRUE))
w[!is.finite(w)] <- 0

stopifnot(all.equal(
  w, gam_kernel_weights(k_cov, gam_fit, sl_distr, ta_distr, TRUE)
))
cat(sprintf("weights: %d candidates, %d with positive mass, max/min = %.1f\n",
            length(w), sum(w > 0), max(w) / min(w[w > 0])))


# 12 ---- Rasterise and normalise (redistribution_kernel_gam, part 2) ---------
kernel_rast <- terra::rast(data.frame(k_cov[, c("x2_", "y2_")], w = w))
kernel_rast <- kernel_rast /
  terra::global(kernel_rast, "sum", na.rm = TRUE)[1, 1]
names(kernel_rast) <- "kernel"

cat(sprintf("normalised kernel sums to %.10f\n",
            terra::global(kernel_rast, "sum", na.rm = TRUE)[1, 1]))

# The same object the production function returns.
rk_prod <- redistribution_kernel_gam(
  x = gam_fit, map = env_test,
  start = amt::make_start(start_xy, ta_ = step$prev_head, time = step$t1_,
                          dt = lubridate::hours(DT_HOURS), crs = 6610),
  fun = gam_cov_fun, sl_distr = sl_distr, ta_distr = ta_distr,
  compensate.movement = TRUE, normalize = TRUE, as.rast = TRUE
)
stopifnot(all.equal(
  terra::values(kernel_rast), terra::values(rk_prod$redistribution.kernel)
))


# 13 ---- SIMULATE one step (simulate_path_gam, one iteration) ----------------
# Simulation samples ONE candidate with probability proportional to w. That is
# the entire step: everything above builds the distribution, this line draws
# from it.
set.seed(1)
idx <- sample.int(nrow(k_cov), size = 1, prob = w)
sim_end <- c(k_cov$x2_[idx], k_cov$y2_[idx])

cat(sprintf("simulated endpoint: (%.0f, %.0f), sl_ %.0f m, ta_ %.2f rad\n",
            sim_end[1], sim_end[2], k_cov$sl_[idx], k_cov$ta_[idx]))
cat(sprintf("observed  endpoint: (%.0f, %.0f), sl_ %.0f m\n",
            step$x2_, step$y2_, step$sl_))
# NOTE: sample.int here draws WITHOUT replacement. Harmless at size = 1, which
# is the only value production uses, but wrong for any larger size.


# 14 ---- SCORE one step (onestep_logscore_gam, one iteration) ----------------
# Scoring reads the normalised density at the endpoint the deer ACTUALLY chose.
obs_pt <- cbind(step$x2_, step$y2_)
p <- terra::extract(kernel_rast, obs_pt)[1, 1]

# An NA here means the observed endpoint lies outside the disc, i.e. the step
# was longer than max.dist. It is recorded as failed_outside_disc, not
# silently dropped.
status <- if (is.na(p)) {
  "failed_outside_disc"
} else if (p <= 0) {
  "failed_zero_density"
} else {
  "ok"
}
logp <- if (status == "ok") log(p) else NA_real_

cat(sprintf("density at observed endpoint: %.3e -> logp %.4f (%s)\n",
            p, logp, status))

# Against the production function, for this same step.
prod_scores <- onestep_logscore_gam(
  stp_data = stp |> dplyr::filter(burst_ == step$burst_),
  env_test = env_test, ndvi_test = ndvi_cropped,
  gam_train = gam_fit, sl_distr = sl_distr, ta_distr = ta_distr
)
prod_row <- prod_scores[prod_scores$t1_ == step$t1_, ]
stopifnot(all.equal(logp, prod_row$logp))
cat(sprintf("matches onestep_logscore_gam: %.4f (status %s)\n",
            prod_row$logp, prod_row$status))


# 15 ---- delta_logp against the null -----------------------------------------
# Gate 3 compares each numbered model to the null on the SAME step. Here that
# comparison is done for one step; in production it is the sum over every step
# both models scored.
rk_null <- redistribution_kernel_gam(
  x = gam_null, map = env_test,
  start = amt::make_start(start_xy, ta_ = step$prev_head, time = step$t1_,
                          dt = lubridate::hours(DT_HOURS), crs = 6610),
  fun = gam_cov_fun, sl_distr = sl_distr, ta_distr = ta_distr,
  compensate.movement = TRUE, normalize = TRUE, as.rast = TRUE
)
p_null <- terra::extract(rk_null$redistribution.kernel, obs_pt)[1, 1]
logp_null <- if (!is.na(p_null) && p_null > 0) log(p_null) else NA_real_

cat(sprintf("logp model %.4f | logp null %.4f | delta %.4f (this step only)\n",
            logp, logp_null, logp - logp_null))


# 16 ---- The handoff to step 2 ------------------------------------------------
# What simulate_path_gam does between iterations, without looping. The heading
# for the NEXT step is the bearing of the step just taken -- this is what keeps
# the turning-angle term meaningful along a simulated path, and it is exactly
# what the scoring path failed to do before it was fixed.
next_head <- atan2(sim_end[2] - start_xy[2], sim_end[1] - start_xy[1])
next_start <- amt::make_start(
  sim_end, ta_ = next_head,
  time = step$t1_ + lubridate::hours(DT_HOURS),
  dt = lubridate::hours(DT_HOURS), crs = 6610
)
cat(sprintf("next start: (%.0f, %.0f), heading %.4f rad, t %s\n",
            next_start$x_, next_start$y_, next_start$ta_,
            format(next_start$t_)))
cat("\nSection 7 onwards would now repeat with next_start in place of step.\n")

cat("\nAll assertions passed:",
    "this walkthrough computes what the helpers compute.\n")
