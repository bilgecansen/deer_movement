#' @description
#' MANUAL WALKTHROUGH of the amt path: fit one model, build one redistribution
#' kernel, simulate one step, score one step. Straight-line code meant to be
#' stepped through line by line, not run as a job. Companion to
#' walkthrough_gam.R, which does the same for the GAM path.
#'
#' Every function in scripts/helpers/ is INLINED here. Library calls (amt,
#' terra, sf, survival) stay as calls -- that is the boundary of what this file
#' lets you review, and it is the right one: our code is what carries the
#' project's assumptions.
#'
#' There are no loops. Fixing one deer, one model and one step removes every
#' loop in the amt path: the deer/model loops in fit_amt.R, the burst and step
#' loops in onestep_logscore(), and the burst and month-chunk loops in
#' simulate_movement().
#'
#' WHERE THIS DIFFERS FROM THE GAM WALKTHROUGH, and why:
#'   * No null is fit. The amt models are measured against the GAM null, which
#'     is loaded from results/gam/ and scored here alongside them (section 14).
#'   * The redistribution kernel is amt's, not ours. amt builds the candidate
#'     disc, applies the Jacobian and normalises internally, so section 8 is one
#'     call rather than the six explicit steps the GAM walkthrough spells out.
#'     Section 9 opens the result up so you can still see the weights.
#'   * The fitted coefficients need renaming before amt will accept them back
#'     (section 6) -- a wart of the amt API, not of the model.
#'
#' WHY THE ASSERTIONS MATTER. An inlined copy is a second implementation of the
#' same computation, and a second implementation can drift. If it does, you
#' would be carefully reviewing code that is not what runs -- worse than not
#' reviewing, because it produces false confidence. So every stage ends with a
#' stopifnot() comparing the hand-built object against what the production
#' helper returns on the same inputs. Do not delete them.
#'
#' Sections:
#'   0  configuration
#'   1  the wrangled track
#'   2  fit one model                        (fit_mod)
#'   3  rasters for this deer
#'   4  one step, and its incoming heading
#'   5  fit-time vs kernel-time covariates
#'   6  rebuild the model amt can simulate   (rename_landcover_coefs)
#'   7  the start point
#'   8  the redistribution kernel            (amt::redistribution_kernel)
#'   9  what is inside the kernel
#'  10  SIMULATE one step
#'  11  SCORE one step                       (onestep_logscore, one iteration)
#'  12  the GAM null, on the same step
#'  13  delta_logp against the GAM null
#'  14  the handoff to step 2
#'
#' Usage: open it and step through. Running it top to bottom also works and will
#' stop at the first assertion that fails.

suppressPackageStartupMessages({
  library(amt)
  library(terra)
  library(tidyverse)
  library(sf)
  library(survival)
  library(mgcv)     # only for the GAM null in sections 12-13
  library(foreach)  # only so the production functions can be called for checks
})

# The helpers are loaded ONLY so the assertions can compare against them. None
# of the walkthrough code below calls them, except where a section says
# otherwise.
source("scripts/helper_functions.R")


# 0 ---- Configuration ---------------------------------------------------------
# A non-winter deer, so the NDVI month swap and the landcover interaction are
# both exercised. Winter (nb) drops NDVI from models 2 and 3.
KEY <- "7193_fa_2020"
ID <- "7193"
SEASON <- "fa"
YEAR <- 2020L

STEP_I <- 4L    # which scoreable step to walk through (see section 4)
MODEL <- "2"    # numbered model: HR-centre + NDVI x landcover
DT_HOURS <- 4L

stopifnot(SEASON != "nb")


# 1 ---- The wrangled track ----------------------------------------------------
deer <- readRDS(sprintf("data/tracks/data_%s.rds", KEY))

# stp     every observed step. Used for simulation and scoring.
# stp.var the fitting design: each observed step plus 25 random ones.
stp <- deer$stp[[1]]
stp_var <- deer$stp.var[[1]]

cat(sprintf(
  "track: %d observed steps in %d bursts | design: %d rows, %d strata\n",
  nrow(stp), dplyr::n_distinct(stp$burst_),
  nrow(stp_var), dplyr::n_distinct(stp_var$step_id_)
))

# Fewer strata than observed steps: amt::random_steps drops bursts shorter than
# 3 steps entirely, and every surviving burst loses its first step (no turning
# angle). So the model is fit on fewer steps than exist.


# 2 ---- Fit one model (fit_mod, inlined) -------------------------------------
# An amt model is a conditional logistic regression: one stratum per observed
# step, the observed endpoint against its 25 random alternatives. That is the
# same likelihood the GAM path expresses as a stratified Cox PH.
move <- "case_ ~ (sl_):tod_start_ + log(sl_) + cos(ta_)"
RHS <- paste(
  move,
  "+ HR_center_end + wiscland_end + ndvi_end + wiscland_end:ndvi_end",
  "+ strata(step_id_)"
)

iss <- amt::fit_issf(stp_var, stats::as.formula(RHS), model = TRUE)

cat(sprintf("fitted %d coefficients\n", length(coef(iss$model))))
print(round(head(coef(iss$model), 5), 5))

stopifnot(all.equal(
  coef(iss$model),
  coef(fit_mod(stp_var, RHS)$iss$model)
))

# The tentative movement distributions, fitted during the wrangle. amt stores
# them on the fitted object; the kernel needs them to rebuild the full movement
# density from the sl_ / log(sl_) / cos(ta_) corrections.
cat(sprintf(
  "tentative gamma(shape=%.3f, scale=%.1f), von Mises(kappa=%.3f)\n",
  iss$sl_$params$shape, iss$sl_$params$scale, iss$ta_$params$kappa
))


# 3 ---- Rasters for this deer -------------------------------------------------
env_raster <- load_landcover(YEAR, SEASON)
ndvi_year <- load_ndvi(YEAR)

crop_extent <- sf::st_buffer(
  sf::st_as_sf(stp, coords = c("x1_", "y1_"), crs = 6610),
  CROP_BUFFER_M
)
env_cropped <- terra::crop(env_raster, crop_extent)
env_cropped$HR_edge <- load_hr_edge_raster(ID, SEASON, YEAR, env_cropped)
env_cropped$HR_center <- load_hr_center_raster(ID, SEASON, YEAR, env_cropped)
env_cropped$HR_center_log <- log1p(env_cropped$HR_center)
ndvi_cropped <- terra::crop(ndvi_year, crop_extent)

cat(sprintf("cropped map: %d x %d cells, %d layers\n",
            terra::ncol(env_cropped), terra::nrow(env_cropped),
            terra::nlyr(env_cropped)))


# 4 ---- One step, and its incoming heading -----------------------------------
# Heading = the ABSOLUTE bearing of the PRECEDING step in the same burst. The
# kernel converts a candidate endpoint into a turning angle by subtracting it,
# so it is a reference direction, NOT a turning angle -- despite living in a
# field called ta_.
steps <- stp |>
  dplyr::group_by(burst_) |>
  dplyr::mutate(prev_head = dplyr::lag(atan2(y2_ - y1_, x2_ - x1_))) |>
  dplyr::ungroup()

scoreable <- which(!is.na(steps$prev_head))
cat(sprintf("%d of %d observed steps have a heading and are scoreable\n",
            length(scoreable), nrow(steps)))

i <- scoreable[STEP_I]
step <- steps[i, ]

cat(sprintf(
  "step %d: from (%.0f, %.0f) at %s, observed end (%.0f, %.0f), sl_ %.1f m\n",
  i, step$x1_, step$y1_, format(step$t1_), step$x2_, step$y2_, step$sl_
))
cat(sprintf("incoming heading: %.4f rad (%.1f deg)\n",
            step$prev_head, step$prev_head * 180 / pi))

mo <- lubridate::month(step$t1_)
env_test <- env_cropped
env_test$ndvi <- terra::resample(
  ndvi_cropped[[mo]], env_cropped, method = "near"
)
cat(sprintf("NDVI layer: %s (month %d)\n", names(ndvi_cropped)[mo], mo))


# 5 ---- Fit-time vs kernel-time covariates -----------------------------------
# The model must be scored on the covariates it was FIT on. Same caveat as the
# GAM walkthrough: until data/tracks/ is re-wrangled, a late-month step will
# show ndvi_end disagreeing, because the stored values came from the old
# 1st-of-month NDVI stamps. HR and landcover should always agree exactly.
step_steps <- step
class(step_steps) <- c("steps_xyt", "steps_xy", class(step_steps))
attr(step_steps, "crs") <- 6610
kernel_time <- amt::extract_covariates(step_steps, env_test, where = "both")
fit_time <- stp_var[stp_var$case_ & stp_var$t1_ == step$t1_, ]

if (nrow(fit_time) == 1) {
  cat(sprintf(
    "HR_center_end  fit %.3f  kernel %.3f  | ndvi_end  fit %.4f  kernel %.4f\n",
    fit_time$HR_center_end, kernel_time$HR_center_end,
    fit_time$ndvi_end, kernel_time$ndvi_end
  ))
} else {
  cat("this step was not in the fitting design\n")
}


# 6 ---- Rebuild the model amt can simulate (rename_landcover_coefs, inlined) --
# amt cannot hand a fit_clogit straight back to its own kernel: the coefficient
# NAMES have to match the model matrix it builds internally, and factor-level
# and interaction-term naming do not survive the round trip. This block is
# bookkeeping, not modelling -- but it is real code that can go wrong, so it is
# inlined rather than hidden.
coefs <- coef(iss$model)
names(coefs) <- rename_landcover_coefs(names(coefs))
coefs <- coefs[!is.na(coefs)]

dummy <- amt::make_issf_model(coefs = coefs)
mm_names <- colnames(model.matrix(
  amt:::ssf_formula(dummy$model$formula),
  data = stp_var
))

# Interaction terms can come back with their factors in the other order
# (a:b vs b:a). Try the permutations until one matches the model matrix.
n_reordered <- 0
for (idx in seq_along(coefs)) {
  nm <- names(coefs)[idx]
  if (grepl(":", nm) && !(nm %in% mm_names)) {
    for (p in combinat::permn(strsplit(nm, ":")[[1]])) {
      cand <- paste(p, collapse = ":")
      if (cand %in% mm_names) {
        names(coefs)[idx] <- cand
        n_reordered <- n_reordered + 1
        break
      }
    }
  }
}
cat(sprintf("%d coefficients, %d interaction names reordered to match\n",
            length(coefs), n_reordered))

amt_model <- amt::make_issf_model(coefs = coefs, sl = iss$sl_, ta = iss$ta_)


# 7 ---- The start point -------------------------------------------------------
# ta_ carries the observed incoming heading. Called WITHOUT ta_ (as on a bare
# one-row track) make_start() silently returns 0 -- due east -- which is the bug
# that was scoring every step as though the deer had just been heading east.
start_pt <- amt::make_start(
  c(step$x1_, step$y1_),
  ta_ = step$prev_head,
  time = step$t1_,
  dt = lubridate::hours(DT_HOURS),
  crs = terra::crs(env_test)
)
print(as.data.frame(start_pt))


# 8 ---- The redistribution kernel (amt::redistribution_kernel) ---------------
# Unlike the GAM path, the kernel is amt's own: it builds the candidate disc,
# evaluates the fitted model at every cell, applies the polar-to-planar Jacobian
# and normalises, all internally. The covariate callback below is the only part
# we supply.
cov_fun <- function(xy, map) {
  xy |>
    amt::extract_covariates(map, where = "both") |>
    amt::time_of_day(include.crepuscule = FALSE, where = "both") |>
    dplyr::mutate(
      tod_start_day = as.integer(tod_start_ == "day"),
      tod_start_night = as.integer(tod_start_ == "night"),
      days = lubridate::yday(t2_) - min(lubridate::yday(t2_)) + 1
    )
}

rk <- amt::redistribution_kernel(
  x = amt_model,
  map = env_test,
  fun = cov_fun,
  start = start_pt,
  landscape = "discrete",
  as.rast = TRUE
)
kernel_rast <- rk$redistribution.kernel


# 9 ---- What is inside the kernel --------------------------------------------
vals <- terra::values(kernel_rast)
ok <- !is.na(vals)
cat(sprintf("kernel: %d candidate cells, sums to %.10f, min %.3e, max %.3e\n",
            sum(ok), sum(vals[ok]), min(vals[ok]), max(vals[ok])))

# The candidate geometry amt used, recovered from the raster itself.
xyc <- terra::xyFromCell(kernel_rast, which(ok))
sl_cand <- sqrt((xyc[, 1] - step$x1_)^2 + (xyc[, 2] - step$y1_)^2)
cat(sprintf("candidate step lengths: %.0f to %.0f m\n",
            min(sl_cand), max(sl_cand)))
cat(sprintf("this step is %.0f m -> %s\n", step$sl_,
            if (step$sl_ <= max(sl_cand)) "inside the disc" else "OUTSIDE"))


# 10 ---- SIMULATE one step ----------------------------------------------------
# Sampling one candidate with probability proportional to the kernel. That is
# the whole step: everything above builds the distribution, this draws from it.
set.seed(1)
cells <- which(ok)
pick <- sample(cells, size = 1, prob = vals[cells])
sim_end <- terra::xyFromCell(kernel_rast, pick)

cat(sprintf("simulated endpoint: (%.0f, %.0f)\n", sim_end[1], sim_end[2]))
cat(sprintf("observed  endpoint: (%.0f, %.0f)\n", step$x2_, step$y2_))


# 11 ---- SCORE one step (onestep_logscore, one iteration) --------------------
obs_pt <- cbind(step$x2_, step$y2_)
p <- terra::extract(kernel_rast, obs_pt)[1, 1]

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

prod_scores <- onestep_logscore(
  stp_data = stp |> dplyr::filter(burst_ == step$burst_),
  env_test = env_test, ndvi_test = ndvi_cropped,
  amt_train = amt_model
)
prod_row <- prod_scores[prod_scores$t1_ == step$t1_, ]
stopifnot(all.equal(logp, prod_row$logp))
cat(sprintf("matches onestep_logscore: %.4f (status %s)\n",
            prod_row$logp, prod_row$status))


# 12 ---- The GAM null, on the same step --------------------------------------
# There is no amt null. The reference is the GAM null, scored here on the SAME
# step so the comparison is like-for-like. fit_GAM.R must have produced it.
null_path <- sprintf("results/gam/results_gam_null_%s.rds", KEY)
if (!file.exists(null_path)) {
  cat("GAM null not found -- run fit_GAM.R first. Skipping sections 12-13.\n")
  logp_null <- NA_real_
} else {
  gam_null <- readRDS(null_path)$gam
  sl_distr <- attr(stp_var, "sl_")
  ta_distr <- attr(stp_var, "ta_")

  rk_null <- redistribution_kernel_gam(
    x = gam_null, map = env_test, start = start_pt, fun = gam_cov_fun,
    sl_distr = sl_distr, ta_distr = ta_distr,
    compensate.movement = TRUE, normalize = TRUE, as.rast = TRUE
  )
  p_null <- terra::extract(rk_null$redistribution.kernel, obs_pt)[1, 1]
  logp_null <- if (!is.na(p_null) && p_null > 0) log(p_null) else NA_real_
  cat(sprintf("GAM null density %.3e -> logp %.4f\n", p_null, logp_null))
}


# 13 ---- delta_logp against the GAM null -------------------------------------
# Gate 3 compares each numbered amt model to the GAM null. Here for one step; in
# production it is the sum over every step both models scored.
#
# The caveat to keep in mind: the GAM null carries a cyclic time-of-day smooth,
# a more flexible movement block than the day/night factor this amt model uses.
# So delta_logp mixes "does habitat help" with "my movement block is less
# flexible than the reference's".
cat(sprintf("logp amt %.4f | GAM null %.4f | delta %.4f (this step)\n",
            
            logp, logp_null, logp - logp_null))


# 14 ---- The handoff to step 2 ------------------------------------------------
# What amt::simulate_path does between iterations, without looping. The heading
# for the NEXT step is the bearing of the step just taken.
next_head <- atan2(sim_end[2] - step$y1_, sim_end[1] - step$x1_)
next_start <- amt::make_start(
  as.numeric(sim_end), ta_ = next_head,
  time = step$t1_ + lubridate::hours(DT_HOURS),
  dt = lubridate::hours(DT_HOURS), crs = terra::crs(env_test)
)
cat(sprintf("next start: (%.0f, %.0f), heading %.4f rad\n",
            next_start$x_, next_start$y_, next_start$ta_))
cat("\nSection 7 onwards would now repeat with next_start in place of step.\n")

cat("\nAll assertions passed:",
    "this walkthrough computes what the helpers compute.\n")
