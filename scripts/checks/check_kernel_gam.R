#' @description
#' Analytic reductions and contract checks for the GAM path: cases where the
#' right answer is known from the maths or from the design of the data,
#' independent of any fitted model, so the code can be falsified rather than
#' merely read.
#'
#' What each check protects against:
#'   K1 movement kernel       a sign / shape / scale slip in gam_movement_kernel
#'   K2 planar Jacobian       the -log(sl_) term in gam_kernel_weights being
#'                            wrong, missing or double-applied. With every
#'                            coefficient zeroed the kernel MUST reduce to the
#'                            tentative gamma x von Mises density in PLANAR
#'                            coordinates = the polar density divided by sl_.
#'                            Ships with a deliberate-break control so we know
#'                            the check can fail.
#'   K3 shift invariance      exp(w - mean(w)) must leave relative weights alone
#'   K4 rasterisation terra::rast(data.frame(x, y, w)) landing weights on
#'                            the wrong cells; normalisation not summing to 1
#'   K5 fit/predict agreement the covariates the kernel builds at simulation and
#'                            scoring time differing from the ones the model was
#'                            FIT on. Silent, and it invalidates every
#'                            downstream number without raising an error.
#'   K6 design contract       prepare_gam_data's stratified Cox layout
#'   K7 metric identities     calc_energy_score / svf_score against themselves
#'
#' Usage: Rscript scripts/checks/check_kernel_gam.R [<id> <season> <year>]

suppressPackageStartupMessages({
  library(mgcv)
  library(amt)
  library(terra)
  library(tidyverse)
  library(sf)
  library(ctmm)
  library(foreach)
})
source("scripts/helper_functions.R")
source("scripts/checks/check_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 3) {
  id <- args[1]
  season <- args[2]
  year <- as.integer(args[3])
} else {
  id <- "7193"
  season <- "fa"
  year <- 2020L
}
key <- sprintf("%s_%s_%d", id, season, year)
cat(sprintf("GAM kernel checks on deer %s\n", key))

deer <- readRDS(sprintf("data/tracks/data_%s.rds", key))
stp <- deer$stp[[1]]
stp_var <- deer$stp.var[[1]]
sl_distr <- attr(stp_var, "sl_")
ta_distr <- attr(stp_var, "ta_")
shape <- sl_distr$params$shape
scale <- sl_distr$params$scale
kappa <- ta_distr$params$kappa

# ---- K1 ---------------------------------------------------------------------
check_section("K1  gam_movement_kernel vs. independent log densities")
set.seed(1)
xy <- tibble(sl_ = runif(500, 1, 1500), ta_ = runif(500, -pi, pi))
phi <- gam_movement_kernel(xy, sl_distr, ta_distr)
lg <- dgamma(xy$sl_, shape = shape, scale = scale, log = TRUE)
lv <- kappa * cos(xy$ta_) - log(2 * pi * besselI(kappa, 0))
check_close("phi - (log gamma + log vonMises) is constant",
            sd(phi - (lg + lv)), 0, tol = 1e-9)

# ---- K2 ---------------------------------------------------------------------
check_section("K2  zero-coefficient kernel == planar gamma x von Mises")
planar <- function(x) {
  dgamma(x$sl_, shape = shape, scale = scale) / x$sl_ *
    exp(kappa * cos(x$ta_))
}
raw <- 0 + gam_movement_kernel(xy, sl_distr, ta_distr) - log(xy$sl_)
w <- exp(raw - mean(raw))
r <- w / planar(xy)
check_close("kernel / analytic planar density is constant (CV)",
            sd(r) / mean(r), 0, tol = 1e-9)

# Deliberate break: without the Jacobian the ratio must NOT be constant. A check
# that has never been seen to fail is not evidence of anything.
raw_nj <- gam_movement_kernel(xy, sl_distr, ta_distr)
r_nj <- exp(raw_nj - mean(raw_nj)) / planar(xy)
cv_nj <- sd(r_nj) / mean(r_nj)
check("control: omitting the Jacobian breaks K2", cv_nj > 1e-6,
      sprintf("CV without Jacobian = %.4g (must be non-zero)", cv_nj))

# ---- K3 ---------------------------------------------------------------------
check_section("K3  normalised weights invariant to a constant shift in eta")
set.seed(3)
e <- rnorm(400)
a <- exp(e - mean(e))
a <- a / sum(a)
b <- exp((e + 17.3) - mean(e + 17.3))
b <- b / sum(b)
check_close(
  "weights unchanged by +17.3 on eta", max(abs(a - b)), 0, tol = 1e-12
)

# ---- Build the runtime map exactly as the runners do ------------------------
env <- load_landcover(year, season)
ndvi <- load_ndvi(year)
crop_extent <- sf::st_buffer(
  sf::st_as_sf(stp, coords = c("x1_", "y1_"), crs = 6610), CROP_BUFFER_M)
envc <- terra::crop(env, crop_extent)
envc$HR_edge <- load_hr_edge_raster(id, season, year, envc)
envc$HR_center <- load_hr_center_raster(id, season, year, envc)
envc$HR_center_log <- log1p(envc$HR_center)
ndvic <- terra::crop(ndvi, crop_extent)

res <- readRDS(sprintf("results/gam/results_gam_%s.rds", key))
ok_i <- which(!vapply(res, function(x) is.character(x$gam), logical(1)))[1]
gam_fit <- res[[ok_i]]$gam

# ---- K4 ---------------------------------------------------------------------
check_section("K4  rasterised kernel: normalisation and cell alignment")
start <- amt::make_start(
  as.numeric(stp[1, c("x1_", "y1_")]), ta_ = 0,
  time = stp$t1_[1], dt = lubridate::hours(4), crs = 6610)
mo <- lubridate::month(stp$t1_[1])
envk <- envc
envk$ndvi <- terra::resample(ndvic[[mo]], envc)

rk <- redistribution_kernel_gam(
  x = gam_fit, map = envk, start = start, fun = gam_cov_fun,
  sl_distr = sl_distr, ta_distr = ta_distr,
  compensate.movement = TRUE, normalize = TRUE, as.rast = TRUE)

if (is.null(rk)) {
  check("kernel built", NA, "redistribution_kernel_gam returned NULL")
} else {
  kr <- rk$redistribution.kernel
  tot <- terra::global(kr, "sum", na.rm = TRUE)[1, 1]
  check_close("normalised kernel sums to 1", tot, 1, tol = 1e-8)
  vals <- terra::values(kr)
  check("kernel has no negative mass", all(vals[!is.na(vals)] >= 0),
        sprintf("min = %.3g", min(vals, na.rm = TRUE)))

  # Independent recomputation at the raster's own cell centres. If rasterisation
  # misplaces weights this correlation drops below 1.
  cells <- which(!is.na(vals))
  set.seed(4)
  cells <- sample(cells, min(2000, length(cells)))
  xyc <- terra::xyFromCell(kr, cells)
  cand <- tibble(
    x1_ = start$x_[1], y1_ = start$y_[1], x2_ = xyc[, 1], y2_ = xyc[, 2],
    t1_ = start$t_, t2_ = start$t_ + start$dt)
  cand$sl_ <- sqrt((cand$x2_ - cand$x1_)^2 + (cand$y2_ - cand$y1_)^2)
  cand$ta_ <- atan2(cand$y2_ - cand$y1_, cand$x2_ - cand$x1_)
  class(cand) <- c("steps_xyt", "steps_xy", class(cand))
  attr(cand, "crs") <- 6610
  wc <- gam_kernel_weights(gam_cov_fun(cand, envk), gam_fit,
                           sl_distr, ta_distr, TRUE)
  keep <- is.finite(wc) & wc > 0
  cr <- suppressWarnings(cor(log(wc[keep]), log(vals[cells][keep])))
  check("raster cell values match recomputed weights", isTRUE(cr > 0.999),
        sprintf("cor(log w, log raster) = %.6f over %d cells", cr, sum(keep)))
}

# ---- K5 ---------------------------------------------------------------------
# The covariates the kernel sees at simulation/scoring time must be the ones the
# model was fit on. Compare gam_cov_fun at the OBSERVED step endpoints against
# the values extract_step_variables stored for those same steps.
check_section("K5  fit-time vs kernel-time covariates at identical locations")
# Match on t1_, never on position. stp holds every observed step; stp.var holds
# only those that survived random-step generation (first-in-burst steps have no
# turning angle, and 2-step bursts drop out entirely), so the two differ in
# length and a positional comparison silently compares different steps.
used <- stp_var |> dplyr::filter(case_) |> dplyr::arrange(t1_)
obs <- stp |> dplyr::filter(t1_ %in% used$t1_) |> dplyr::arrange(t1_)
check("used rows matched to observed steps by t1_",
      nrow(used) == nrow(obs) &&
        isTRUE(all.equal(as.numeric(used$t1_), as.numeric(obs$t1_))),
      sprintf("%d fit steps of %d observed", nrow(used), nrow(stp)))

# Mirror the runners: swap in the NDVI layer for month(t1_), per month.
obs$.row <- seq_len(nrow(obs))
kern_cov <- purrr::map_dfr(split(obs, lubridate::month(obs$t1_)), function(g) {
  m <- lubridate::month(g$t1_[1])
  e <- envc
  e$ndvi <- terra::resample(ndvic[[m]], envc)
  gg <- g
  class(gg) <- c("steps_xyt", "steps_xy", class(gg))
  attr(gg, "crs") <- 6610
  out <- gam_cov_fun(gg, e)
  out$.row <- g$.row
  out
}) |> dplyr::arrange(.row)

for (v in c("HR_center_end", "HR_edge_end", "ndvi_end")) {
  if (!v %in% names(kern_cov) || !v %in% names(used)) {
    check(sprintf("%s present both sides", v), NA, "missing")
    next
  }
  A <- used[[v]]
  B <- kern_cov[[v]]
  ok <- is.finite(A) & is.finite(B)
  d <- abs(A[ok] - B[ok])
  rel <- mean(d > 1e-6 * pmax(abs(A[ok]), 1))
  passed <- rel == 0
  check(sprintf("%s agrees fit vs kernel", v), passed,
        sprintf("%.1f%% of %d steps differ; max |diff| = %.4g",
                100 * rel, sum(ok), if (length(d)) max(d) else 0))
  if (!passed && v == "ndvi_end") {
    cat("        ^ expected until data/tracks/ is re-wrangled: these\n",
        "         files were built with the old 1st-of-month NDVI\n",
        "         stamps. K5b below tests the current stamps directly\n",
        "         and is the authoritative check.\n", sep = "")
  }
}
agree_lc <- mean(as.character(used$wiscland_end) ==
                   as.character(kern_cov$wiscland_end), na.rm = TRUE)
check("wiscland_end agrees fit vs kernel", isTRUE(agree_lc == 1),
      sprintf("%.2f%% agree", 100 * agree_lc))
check("tod_ agrees fit vs kernel",
      isTRUE(max(abs(used$tod_ - kern_cov$tod_)) < 1e-9),
      sprintf("max |diff| = %.3g", max(abs(used$tod_ - kern_cov$tod_))))

# K5b: the same question asked WITHOUT any stored data, so it validates the
# NDVI timestamps themselves rather than whatever is currently in data/tracks/.
# Layer m is filled with the constant m, so the extracted value names the layer
# amt chose; we then ask whether it equals the step's calendar month, which is
# what the kernel indexes by. Stamped on the 1st this mismatched 171/366 days.
check_section("K5b fit-time layer selection == kernel-time calendar month")
probe <- terra::rast(nrows = 4, ncols = 4, xmin = 0, xmax = 1e3,
                     ymin = 0, ymax = 1e3, nlyrs = 12, crs = "EPSG:6610")
terra::values(probe) <- rep(1:12, each = 16)
names(probe) <- month.abb
terra::time(probe) <- ndvi_layer_times(year)
days <- seq(as.POSIXct(sprintf("%d-01-01 12:00:00", year), tz = "UTC"),
            as.POSIXct(sprintf("%d-12-31 12:00:00", year), tz = "UTC"),
            by = "day")
ps <- tibble(x1_ = 500, y1_ = 500, x2_ = 500, y2_ = 500,
             t1_ = days, t2_ = days + 4 * 3600)
class(ps) <- c("steps_xyt", "steps_xy", class(ps))
attr(ps, "crs") <- 6610
sel <- amt::extract_covariates_var_time(
  ps, probe, max_time = lubridate::days(31),
  when = "any", where = "both", name_covar = "ndvi")$ndvi_end
mism <- which(sel != lubridate::month(days))
lbl <- format(days[mism], "%b %d")

# Known irreducible residual: because months differ in length, the midpoint rule
# cannot put the nearest-layer boundary exactly on the month boundary next to
# February. Measured over 2016-2022 this is always and only Jan 31 / Mar 01
# (1-2 days a year, vs 171 under the old 1st-of-month stamps). Both fall in
# winter, where the models use landcover rather than NDVI, so the practical
# effect on NDVI-using models is nil. Anything OUTSIDE those two days is a
# genuine regression and must fail.
allowed <- c("Jan 31", "Mar 01")
unexpected <- setdiff(lbl, allowed)
check("NDVI layer choice agrees outside the known Feb-adjacent days",
      length(unexpected) == 0,
      sprintf("%d/%d days differ (%s); unexpected: %s",
              length(mism), length(days),
              if (length(lbl)) paste(lbl, collapse = ", ") else "none",
              if (length(unexpected)) {
                paste(unexpected, collapse = ", ")
              } else {
                "none"
              }))

# ---- K6 ---------------------------------------------------------------------
check_section("K6  prepare_gam_data stratified Cox contract")
gd <- prepare_gam_data(stp_var)

# HAZARD (amt, not our code): amt registers group_by.steps_xy and
# summarise.steps_xy, and the pair SILENTLY drops the grouping -- on a steps
# object, group_by(x) |> summarise(...) returns ONE aggregate row instead of one
# row per group, with no warning. group_by |> slice_tail is unaffected, which is
# why overlap_ud() is fine. Always as_tibble() a steps object before grouping.
n_grp <- dplyr::n_distinct(gd$stratum)
raw_grp <- nrow(
  gd |>
    dplyr::group_by(stratum) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop")
)
if (raw_grp != n_grp) {
  cat(sprintf(
    paste("  NOTE: amt's summarise.steps_xy dropped grouping",
          "(%d rows, not %d) -- expected; as_tibble() first.\n"),
    raw_grp, n_grp))
}
per <- tibble::as_tibble(gd) |> dplyr::group_by(stratum) |>
  dplyr::summarise(
    n_used = sum(obs == 1), n_avail = sum(obs == 0), .groups = "drop"
  )
check("exactly one used point per stratum", all(per$n_used == 1),
      sprintf("range %d..%d", min(per$n_used), max(per$n_used)))
check("constant number of available points per stratum",
      dplyr::n_distinct(per$n_avail) == 1,
      sprintf("values: %s", paste(sort(unique(per$n_avail)), collapse = ", ")))
check("event times constant", dplyr::n_distinct(gd$times) == 1,
      sprintf("%d distinct", dplyr::n_distinct(gd$times)))
check("obs is 0/1 only", all(gd$obs %in% c(0, 1)),
      paste(sort(unique(gd$obs)), collapse = ", "))
check("stratum count equals observed step count",
      nrow(per) == sum(stp_var$case_),
      sprintf("%d strata vs %d used rows", nrow(per), sum(stp_var$case_)))
check("as_tibble() restores correct grouping on steps objects",
      nrow(per) == n_grp,
      sprintf("%d groups vs %d distinct strata", nrow(per), n_grp))

# ---- K8 ---------------------------------------------------------------------
# The scoring code must condition on the deer's ACTUAL incoming heading. The
# kernel converts a candidate endpoint to a turning angle via ta_ = bearing -
# start$ta_, so start$ta_ is a reference direction, not a turning angle.
# amt::make_start() on a bare one-row track returns ta_ = 0 (due east); scoring
# every step from such a start evaluates each one as though the deer had just
# been travelling east, which scrambles cos(ta_) by up to ~0.8 log units per
# step and shifts delta_logp by roughly the size of the gate-3 threshold.
check_section("K8  scoring conditions on the observed incoming heading")

bare <- stp[1, c("x1_", "y1_", "t1_")] |>
  amt::make_track(.x = x1_, .y = y1_, .t = t1_, crs = 6610) |>
  amt::make_start()
check("amt::make_start() on one point still defaults to ta_ = 0 (the trap)",
      isTRUE(bare$ta_[1] == 0),
      sprintf("ta_ = %s", bare$ta_[1]))

# First-in-burst steps have no preceding step, so no heading exists; they must
# be skipped rather than scored from an invented direction.
heads <- stp |> dplyr::group_by(burst_) |>
  dplyr::mutate(ph = dplyr::lag(atan2(y2_ - y1_, x2_ - x1_))) |>
  dplyr::ungroup()
n_first <- sum(is.na(heads$ph))

envk2 <- envc
envk2$ndvi <- terra::resample(ndvic[[lubridate::month(stp$t1_[1])]], envc)
ls_out <- onestep_logscore_gam(stp, envk2, ndvic, gam_fit,
                               sl_distr = sl_distr, ta_distr = ta_distr)
check("every first-in-burst step is skipped, every other step scored",
      sum(is.na(ls_out$logp)) == n_first &&
        sum(!is.na(ls_out$logp)) == nrow(stp) - n_first,
      sprintf("%d scored + %d skipped = %d steps (%d are first-in-burst)",
              sum(!is.na(ls_out$logp)), sum(is.na(ls_out$logp)),
              nrow(ls_out), n_first))

# Direct: scoring with the true heading must differ from scoring with ta_ = 0.
# If these ever agree, the heading is being ignored again.
i2 <- which(!is.na(heads$ph))[1]
lp_at <- function(ta_start) {
  st <- amt::make_start(c(heads$x1_[i2], heads$y1_[i2]), ta_ = ta_start,
                        time = heads$t1_[i2],
                        dt = lubridate::hours(4), crs = 6610)
  k <- tryCatch(redistribution_kernel_gam(
    x = gam_fit, map = envk2, start = st, fun = gam_cov_fun,
    sl_distr = sl_distr, ta_distr = ta_distr,
    compensate.movement = TRUE, as.rast = TRUE)$redistribution.kernel,
    error = function(e) NULL)
  if (is.null(k)) return(NA_real_)
  v <- terra::extract(k, cbind(heads$x2_[i2], heads$y2_[i2]))
  log(as.numeric(v[1, ncol(v)]))
}
d_head <- lp_at(heads$ph[i2]) - lp_at(0)
check("the start's heading actually changes the kernel (not ignored)",
      is.finite(d_head) && abs(d_head) > 1e-9,
      sprintf("log p differs by %.4g between true heading and ta_ = 0", d_head))

# ---- K7 ---------------------------------------------------------------------
check_section("K7  metric identities")
sim_self <- purrr::map_dfr(1:5, function(k)
  tibble(x_ = obs$x1_, y_ = obs$y1_, t_ = obs$t1_, nsim = k))
check_close("energy score of observed vs itself is 0",
            calc_energy_score(obs, sim_self), 0, tol = 1e-8)

quit(status = if (check_summary() > 0) 1L else 0L)
