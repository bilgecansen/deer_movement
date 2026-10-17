#' @description
#' Fit the hierarchical GAM of one season: every training deer-year pooled
#' (prep_pool_gam.R) into one step-selection model across deer, fitted as a
#' stratified Cox model in mgcv (cox.ph, Klappstein et al. 2024) by REML.
#'
#' Every term has a population effect and a deviation for each deer (the
#' animal, across its years), shrunk by a penalty estimated from the data:
#'   * straight lines: log(sl_), cos(ta_), the habitat columns and the
#'     landcover indicators (forest is the reference). The population slope
#'     is a parametric term; the deer's deviation is
#'     s(animal, by = x, bs = "re"), a normal random slope.
#'   * the step-length rate by time of day: sl_ times a cyclic curve in
#'     tod_, s(tod_, bs = "cc", by = sl_), with a cyclic deviation curve per
#'     deer sharing one penalty (bs = "fs", Pedersen et al. 2019, model GS).
#'   * distance to the home-range centre: a curve, s(HR_center_end), with a
#'     deviation curve per deer (bs = "fs").
#' The random points were drawn from each deer-year's own gamma and von
#' Mises (stp.var), so the movement terms correct those distributions. A
#' deer's own intercept in the HR deviation curve is constant within a step
#' and cancels; only its penalty holds it at 0.
#'
#' Inputs:  data/gam/pooled_<season>.rds (prep_pool_gam.R)
#' Output:  results/gam/hier_<season>.rds — the fitted gam, the formula, the
#'          run time, and the table's size
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
SEASON <- "fa"
K_TOD <- 10L
K_HR <- 5L

# Load packages ---------------------------------------------------------------
library(mgcv)

# helper functions
source("scripts/helper_functions.R")

pool <- readRDS(sprintf("data/gam/pooled_%s.rds", SEASON))
d <- pool$pooled
d$animal <- factor(d$animal)
d$times <- 1
lc_cols <- paste0(LANDCOVER_LEVELS[-1], "_end")
lin_vars <- c("log_sl", "cos_ta", setdiff(pool$hab_vars, "HR_center_end"),
              lc_cols)

rhs <- c(
  lin_vars,
  sprintf("s(tod_, bs = 'cc', k = %d, by = sl_)", K_TOD),
  sprintf("s(tod_, animal, bs = 'fs', xt = 'cc', k = %d, by = sl_)", K_TOD),
  sprintf("s(HR_center_end, k = %d)", K_HR),
  sprintf("s(HR_center_end, animal, bs = 'fs', k = %d)", K_HR),
  sprintf("s(animal, by = %s, bs = 're')", lin_vars)
)
f <- stats::as.formula(paste("cbind(times, stratum) ~",
                             paste(rhs, collapse = " + ")))
cat(sprintf("%s: %d rows, %d steps, %d deer\n", SEASON, nrow(d),
            max(d$stratum), nlevels(d$animal)))
print(f)

t0 <- Sys.time()
fit <- mgcv::gam(
  f,
  data = d,
  family = mgcv::cox.ph(),
  weights = case_,
  method = "REML",
  knots = list(tod_ = c(0, 24)),
  control = mgcv::gam.control(trace = TRUE)
)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

dir.create("results/gam", showWarnings = FALSE, recursive = TRUE)
out_path <- sprintf("results/gam/hier_%s.rds", SEASON)
saveRDS(list(fit = fit, formula = f, minutes = elapsed, rows = nrow(d),
             steps = max(d$stratum), deer = nlevels(d$animal)),
        out_path)
cat(sprintf("\nfit in %.1f min, %d coefficients -> %s\n", elapsed,
            length(coef(fit)), out_path))
print(fit$outer.info[c("conv", "iter", "grad")])
print(summary(fit))
