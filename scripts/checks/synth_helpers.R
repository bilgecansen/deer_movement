# Synthetic-world generator ---------------------------------------------------
#
# Builds a landscape with a KNOWN selection function and generates tracks from
# it, so the fitting path can be checked against ground truth.
#
# The one design decision that makes this test worth anything: the generator
# samples in POLAR coordinates by rejection sampling, and never touches the
# planar redistribution kernel. Concretely, it draws
#     sl ~ Gamma(shape, scale),  ta ~ vonMises(kappa)
# and accepts the proposed step with probability exp(f_true(covariate at the
# endpoint)). That is the step-selection generative model stated directly, with
# no Jacobian anywhere in it.
#
# If instead we generated with redistribution_kernel_gam(), a wrong polar ->
# planar Jacobian would be applied at generation AND at fitting, cancel out, and
# the test would pass while the bug remained. Rejection sampling in polar space
# is what keeps the generator independent of the thing under test.

#' Continuous covariate field on [0, 1].
#'
#' Wavelengths are deliberately SHORT (900-1500 m) relative to the 99th-percentile
#' step length (~1 km). The covariate has to vary within a single step's reach or
#' the control points in a stratum all carry nearly the same value, the design has
#' no within-stratum contrast, and nothing is identifiable. A first attempt using
#' broad Gaussian bumps (sigma ~2.5 km) failed exactly this way: the walk settled
#' on a plateau and every realised endpoint had cov_ = 1.
#'
#' Not monotonic in space, so the covariate is not a proxy for
#' distance-from-start and cannot be recovered by accident.
synth_field <- function(x, y, cx, cy) {
  a <- sin(2 * pi * (x - cx) / 900) * cos(2 * pi * (y - cy) / 1100)
  b <- sin(2 * pi * ((x - cx) + (y - cy)) / 1500)
  as.numeric((0.6 * a + 0.4 * b + 1) / 2)
}

#' Build the synthetic landscape as a SpatRaster.
#' Carries the layers the production covariate functions expect:
#'   cov      the continuous covariate carrying the selection signal
#'   wiscland a constant categorical layer (gam_cov_fun factors it; unused here)
#' @param res_m cell size in metres; 30 matches the real landcover grid.
synth_landscape <- function(cx = 512437, cy = 289293, half = 12000, res_m = 30) {
  r <- terra::rast(
    xmin = cx - half, xmax = cx + half,
    ymin = cy - half, ymax = cy + half,
    resolution = res_m, crs = "EPSG:6610"
  )
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  # Layer is named "cov" (not "cov_"): amt::extract_covariates appends
  # "_start"/"_end", so "cov" yields cov_start / cov_end, matching the rest of
  # the pipeline's naming. "cov_" would yield "cov__end".
  cov_ <- r; terra::values(cov_) <- synth_field(xy[, 1], xy[, 2], cx, cy)
  names(cov_) <- "cov"
  wl <- r; terra::values(wl) <- 1L; names(wl) <- "wiscland"
  levels(wl) <- data.frame(value = 1L, wiscland = "forest")
  c(cov_, wl)
}

#' Twelve-layer dummy NDVI stack, so the production scoring functions -- which
#' swap NDVI by month unconditionally -- can run against the synthetic world.
synth_ndvi <- function(template, year) {
  n <- terra::rast(template[[1]])
  terra::values(n) <- 0
  s <- terra::rast(rep(list(n), 12))
  names(s) <- month.abb
  terra::time(s) <- ndvi_layer_times(year)
  s
}

#' Generate one track by rejection sampling in polar coordinates.
#'
#' @param f_true  function(cov) -> log selection, with max 0 over [0, 1] so the
#'   acceptance probability exp(f_true) is a valid probability.
#' @return tibble of steps in amt "steps" layout (x1_, y1_, x2_, y2_, t1_, t2_,
#'   sl_, ta_, burst_), plus cov_end, the covariate at the realised endpoint.
synth_track <- function(land, f_true, n_steps = 1000,
                        shape = 1.5637, scale = 145.87, kappa = 0.5591,
                        t0 = as.POSIXct("2020-05-01 00:00:00", tz = "UTC"),
                        dt_h = 4, seed = 1, batch = 4000) {
  set.seed(seed)
  # unname(): terra::ext() yields NAMED values, and c(x1_ = x, ...) would then
  # produce a column called "x1_.xmin" instead of "x1_".
  e <- as.vector(terra::ext(land))
  cx <- unname((e[["xmin"]] + e[["xmax"]]) / 2)
  cy <- unname((e[["ymin"]] + e[["ymax"]]) / 2)
  e <- list(xmin = unname(e[["xmin"]]), xmax = unname(e[["xmax"]]),
            ymin = unname(e[["ymin"]]), ymax = unname(e[["ymax"]]))
  margin <- 2500  # keep the walk clear of the edge so the disc is never clipped

  x <- cx; y <- cy; head_ <- runif(1, -pi, pi)
  out <- vector("list", n_steps)

  for (i in seq_len(n_steps)) {
    accepted <- FALSE
    while (!accepted) {
      sl <- rgamma(batch, shape = shape, scale = scale)
      ta <- circular_rvm(batch, kappa)
      ang <- head_ + ta
      x2 <- x + sl * cos(ang); y2 <- y + sl * sin(ang)
      inside <- x2 > e$xmin + margin & x2 < e$xmax - margin &
        y2 > e$ymin + margin & y2 < e$ymax - margin
      cv <- rep(NA_real_, batch)
      cv[inside] <- terra::extract(land[["cov"]],
                                   cbind(x2[inside], y2[inside]))[, 1]
      p <- ifelse(is.na(cv), 0, exp(f_true(cv)))
      hit <- which(runif(batch) < p)
      if (length(hit)) {
        j <- hit[1]
        out[[i]] <- c(x1_ = x, y1_ = y, x2_ = x2[j], y2_ = y2[j],
                      sl_ = sl[j], ta_ = ta[j], cov_end = cv[j])
        x <- x2[j]; y <- y2[j]; head_ <- ang[j]
        accepted <- TRUE
      }
    }
  }
  d <- tibble::as_tibble(do.call(rbind, out))
  d$t1_ <- t0 + (seq_len(n_steps) - 1) * dt_h * 3600
  d$t2_ <- d$t1_ + dt_h * 3600
  d$burst_ <- 1L
  d
}

#' von Mises draws with mu = 0 (Best & Fisher rejection algorithm), so the
#' generator does not depend on circular:: being installed.
circular_rvm <- function(n, kappa) {
  if (kappa < 1e-8) return(runif(n, -pi, pi))
  a <- 1 + sqrt(1 + 4 * kappa^2)
  b <- (a - sqrt(2 * a)) / (2 * kappa)
  r <- (1 + b^2) / (2 * b)
  out <- numeric(n); k <- 0L
  while (k < n) {
    m <- n - k
    u1 <- runif(m); u2 <- runif(m); u3 <- runif(m)
    z <- cos(pi * u1)
    f <- (1 + r * z) / (r + z)
    c_ <- kappa * (r - f)
    ok <- (c_ * (2 - c_) - u2 > 0) | (log(c_ / u2) + 1 - c_ >= 0)
    if (any(ok)) {
      th <- sign(u3[ok] - 0.5) * acos(f[ok])
      take <- min(length(th), n - k)
      out[(k + 1):(k + take)] <- th[seq_len(take)]
      k <- k + take
    }
  }
  out
}

#' Build the case/control design: each observed step plus `n_ctrl` control steps
#' drawn from the SAME tentative gamma / von Mises proposal the generator used.
#' This is what amt::random_steps does, written out so the synthetic test does
#' not depend on the wrangle (option A).
synth_design <- function(trk, land, n_ctrl = 25,
                         shape = 1.5637, scale = 145.87, kappa = 0.5591,
                         seed = 99) {
  set.seed(seed)
  n <- nrow(trk)
  head_ <- atan2(trk$y2_ - trk$y1_, trk$x2_ - trk$x1_) - trk$ta_

  ctrl <- purrr::map_dfr(seq_len(n), function(i) {
    sl <- rgamma(n_ctrl, shape = shape, scale = scale)
    ta <- circular_rvm(n_ctrl, kappa)
    ang <- head_[i] + ta
    tibble::tibble(
      step_id_ = i, case_ = FALSE,
      x1_ = trk$x1_[i], y1_ = trk$y1_[i],
      x2_ = trk$x1_[i] + sl * cos(ang), y2_ = trk$y1_[i] + sl * sin(ang),
      sl_ = sl, ta_ = ta, t1_ = trk$t1_[i], t2_ = trk$t2_[i], burst_ = 1L
    )
  })
  used <- tibble::tibble(
    step_id_ = seq_len(n), case_ = TRUE,
    x1_ = trk$x1_, y1_ = trk$y1_, x2_ = trk$x2_, y2_ = trk$y2_,
    sl_ = trk$sl_, ta_ = trk$ta_, t1_ = trk$t1_, t2_ = trk$t2_, burst_ = 1L
  )
  d <- dplyr::bind_rows(used, ctrl) |> dplyr::arrange(step_id_, dplyr::desc(case_))
  d$cov_end <- terra::extract(land[["cov"]], cbind(d$x2_, d$y2_))[, 1]
  d$tod_ <- lubridate::hour(d$t1_) + lubridate::minute(d$t1_) / 60
  d |> dplyr::filter(!is.na(cov_end))
}
