#' @description
#' Tier 1 contract checks for the wrangle stage: everything
#' make_random_pt_extraction() and extract_step_variables() write into
#' data/tracks/data_<key>.rds, re-derived independently and compared against
#' what was stored.
#'
#' The wrangle has no "expected output" to diff against, but almost every
#' column it produces is checkable, because it is either
#'   (a) a deterministic function of other stored columns (sl_ from the
#'       coordinates, tod_ from t1_, HR_center_log from HR_center, the
#'       per-class indicators from the categorical band), or
#'   (b) a raster value at a location that is itself stored (x2_, y2_), so it
#'       can be looked up again straight from the source file.
#'
#' On (b) and independence: the re-extraction here reads the ORIGINAL annual
#' rasters in library/ and the per-deer files in data/HR, at the stored
#' coordinates, with no cropping, no resampling and no season/year plumbing.
#' It shares terra's cell lookup with amt::extract_covariates, so it does not
#' test that; it tests everything wrapped around it -- the crop window, the
#' resample onto the env grid, the seasonal band choice, the year, the NDVI
#' layer choice, and the column naming and row order of a nine-call pipe. That
#' is where this pipeline's bugs have actually been: the NDVI mid-month bug
#' was exactly a layer-choice bug, and this check would have named it.
#'
#' Checks are grouped W1..W6 (see the section banners). Results are counted
#' across deer and reported once per contract, with the offending deer named,
#' rather than one PASS line per deer per check.
#'
#' Usage: Rscript scripts/checks/check_wrangle.R [n_deer] [design]
#'   n_deer — how many deer to check, or "all" (default 25). Deer are taken in
#'            a fixed shuffled order so a partial run is not always the same
#'            alphabetical prefix.
#'   design — "amt" (stp.var, 25 random points, gamma/von Mises proposal),
#'            "nonp" (stp.var.nonp, 50 points, uniform-disc proposal), or
#'            "both". Default "amt".

suppressPackageStartupMessages({
  library(terra)
  library(tidyverse)
  library(amt)
  library(sf)
})
source("scripts/helper_functions.R")
source("scripts/checks/check_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
n_arg <- if (length(args) >= 1) args[1] else "25"
design_arg <- if (length(args) >= 2) args[2] else "amt"
designs <- if (design_arg == "both") c("amt", "nonp") else design_arg

# Tunables ---------------------------------------------------------------------
# TOL_M: metres. Coordinate-derived quantities (sl_) should agree to floating
#   point; a micron is already many orders of magnitude tighter than anything
#   that could matter and still far above double round-off on ~1e5 m easting.
# TOL_RAD: radians, for the turning-angle identity.
# TOL_REL: relative tolerance for raster values, which are float32 on disk.
TOL_M <- 1e-6
TOL_RAD <- 1e-8
TOL_REL <- 1e-6

DESIGN_SPEC <- list(
  amt = list(col = "stp.var", n_pts = 25L),
  nonp = list(col = "stp.var.nonp", n_pts = 50L)
)

# Which covariates actually matter, per season ---------------------------------
# An NA in a covariate a model references silently drops that row at fit time
# (mgcv's default na.action), so "which columns must be clean" is exactly "which
# columns the formulas name". That set is season-dependent, and it changes every
# time a model is added -- so read it from the fitting scripts themselves rather
# than restating it here, where it would rot into a check that passes because it
# is testing last month's model set.

#' Evaluate only the top-level function and scalar-constant definitions in a
#' script, skipping its library() calls, file reads and the work it does when
#' run. Lets a check import fit_GAM.R's make_formulas() without fitting a GAM.
load_defs <- function(path) {
  env <- new.env(parent = globalenv())
  for (e in parse(path)) {
    if (!is.call(e) || length(e) != 3L) {
      next
    }
    # e[[1]] is a call, not a name, for things like pkg::fn(a, b) -- which is
    # also length 3. Test for a name first or as.character() returns the
    # operator's three parts and the comparison below is no longer scalar.
    if (!is.name(e[[1]]) || !as.character(e[[1]]) %in% c("<-", "=")) {
      next
    }
    if (!is.name(e[[2]])) {
      next
    }
    rhs <- e[[3]]
    is_fn <- is.call(rhs) && identical(as.character(rhs[[1]]), "function")
    if (is_fn || is.atomic(rhs)) {
      eval(e, envir = env)
    }
  }
  env
}

GAM_DEFS <- load_defs("scripts/gam/fit_GAM.R")
AMT_DEFS <- load_defs("scripts/amt/fit_amt.R")
for (nm in c("make_formulas", "make_null_formula", "K_TOD", "K_NDVI")) {
  if (is.null(GAM_DEFS[[nm]])) {
    stop(
      "check_wrangle could not import '", nm, "' from scripts/gam/fit_GAM.R. ",
      "The check derives its covariate set from the real model formulas; ",
      "fix the import rather than hardcoding a set here."
    )
  }
}

#' Every variable named by any model this season would actually fit -- the GAM
#' null, the four numbered GAMs, and the four amt models.
model_vars <- function(season) {
  f <- c(
    GAM_DEFS$make_null_formula(GAM_DEFS$K_TOD),
    GAM_DEFS$make_formulas(GAM_DEFS$K_TOD, GAM_DEFS$K_NDVI, season),
    AMT_DEFS$make_formulas(season)
  )
  vars <- lapply(f, function(s) {
    rhs <- sub("^\\s*[^~]*~", "", s)
    all.vars(stats::as.formula(paste("~", rhs)))
  })
  sort(unique(unlist(vars)))
}

# Raster cache -----------------------------------------------------------------
# Annual landcover and NDVI are shared across deer; opening them once per
# (year, season) rather than once per deer is the difference between a 3-minute
# and a 20-minute full run. SpatRasters are lazy, so this caches file handles,
# not pixels.
.RC <- new.env(parent = emptyenv())

cached <- function(k, fn) {
  if (is.null(.RC[[k]])) {
    .RC[[k]] <- fn()
  }
  .RC[[k]]
}

lc_band <- function(year, season) {
  band_for <- c(
    br = "breeding", nb = "non_breeding",
    fa = "fawning", pf = "post_fawning"
  )
  cached(sprintf("lc_%d_%s", year, season), function() {
    terra::rast(sprintf("library/landcover/landcover_%d.tif", year))[[
      band_for[[season]]
    ]]
  })
}

ndvi_month <- function(year, mo) {
  cached(sprintf("ndvi_%d_%02d", year, mo), function() {
    terra::rast(sprintf("library/ndvi/ndvi_%d_%02d.tif", year, mo))
  })
}

# Small utilities --------------------------------------------------------------

#' Wrap an angle to (-pi, pi].
wrap_pi <- function(a) {
  ((a + pi) %% (2 * pi)) - pi
}

#' Extract one layer at a matrix of coordinates, returning a bare vector.
#' terra::extract on a points matrix returns a one-column data.frame; going
#' through [[1]] rather than positional indexing keeps this correct if terra
#' ever adds an ID column.
xtract <- function(r, xy) {
  terra::extract(r, xy)[[1]]
}

#' NA-aware equality: TRUE where both are NA, or both are equal and non-NA.
#' Plain == returns NA for either side missing, and sum(NA) then hides the
#' comparison instead of failing it.
identical_na <- function(a, b) {
  (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
}

#' Record one contract result for one deer.
#'
#' n_bad / n_tot are counts, so the aggregate across deer is a plain sum and a
#' contract that could not be evaluated for a deer (n_tot = 0) neither passes
#' nor fails on that deer's account.
#'
#' `worst` is the largest discrepancy this deer saw, reported even when the
#' check passes -- per check_helpers.R, a check that only ever prints PASS is
#' untested itself, and "max |diff| = 2e-12 over 54,288 rows" is the difference
#' between a check that compared something and one that compared nothing.
#' `extra` counts a tolerated or context sub-population and is summed across
#' deer. `note` is free text for the failing case only.
rec <- function(acc, name, n_bad, n_tot, note = "",
                worst = NA_real_, worst_lab = "max |diff|",
                extra = 0, extra_lab = "") {
  acc[[length(acc) + 1L]] <- tibble::tibble(
    name = name,
    n_bad = as.numeric(n_bad),
    n_tot = as.numeric(n_tot),
    note = note,
    worst = as.numeric(worst),
    worst_lab = worst_lab,
    extra = as.numeric(extra),
    extra_lab = extra_lab
  )
  acc
}

#' Largest element, or NA for an empty / all-NA vector. max() on an empty
#' numeric returns -Inf with a warning, which would read as a real measurement.
safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else max(x)
}

# Per-deer evaluation ----------------------------------------------------------

check_one_deer <- function(path, spec) {
  key <- gsub("^data_(.*)\\.rds$", "\\1", basename(path))
  d <- readRDS(path)
  acc <- list()

  if (!spec$col %in% names(d) || is.null(d[[spec$col]][[1]])) {
    return(list(key = key, rows = rec(list(), "design column present", 1, 1)))
  }

  id <- d$id
  season <- d$season
  year <- as.integer(d$year)

  # as_tibble first, always. HAZARD (amt, not our code): amt registers
  # group_by.steps_xy / summarise.steps_xy and the pair SILENTLY drops the
  # grouping, returning one aggregate row instead of one row per stratum. Every
  # per-stratum count below would be wrong, and quietly so.
  sv <- tibble::as_tibble(d[[spec$col]][[1]])
  stp <- tibble::as_tibble(d$stp[[1]])
  xy_end <- cbind(sv$x2_, sv$y2_)
  is_case <- sv$case_

  # ---- W1  geometry and stratum contract ------------------------------------
  # Pure functions of the stored coordinates: no raster, no model, no amt.

  sl_recomp <- sqrt((sv$x2_ - sv$x1_)^2 + (sv$y2_ - sv$y1_)^2)
  acc <- rec(
    acc, "W1.1 sl_ equals distance between stored endpoints",
    sum(abs(sl_recomp - sv$sl_) > TOL_M), nrow(sv),
    worst = safe_max(abs(sl_recomp - sv$sl_)), worst_lab = "max |diff| m"
  )

  # The turning angle is measured from the INCOMING heading, which is not
  # stored. Recover it from the observed step of each stratum
  # (head_in = bearing(case) - ta_case), then every control's ta_ must be its
  # own bearing measured from that same head_in. This is the identity the
  # ta_ = 0 heading bug violated downstream, asserted here at the source.
  bearing <- atan2(sv$y2_ - sv$y1_, sv$x2_ - sv$x1_)
  head_in <- tibble::tibble(
    step_id_ = sv$step_id_, b = bearing, ta = sv$ta_, case_ = is_case
  ) |>
    dplyr::filter(case_) |>
    dplyr::transmute(step_id_, head_in = wrap_pi(b - ta))
  ta_chk <- tibble::tibble(step_id_ = sv$step_id_, b = bearing, ta = sv$ta_) |>
    dplyr::left_join(head_in, by = "step_id_") |>
    dplyr::mutate(gap = abs(wrap_pi(b - head_in - ta)))
  gap_ok <- is.finite(ta_chk$gap)
  acc <- rec(
    acc, "W1.2 ta_ is the bearing measured from the stratum's heading",
    sum(ta_chk$gap[gap_ok] > TOL_RAD), sum(gap_ok),
    worst = safe_max(ta_chk$gap[gap_ok]), worst_lab = "max |diff| rad"
  )

  strata <- sv |>
    dplyr::group_by(step_id_) |>
    dplyr::summarise(
      n_case = sum(case_),
      n_ctrl = sum(!case_),
      n_x1 = dplyr::n_distinct(x1_),
      n_y1 = dplyr::n_distinct(y1_),
      n_t1 = dplyr::n_distinct(t1_),
      n_t2 = dplyr::n_distinct(t2_),
      .groups = "drop"
    )
  acc <- rec(
    acc, "W1.3 exactly one observed step per stratum",
    sum(strata$n_case != 1), nrow(strata),
    sprintf("case counts seen: %s",
            paste(sort(unique(strata$n_case)), collapse = ", "))
  )
  acc <- rec(
    acc, sprintf("W1.4 exactly %d random points per stratum", spec$n_pts),
    sum(strata$n_ctrl != spec$n_pts), nrow(strata),
    sprintf("control counts seen: %s",
            paste(sort(unique(strata$n_ctrl)), collapse = ", "))
  )
  acc <- rec(
    acc, "W1.5 stratum shares one start point and one time interval",
    sum(strata$n_x1 != 1 | strata$n_y1 != 1 |
          strata$n_t1 != 1 | strata$n_t2 != 1),
    nrow(strata)
  )

  # The observed step in stp.var must be the observed step in stp -- same
  # geometry, not merely a plausible one. Matched on t1_, never positionally:
  # stp holds every step, stp.var only those that survived random-step
  # generation, so the two differ in length.
  used <- sv |> dplyr::filter(is_case) |> dplyr::arrange(t1_)
  obs <- stp |> dplyr::filter(t1_ %in% used$t1_) |> dplyr::arrange(t1_)
  if (nrow(used) == nrow(obs) && nrow(used) > 0) {
    dgeom <- pmax(
      abs(used$x1_ - obs$x1_), abs(used$y1_ - obs$y1_),
      abs(used$x2_ - obs$x2_), abs(used$y2_ - obs$y2_)
    )
    acc <- rec(
      acc, "W1.6 observed steps carried through unchanged from stp",
      sum(dgeom > TOL_M), nrow(used),
      worst = safe_max(dgeom), worst_lab = "max |diff| m"
    )
  } else {
    acc <- rec(
      acc, "W1.6 observed steps carried through unchanged from stp",
      1, 1, sprintf("%d used vs %d matched in stp", nrow(used), nrow(obs))
    )
  }

  # ---- W2  derived-column identities ----------------------------------------
  # Every one of these is a column the wrangle computed from another column it
  # also stored, so the two must agree exactly or one of them is stale.

  tod_recomp <- lubridate::hour(sv$t1_) + lubridate::minute(sv$t1_) / 60
  acc <- rec(
    acc, "W2.1 tod_ equals the decimal hour of t1_",
    sum(abs(tod_recomp - sv$tod_) > 1e-9), nrow(sv),
    worst = safe_max(abs(tod_recomp - sv$tod_)), worst_lab = "max |diff| h"
  )

  acc <- rec(
    acc, "W2.2 tod_start_day / _night agree with the tod_start_ factor",
    sum(
      sv$tod_start_day != as.integer(sv$tod_start_ == "day") |
        sv$tod_start_night != as.integer(sv$tod_start_ == "night")
    ),
    nrow(sv)
  )
  # day and night are complementary here because include.crepuscule = FALSE.
  acc <- rec(
    acc, "W2.3 tod_start_day and _night partition every row",
    sum(sv$tod_start_day + sv$tod_start_night != 1), nrow(sv)
  )

  for (sfx in c("start", "end")) {
    a <- sv[[sprintf("HR_center_log_%s", sfx)]]
    b <- log1p(sv[[sprintf("HR_center_%s", sfx)]])
    ok <- is.finite(a) & is.finite(b)
    acc <- rec(
      acc, sprintf("W2.4 HR_center_log_%s == log1p(HR_center_%s)", sfx, sfx),
      sum(abs(a[ok] - b[ok]) > TOL_REL * pmax(abs(b[ok]), 1)), sum(ok),
      worst = safe_max(abs(a[ok] - b[ok]))
    )
  }

  days_recomp <- lubridate::yday(sv$t2_) - min(lubridate::yday(sv$t2_)) + 1
  acc <- rec(
    acc, "W2.5 days is yday(t2_) offset from the deer's first day",
    sum(days_recomp != sv$days), nrow(sv)
  )

  # Per-class indicators vs the categorical band they were built from. forest
  # is the reference level and open_water the exclusion class; neither has an
  # indicator layer, so for both the whole indicator block reads as zeros.
  pred_levels <- setdiff(LANDCOVER_LEVELS, "forest")
  ind_bad <- 0
  ind_tot <- 0
  for (lv in pred_levels) {
    col <- sprintf("%s_end", lv)
    if (!col %in% names(sv)) {
      next
    }
    expect <- as.integer(!is.na(sv$wiscland_end) &
                           as.character(sv$wiscland_end) == lv)
    ind_bad <- ind_bad + sum(sv[[col]] != expect, na.rm = TRUE)
    ind_tot <- ind_tot + nrow(sv)
  }
  acc <- rec(
    acc, "W2.6 per-class indicators agree with wiscland_end",
    ind_bad, ind_tot
  )

  ind_sum <- rowSums(
    as.matrix(sv[, sprintf("%s_end", pred_levels), drop = FALSE])
  )
  acc <- rec(
    acc, "W2.7 indicators mutually exclusive (row sum 0 or 1)",
    sum(!ind_sum %in% c(0, 1)), nrow(sv),
    sprintf("row sums seen: %s",
            paste(sort(unique(ind_sum)), collapse = ", "))
  )

  # ---- W3  independent re-extraction from the source rasters ----------------
  # Everything above stays inside the stored table. These go back to the files.

  lc <- lc_band(year, season)
  hr_files <- c(
    HR_bin = sprintf("data/HR/HRbin_%s.tif", key),
    HR_edge = sprintf("data/HR/HRedge_%s.tif", key),
    HR_center = sprintf("data/HR/HRcenter_%s.tif", key)
  )

  # Grid alignment is a precondition for every exact comparison in W3: the
  # wrangle resamples the HR rasters onto the landcover grid, and resampling is
  # an identity only when the grids already agree. If they do not, the stored
  # HR values are an interpolation of the source and W3.3 SHOULD differ -- so
  # assert alignment first, or a later failure is unreadable.
  src_hr <- terra::rast(hr_files[["HR_center"]])
  aligned <- isTRUE(all.equal(terra::res(src_hr), terra::res(lc))) &&
    isTRUE(all.equal(terra::origin(src_hr), terra::origin(lc)))
  acc <- rec(
    acc, "W3.1 HR rasters share the landcover grid (resample is identity)",
    as.integer(!aligned), 1,
    sprintf("res %s vs %s; origin %s vs %s",
            paste(terra::res(src_hr), collapse = "x"),
            paste(terra::res(lc), collapse = "x"),
            paste(terra::origin(src_hr), collapse = ","),
            paste(terra::origin(lc), collapse = ","))
  )

  direct_lc <- as.character(xtract(lc, xy_end))
  stored_lc <- as.character(sv$wiscland_end)
  # open_water is absent from LANDCOVER_LEVELS, so the wrangle's factor() turns
  # a water endpoint into NA. That is the expected mapping, not a mismatch --
  # compare against it explicitly rather than dropping NAs and calling it a
  # pass.
  expect_lc <- ifelse(direct_lc %in% LANDCOVER_LEVELS, direct_lc, NA_character_)
  acc <- rec(
    acc, "W3.2 wiscland_end == the source band at (x2_, y2_)",
    sum(!identical_na(stored_lc, expect_lc)), nrow(sv)
  )

  n_out_extent <- 0
  for (nm in names(hr_files)) {
    r <- terra::rast(hr_files[[nm]])
    direct <- xtract(r, xy_end)
    stored <- sv[[sprintf("%s_end", nm)]]
    inside <- !is.na(direct)
    n_out_extent <- max(n_out_extent, sum(!inside))
    gap <- abs(stored[inside] - direct[inside])
    acc <- rec(
      acc, sprintf("W3.3 %s_end == the source raster at (x2_, y2_)", nm),
      sum(gap > TOL_REL * pmax(abs(direct[inside]), 1)), sum(inside),
      worst = safe_max(gap)
    )
  }
  # Endpoints beyond the source HR extent are not an error -- the loaders fill
  # them -- but every one of them gets the SAME constant (min for edge, max for
  # centre), a flat plateau the model reads as real signal. CROP_BUFFER_M is
  # documented to sit inside the HR rasters' own buffer, so the count should be
  # zero and a non-zero count means that invariant has slipped.
  acc <- rec(
    acc, "W3.4 no endpoint falls outside the per-deer HR raster extent",
    n_out_extent, nrow(sv)
  )

  # NDVI. load_ndvi() stacks the twelve layers of the deer's `year` field and
  # nothing else, so a step in any other calendar year has no layer within
  # extract_covariates_var_time's 31-day window and comes back NA. Check that
  # first: it is a precondition for the two comparisons after it, and for the
  # non-breeding season it is where the interesting failure lives.
  step_year <- lubridate::year(sv$t1_)
  in_year <- step_year == year
  acc <- rec(
    acc, "W3.5 every step falls inside the year whose NDVI stack was loaded",
    sum(!in_year), nrow(sv),
    sprintf("years present: %s (loaded %d)",
            paste(sort(unique(step_year)), collapse = ", "), year)
  )

  # Reproduce amt's layer choice from the timestamps rather than assuming it,
  # then read that month's file directly. This tests the whole chain -- stamp,
  # selection rule, crop -- in one comparison. Restricted to in-year steps,
  # since out-of-year steps have no correct answer to compare against.
  lt <- ndvi_layer_times(year)
  m_star <- rep(NA_integer_, nrow(sv))
  m_star[in_year] <- vapply(
    sv$t2_[in_year],
    function(tt) which.min(abs(as.numeric(difftime(tt, lt, units = "secs")))),
    integer(1)
  )
  ndvi_direct <- rep(NA_real_, nrow(sv))
  for (mo in sort(unique(m_star[in_year]))) {
    sel <- in_year & m_star == mo
    ndvi_direct[sel] <- xtract(
      ndvi_month(year, mo), xy_end[sel, , drop = FALSE]
    )
  }
  okn <- in_year & is.finite(ndvi_direct) & is.finite(sv$ndvi_end)
  acc <- rec(
    acc, "W3.6 ndvi_end == the nearest-in-time source layer at (x2_, y2_)",
    sum(abs(sv$ndvi_end[okn] - ndvi_direct[okn]) >
          TOL_REL * pmax(abs(ndvi_direct[okn]), 1)),
    sum(okn),
    worst = safe_max(abs(sv$ndvi_end[okn] - ndvi_direct[okn]))
  )
  # The layer amt picks at fit time must be the layer the kernel indexes at
  # simulation / scoring time, which is month(t1_). This is the mid-month bug
  # asked on real steps instead of a synthetic probe (check K5b).
  #
  # A residual survives the mid-month restamp, and it is confined to the last
  # hours of a month. Two things stack: the stamps are 00:00 UTC while the
  # timestamps are America/Chicago, and months differ in length, so the
  # midpoint-to-midpoint boundary lands up to ~12 h off the month boundary
  # (measured: the May/June boundary for 2017 falls at 31 May 13:00 CDT, so
  # that day's last two steps take June's layer). At a 4-hour fix rate that is
  # 2-3 steps per month boundary.
  #
  # So the contract is tightened rather than relaxed: a mismatch is tolerated
  # ONLY on the first or last calendar day of a month, and fails anywhere else.
  # The original 1st-of-month bug mismatched uniformly across days 16-31, so it
  # would still fail this loudly -- what is excluded is only the boundary
  # hours, whose mechanism is known and whose size is bounded.
  mism <- in_year & m_star != lubridate::month(sv$t1_)
  dom <- lubridate::mday(sv$t1_)
  last_dom <- lubridate::mday(
    lubridate::ceiling_date(sv$t1_, "month") - lubridate::days(1)
  )
  at_boundary <- dom == 1 | dom == last_dom
  acc <- rec(
    acc, "W3.7 fit-time NDVI layer == calendar month, away from month ends",
    sum(mism & !at_boundary), sum(in_year),
    extra = sum(mism & at_boundary),
    extra_lab = "month-end mismatches tolerated"
  )

  # ---- W4  water exclusion and NA propagation -------------------------------

  ctrl_water <- sum(direct_lc[!is_case] == "open_water", na.rm = TRUE)
  acc <- rec(
    acc, "W4.1 no random point falls on open_water",
    ctrl_water, sum(!is_case)
  )
  # Observed steps are NOT filtered for water (only the controls are), so a
  # deer standing at a water-classified cell keeps its step -- and picks up
  # wiscland_end = NA, which mgcv then drops. Counted, not asserted away.
  case_water <- sum(direct_lc[is_case] == "open_water", na.rm = TRUE)
  acc <- rec(
    acc, "W4.2 no observed step ends on open_water",
    case_water, sum(is_case)
  )

  # NA in a covariate this season's models reference. Season-aware on purpose:
  # winter formulas name no NDVI term, so an all-NA ndvi_end costs a winter fit
  # nothing -- while the same column being NA for a fawning deer would silently
  # delete the row. Checking one fixed list for every season would either cry
  # wolf on winter or go quiet on the case that matters.
  covs <- intersect(model_vars(season), names(sv))
  for (v in covs) {
    acc <- rec(
      acc, sprintf("W4.3 %s has no NA (in seasons that use it)", v),
      sum(is.na(sv[[v]])), nrow(sv)
    )
  }

  # ---- W5  the proposal's own support ---------------------------------------

  acc <- rec(
    acc, "W5.1 all step lengths finite and strictly positive",
    sum(!is.finite(sv$sl_) | sv$sl_ <= 0), nrow(sv)
  )
  ta_fin <- sv$ta_[is.finite(sv$ta_)]
  acc <- rec(
    acc, "W5.2 all turning angles inside (-pi, pi]",
    sum(ta_fin <= -pi | ta_fin > pi), length(ta_fin),
    worst = safe_max(abs(ta_fin)), worst_lab = "max |ta_|"
  )
  if (identical(spec$col, "stp.var.nonp")) {
    # The uniform-disc proposal has a hard radius: R = the longest observed
    # step. A control beyond it is outside the support the design claims.
    R <- max(stp$sl_, na.rm = TRUE)
    acc <- rec(
      acc, "W5.3 uniform-disc controls inside R = max observed step",
      sum(sv$sl_[!is_case] > R + TOL_M), sum(!is_case),
      worst = safe_max(sv$sl_[!is_case] - R),
      worst_lab = sprintf("longest control minus R (R = %.0f m)", R)
    )
  }

  # ---- W6  what a fit would silently lose -----------------------------------
  # mgcv drops incomplete cases after prepare_gam_data has already numbered the
  # strata. Losing the observed step costs the stratum its event and it drops
  # out of the partial likelihood entirely; losing controls shrinks the risk
  # set, which biases the case's apparent preference upward. Neither warns.
  complete <- stats::complete.cases(sv[, covs, drop = FALSE])
  lost <- sv |>
    dplyr::mutate(.ok = complete) |>
    dplyr::group_by(step_id_) |>
    dplyr::summarise(
      case_lost = any(case_ & !.ok),
      ctrl_lost = sum(!case_ & !.ok),
      .groups = "drop"
    )
  acc <- rec(
    acc, "W6.1 no stratum loses its observed step to an NA covariate",
    sum(lost$case_lost), nrow(lost)
  )
  acc <- rec(
    acc, "W6.2 no stratum loses random points to an NA covariate",
    sum(lost$ctrl_lost > 0), nrow(lost),
    extra = sum(lost$ctrl_lost), extra_lab = "control rows lost"
  )

  list(key = key, rows = acc)
}

# Driver -----------------------------------------------------------------------

files <- sort(list.files("data/tracks", "^data_.*\\.rds$", full.names = TRUE))
if (length(files) == 0) {
  cat("No files in data/tracks/ -- nothing to check.\n")
  quit(status = 0L)
}
set.seed(1L)
files <- files[sample.int(length(files))]
if (n_arg != "all") {
  files <- head(files, as.integer(n_arg))
}

for (dz in designs) {
  spec <- DESIGN_SPEC[[dz]]
  if (is.null(spec)) {
    stop("Unknown design '", dz, "' (expected 'amt', 'nonp' or 'both')")
  }
  check_section(sprintf(
    "Wrangle contracts -- %s design (%s), %d deer",
    dz, spec$col, length(files)
  ))

  out <- purrr::map(files, function(f) {
    res <- tryCatch(
      check_one_deer(f, spec),
      error = function(e) {
        list(
          key = gsub("^data_(.*)\\.rds$", "\\1", basename(f)),
          rows = rec(list(), "W0.0 deer evaluated without error", 1, 1,
                     conditionMessage(e))
        )
      }
    )
    cat(".")
    utils::flush.console()
    res
  })
  cat("\n")

  tally <- purrr::map_dfr(out, function(r) {
    dplyr::bind_rows(r$rows) |> dplyr::mutate(key = r$key)
  })

  # One line per contract, totalled over deer. Sorted by the W-number so the
  # order is the same whichever deer happened to be checked first -- the
  # season-dependent W4.3 lines exist only for the deer they apply to.
  for (nm in sort(unique(tally$name))) {
    sub <- dplyr::filter(tally, name == nm)
    bad_deer <- dplyr::filter(sub, n_bad > 0)
    n_bad <- sum(sub$n_bad)
    n_tot <- sum(sub$n_tot)
    notes <- unique(bad_deer$note[nzchar(bad_deer$note)])

    # The observed measurement, shown pass or fail: the worst discrepancy any
    # deer saw, and any tolerated sub-population, so a green line still states
    # what it measured rather than only that it was happy.
    worst <- safe_max(sub$worst)
    stats <- character(0)
    if (!is.na(worst)) {
      stats <- c(stats, sprintf("%s = %.3g", sub$worst_lab[1], worst))
    }
    if (sum(sub$extra) > 0) {
      stats <- c(stats, sprintf(
        "%s %s", format(sum(sub$extra), big.mark = ","), sub$extra_lab[1]
      ))
    }
    if (n_bad > 0 && length(notes)) {
      stats <- c(stats, notes[1])
    }
    stat_str <- if (length(stats)) {
      sprintf(" [%s]", paste(stats, collapse = "; "))
    } else {
      ""
    }

    detail <- if (n_tot == 0) {
      "nothing to compare"
    } else if (n_bad == 0) {
      sprintf("%s rows OK%s", format(n_tot, big.mark = ","), stat_str)
    } else {
      sprintf(
        "%s/%s rows in %d/%d deer (%s)%s",
        format(n_bad, big.mark = ","), format(n_tot, big.mark = ","),
        nrow(bad_deer), nrow(sub),
        paste(head(bad_deer$key, 3), collapse = ", "),
        stat_str
      )
    }
    check(nm, if (n_tot == 0) NA else n_bad == 0, detail)

    # The stored NDVI values, unlike everything else here, depend on WHEN
    # data/tracks/ was built: files written before the mid-month restamp chose
    # the following month's layer for any step past mid-month. Say so rather
    # than leaving a red line that looks like a live defect.
    if (n_bad > 0 && grepl("^W3\\.6", nm)) {
      cat("        ^ expected until data/tracks/ is re-wrangled: these files\n",
          "         were built with the old 1st-of-month NDVI stamps. W3.7\n",
          "         tests the current stamp rule directly and is unaffected.\n",
          sep = "")
    }
  }
}

quit(status = if (check_summary() > 0) 1L else 0L)
