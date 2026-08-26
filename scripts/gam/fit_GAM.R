#' @description
#' Fit deer movement GAMs (SSFs with penalised smooths) for every deer with a
#' wrangled track file in data/tracks/. This is the GAM analogue of fit_amt.R.
#'
#' Following Klappstein et al. (2024, Methods Ecol Evol), an (integrated) SSF is
#' fit as a stratified Cox proportional-hazards model in mgcv:
#'   gam(cbind(times, stratum) ~ <smooth/linear terms>,
#'       data = ..., family = cox.ph, weights = obs)
#' where `times` is a constant event time, `stratum` groups each observed step
#' with its random points, and `obs` is the used/available (1/0) indicator
#' passed as the cox.ph censoring weight. The reshape is done by
#' prepare_gam_data(); the per-model fit by fit_gam_mod() (both in
#' helper_functions.R).
#'
#' Two things are fit and saved per deer, deliberately kept in separate files:
#'   * the NULL model (movement + s(HR_center_end)) — the reference every
#'     numbered model is scored against downstream — to
#'     results/gam/results_gam_null_<key>.rds;
#'   * the four numbered candidate models, to
#'     results/gam/results_gam_<key>.rds.
#' Keeping them apart means the null can be refit or swapped without touching
#' the candidate set (and vice versa), and downstream code never has to know
#' which slot "the null" happens to occupy.
#'
#' For each (id, season, year):
#'   * If both output rds files already exist and `overwrite` is FALSE, skip.
#'   * Otherwise read the wrangled data, reshape the chosen design's covariate
#'     table (`gam_input`) into mgcv cox.ph form, fit the null formula and all
#'     candidate formulas with fit_gam_mod() (which returns failure-status
#'     objects rather than erroring per-model), and save both files.
#' Errors are caught per deer so a single failure does not halt the loop.
#'
#' Random points: two designs are wrangled per deer.
#'   * "stp.var.nonp" — random steps sampled uniformly over a disc of radius
#'     ~ max step length (wrangle_deer_mvt.R, model = "nonp"). Required for
#'     *non-parametric* movement smooths s(sl_) / s(ta_) (Klappstein et al.
#'     2024, Section 3.1); also valid for parametric kernels.
#'   * "stp.var" — gamma / von Mises (amt-style). Usable for GAMs with a
#'     *parametric* movement kernel, and more efficient there.
#' Pick via `gam_input` below.
#'
#' Configuration: edit `overwrite` / `gam_input` below before running.

# Configuration ---------------------------------------------------------------
overwrite <- T # set to TRUE to refit deer that already have output

# Which per-deer design to fit:
#   "stp.var"      — gamma / von Mises (parametric movement kernel; default)
#   "stp.var.nonp" — uniform-disc (non-parametric movement kernel)
gam_input <- "stp.var"

# Basis / penalty settings ----------------------------------------------------
# K_TOD   max basis dimension for the cyclic time-of-day movement smooths,
#         capped per deer at (distinct tod values - 1) so the 'cc' smooth stays
#         identifiable. Raised 8 -> 10: at 8 the tod smooth was k-bound (edf ~6
#         vs k'=7) on 84 deer; a K_TOD sweep showed edf stabilises at ~6 by k=10
#         (flags clear), while k>=12 over-parameterises relative to the real ~12
#         distinct tod positions and triggers convergence step-failures. The
#         per-deer cap uses the RAW distinct-tod count (~24-50; inflated by
#         non-4h data gaps shifting the fix phase + ~1-min timestamp jitter), so
#         it never binds -- the true resolution ceiling is ~12 (see MOVE comment
#         and doc gam_modeling_decisions.md sec 4).
# K_NDVI basis dimension (per landcover class) of the NDVI 'fs' smooths. NDVI is
#         continuous with rich support, so it is NOT tod-constrained; the
#         k-audit reports whether it is ever the binding constraint (it has wide
#         headroom at 5, so bumping it is unnecessary unless the audit
#         says otherwise).
# SELECT  double-penalty shrinkage (Marra & Wood 2011). FALSE for the production
#         fit: with a fixed, deliberate model set, between-model filtering is
#         the selection step, so we don't also want per-deer within-model term
#         removal (it makes "model 4" a different effective structure on
#         different deer). REML smoothing penalties still regularise each fit;
#         dropping the extra null-space penalty is also faster and converges
#         more reliably. (TRUE was useful during exploration.)
K_TOD <- 10L
K_NDVI <- 5L
SELECT <- FALSE

# Load packages ---------------------------------------------------------------
library(mgcv)
library(tidyverse)
library(furrr)

# helper functions
source("scripts/helper_functions.R")

# Formulas --------------------------------------------------------------------
# make_null_formula() builds the NULL model and make_formulas() the four
# numbered candidates, both for a given cyclic-tod basis (k_tod, capped per
# deer) and NDVI-by-landcover 'fs' basis (k_ndvi). fit_gam_mod prepends the
# "cbind(times, stratum) ~" Cox-PH response.
#
# NULL model = movement + s(HR_center_end). Distance to the home-range centre is
# the minimal "this deer has a home range" description: any candidate model must
# beat it to count as evidence that habitat (not just central-place tethering)
# drives selection. It is fit and stored on its own (see header) and is NOT one
# of the numbered models.
#
# The HR_edge model (formerly numbered 2) is retired: HR_center is the null, and
# the two are near-redundant descriptions of the same tethering. The HR_edge_end
# covariate itself is untouched upstream — it is still computed in the wrangle
# and available in the covariate table, just not fit here.
#
# MOVE = parametric gamma / von Mises movement kernel. The three movement
# covariates enter as parametric main effects (sl_ + log(sl_) + cos(ta_)), which
# carry the baseline correction to the tentative kernel (sl_/log(sl_) -> gamma
# rate/shape, cos(ta_) -> von Mises concentration). Step length sl_ additionally
# enters as a single cyclic-spline interaction with time of day (zebra model;
# Klappstein et al. 2024): the by= smooth adds the (shrinkable) time-of-day
# modulation of movement rate. k_tod caps wiggliness below each deer's
# distinct-tod count (the RAW count is ~24-50, inflated by non-4h data gaps that
# shift the fix phase off the grid plus ~1-min timestamp jitter; the meaningful
# resolution is ~12). Landcover interactions use a global smooth plus
# hierarchical "fs" deviations (Pedersen et al. 2019 "Model GS": one shrunk
# curve per class around a shared global response) rather than independent by=
# smooths, which are unidentifiable when a deer visits only a few classes. The
# global smooth mirrors the amt main effect and lets sparse classes pool toward
# the shared response instead of toward zero.
#
# Models 2 & 3 (breeding) put NDVI in exactly this GS form over ALL landcover
# classes: a global s(ndvi_end) plus an fs interaction s(ndvi_end,
# wiscland_end). The fs carries each class's own intercept AND NDVI curve, so
# per-class landcover selection is part of the GS term -- no separate random
# intercept is needed (Pedersen et al. 2019; adding one would double-count the
# intercept). Where NDVI is uninformative (e.g. forest / wetlands / developed)
# the shared penalty just shrinks that class's deviation toward the global. GS
# also lets us predict classes a deer never visited (Model GS supports
# unobserved levels; GI does not). Model 4 uses the same GS form for HR_center
# (s(HR_center_end) + hc_fs).

# Parametric movement kernel, shared by the null and every numbered model.
make_move <- function(k_tod) {
  sprintf(
    paste0(
      "sl_ + log(sl_) + cos(ta_) + ",
      "s(tod_, bs = 'cc', k = %1$d, by = sl_)"
    ),
    k_tod
  )
}

# NULL model: movement + distance to home-range centre.
make_null_formula <- function(k_tod) {
  paste(make_move(k_tod), "+ s(HR_center_end)")
}

make_formulas <- function(k_tod, k_ndvi, season) {
  move <- make_move(k_tod)
  fs <- sprintf("s(ndvi_end, wiscland_end, bs = 'fs', k = %d)", k_ndvi)
  hc_fs <- sprintf("s(HR_center_end, wiscland_end, bs = 'fs', k = %d)", k_ndvi)

  # Always-fit structural models
  f1 <- move # 1 movement only
  f4 <- paste(move, "+ s(HR_center_end) +", hc_fs) # 4 HR-center x landcover

  # Resource models for slots 2 & 3 — season-dependent. Non-breeding (winter)
  # deer have no valid MODIS NDVI (snow / dormancy -> all-or-mostly NA), so the
  # NDVI models (2 & 3) can't be fit; we substitute the landcover-only models
  # and record them in slots 2 & 3. Either way the object holds 4 models
  # and slots 2/3 always mean "the resource-selection model".
  if (season == "nb") {
    # Winter NDVI is unusable; the landcover-only substitutes use a penalised
    # random effect on landcover, s(wiscland_end, bs = 're'), rather than a
    # fixed factor. It shrinks sparse classes and avoids the near-separation
    # that makes the raw factor slow / non-convergent in cox.ph (validated:
    # ~17x faster, no step failures).
    f2 <- paste(move, "+ s(HR_center_end) + s(wiscland_end, bs = 're')")
    f3 <- paste(move, "+ s(wiscland_end, bs = 're')")
  } else {
    # NDVI x landcover, GS form: global s(ndvi_end) + per-class fs deviations
    # (shared penalty). The fs carries each class's intercept, so landcover
    # selection rides along with the NDVI smooth -- no separate re needed.
    f2 <- paste(move, "+ s(HR_center_end) + s(ndvi_end) +", fs)
    f3 <- paste(move, "+ s(ndvi_end) +", fs)
  }

  stats::setNames(c(f1, f2, f3, f4), c("1", "2", "3", "4"))
}

# Discover deer from data/tracks/ ---------------------------------------------
track_files <- list.files(
  "data/tracks",
  pattern = "^data_.*\\.rds$",
  full.names = TRUE
)

keys <- gsub("^data_(.*)\\.rds$", "\\1", basename(track_files))

cat(sprintf("Found %d wrangled deer in data/tracks/\n", length(track_files)))

# Output dir
dir.create("results/gam", showWarnings = FALSE, recursive = TRUE)

# Loop (parallel across deer) -------------------------------------------------
# fit_GAM is low-memory tabular work (no rasters), so deer are independent,
# CPU-bound units; parallelising across deer gives the best load-balancing (a
# free worker grabs the next deer rather than waiting on a deer's slowest model)
# and avoids re-shipping each deer's table to model-level workers. Each worker
# reads one deer, fits all candidate models, saves its rds, and returns its
# status + k/shrinkage diagnostics; the main process reduces them below.
process_deer <- function(i) {
  key <- keys[i]
  out_path <- sprintf("results/gam/results_gam_%s.rds", key)
  null_path <- sprintf("results/gam/results_gam_null_%s.rds", key)

  # Skip: both outputs already exist and we're not overwriting
  if (!overwrite && file.exists(out_path) && file.exists(null_path)) {
    return(list(status = "skip", key = key, audit = NULL))
  }

  # Wrapped so a per-deer error (data load, missing columns) doesn't halt the
  # run; fit_gam_mod itself already returns failure-status objects per model.
  tryCatch(
    {
      one_deer <- readRDS(track_files[i])
      gam_data <- prepare_gam_data(one_deer[[gam_input]][[1]])

      # Per-deer cyclic-tod basis cap: k must stay below the number of distinct
      # times of day this deer was fixed at (cc identifiability). Uses the RAW
      # distinct-tod count (~24-50, inflated by gaps/jitter), so it never binds
      # at K_TOD = 10; the meaningful tod resolution is ~12.
      n_tod <- length(unique(gam_data$tod_))
      k_tod <- max(3L, min(K_TOD, n_tod - 1L))

      # Null model first, saved on its own. It shares the data prep and k_tod
      # cap with the candidates but nothing else, so downstream code loads
      # exactly one reference model without indexing into the candidate list.
      results_gam_null <- fit_gam_mod(
        gam_data,
        make_null_formula(k_tod),
        select = SELECT
      )
      saveRDS(results_gam_null, null_path)

      # Season picks the slot-2/3 resource models (NDVI vs landcover; see
      # make_formulas). results_gam is named by model number ("1".."4").
      season <- strsplit(key, "_")[[1]][2]
      formulas <- make_formulas(k_tod, K_NDVI, season)

      results_gam <- purrr::map(
        formulas,
        function(f) fit_gam_mod(gam_data, f, select = SELECT)
      )

      saveRDS(results_gam, out_path)

      # Per-smooth k / shrinkage diagnostics for the end-of-run audit. imap's
      # index is the model-number name ("1".."4"); the null is tagged "null" so
      # the audit covers every fit this script produced.
      audit <- purrr::imap_dfr(
        c(list(null = results_gam_null), results_gam),
        function(r, model_no) {
          sd <- r$smooth_diag
          if (is.null(sd)) {
            return(NULL)
          }
          sd$key <- key
          sd$model <- model_no
          sd$k_tod <- k_tod
          sd
        }
      )

      list(status = "done", key = key, audit = audit)
    },
    error = function(e) {
      list(status = "fail", key = key, msg = conditionMessage(e), audit = NULL)
    }
  )
}

start_time <- Sys.time()

n_workers <- max(1L, parallel::detectCores() - 1L)
cat(sprintf(
  "Fitting %d deer on %d workers...\n",
  length(track_files),
  n_workers
))
future::plan(future::multisession, workers = n_workers)

results <- furrr::future_map(
  seq_along(track_files),
  process_deer,
  .options = furrr::furrr_options(packages = "mgcv", seed = TRUE),
  .progress = TRUE
)

future::plan(future::sequential)

# Reduce ----------------------------------------------------------------------
status <- vapply(results, function(r) r$status, character(1))
n_done <- sum(status == "done")
n_skipped <- sum(status == "skip")
n_failed <- sum(status == "fail")
audit_rows <- lapply(results[status == "done"], function(r) r$audit)

for (r in results[status == "fail"]) {
  cat(sprintf("[fail] %s: %s\n", r$key, r$msg))
}

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf(
  "Done: %d   Skipped: %d   Failed: %d   Elapsed: %.1f min\n",
  n_done,
  n_skipped,
  n_failed,
  elapsed
))

# k / shrinkage audit ---------------------------------------------------------
# Aggregated edf vs k' (and select=TRUE shrinkage outcomes) over every smooth
# fit this run. "k-bound" rows are the signal to raise k; "near-linear" /
# "removed" show where shrinkage simplified a smooth (e.g. a movement parameter
# with no time-of-day effect, or a habitat term with no support).
if (length(audit_rows)) {
  audit <- dplyr::bind_rows(audit_rows)
  saveRDS(audit, "results/gam/k_audit_gam.rds")

  cat(sprintf("\n=== k / shrinkage audit (%d smooth fits) ===\n", nrow(audit)))
  cat("status counts:\n")
  print(table(audit$status))

  kb <- audit[audit$status == "k-bound", , drop = FALSE]
  if (nrow(kb)) {
    cat(
      "\nk-bound smooths (edf >= 0.8 * k' -> consider higher k), by smooth:\n"
    )
    print(sort(table(kb$smooth), decreasing = TRUE))
  } else {
    cat(
      "\nNo smooth was k-bound: the chosen k values are ample dataset-wide.\n"
    )
  }

  simp <- audit[audit$status %in% c("near-linear", "removed"), , drop = FALSE]
  if (nrow(simp)) {
    cat(
      "\nselect = TRUE simplified these smooths (linear / removed):\n"
    )
    print(table(smooth = simp$smooth, status = simp$status))
  }
  cat("\nFull per-fit audit -> results/gam/k_audit_gam.rds\n")
}
