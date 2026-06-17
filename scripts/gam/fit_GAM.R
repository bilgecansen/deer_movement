#' @description
#' Fit deer movement GAMs (SSFs with penalised smooths) for every deer with a
#' wrangled track file in data/tracks/. This is the GAM analogue of fit_issf.R.
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
#' For each (id, season, year):
#'   * If the output rds already exists and `overwrite` is FALSE, skip.
#'   * Otherwise read the wrangled data, reshape the chosen design's covariate
#'     table (`gam_input`) into mgcv cox.ph form, fit all candidate formulas
#'     with fit_gam_mod() (which returns failure-status objects rather than
#'     erroring per-model), and save the list to results/results_gam_<key>.rds.
#' Errors are caught per deer so a single failure does not halt the loop.
#'
#' Random points: two designs are wrangled per deer.
#'   * "stp.var.nonp" — random steps sampled uniformly over a disc of radius
#'     ~ max step length (wrangle_deer_mvt.R, model = "nonp"). Required for
#'     *non-parametric* movement smooths s(sl_) / s(ta_) (Klappstein et al.
#'     2024, Section 3.1); also valid for parametric kernels.
#'   * "stp.var" — gamma / von Mises (iSSF-style). Usable for GAMs with a
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
#         identifiable. With 4-hour fixes distinct tod ~12 (min ~10), so the cap
#         is normally a no-op at 8 (see the MOVE comment below).
# K_NDVI  basis dimension (per landcover class) of the NDVI 'fs' smooths. NDVI is
#         continuous with rich support, so it is NOT tod-constrained; the k-audit
#         reports whether it is ever the binding constraint (it has wide headroom
#         at 5, so bumping it is unnecessary unless the audit says otherwise).
# SELECT  double-penalty shrinkage (Marra & Wood 2011). FALSE for the production
#         fit: with a fixed, deliberate model set, between-model filtering is the
#         selection step, so we don't also want per-deer within-model term removal
#         (it makes "model 4" a different effective structure on different deer).
#         REML smoothing penalties still regularise each fit; dropping the extra
#         null-space penalty is also faster and converges more reliably. (TRUE was
#         useful during exploration.)
K_TOD <- 8L
K_NDVI <- 5L
SELECT <- FALSE

# Load packages ---------------------------------------------------------------
library(mgcv)
library(tidyverse)
library(furrr)

# helper functions
source("scripts/helper_functions.R")

# Formulas --------------------------------------------------------------------
# make_formulas() builds the candidate RHS vector for a given cyclic-tod basis
# (k_tod, capped per deer) and NDVI-by-landcover 'fs' basis (k_ndvi). fit_gam_mod
# prepends the "cbind(times, stratum) ~" Cox-PH response. Each model is the GAM
# translation of the corresponding iSSF formula in fit_issf.R (numbered 1-6 here;
# model 6 is fit_issf.R's model 8, which still uses the older 1-5,8 numbering).
#
# MOVE = parametric gamma / von Mises movement kernel. The three movement
# covariates enter both as parametric main effects (sl_ + log(sl_) + cos(ta_))
# AND as cyclic-spline interactions with time of day (zebra model; Klappstein
# et al. 2024). The main effects carry the baseline correction to the tentative
# kernel (sl_/log(sl_) -> gamma rate/shape, cos(ta_) -> von Mises concentration);
# the by= smooths add the (shrinkable) time-of-day modulation.
# k_tod caps wiggliness below each deer's distinct-tod count
# (4-hour fixes -> ~12 positions, min ~10). Landcover interactions use a global
# smooth plus hierarchical "fs" deviations (Pedersen et al. 2019 "Model GS":
# one shrunk curve per class around a shared global response) rather than
# independent by= smooths, which are unidentifiable when a deer visits only a few
# classes. The global smooth mirrors the iSSF main effect and lets sparse classes
# pool toward the shared response instead of toward zero.
#
# NDVI (models 4 & 5) is biologically meaningful only on the open-canopy classes
# in NDVI_VEG_CLASSES (corn / soybeans / alfalfa_hay / small_grains / other_ag /
# grassland); on forest, wetlands and developed it is noise. So its hierarchical
# block is the GS structure RESTRICTED to those classes: the numeric 0/1 is_veg
# indicator (prepare_gam_data) is passed as `by=` to both the global veg smooth
# s(ndvi_veg, by = is_veg) and the per-class deviations
# s(ndvi_veg, veg_class, bs = 'fs', by = is_veg), zeroing each smooth exactly on
# non-veg steps (by= multiplies the whole smooth). veg_class carries only the 6
# veg levels (non-veg rows parked on level 1, made inert by is_veg = 0). The
# class-level intercept stays a single random effect over ALL classes,
# s(wiscland_end, bs = 're') = alpha_land[k] ~ N(0, sigma_land) (the SSF has no
# global intercept, so its mean is absorbed by the baseline hazard). HR_center
# (model 6) keeps the all-class fs (hc_fs); it is meaningful everywhere.

make_formulas <- function(k_tod, k_ndvi, season) {
  move <- sprintf(
    paste0(
      "sl_ + log(sl_) + cos(ta_) + ",
      "s(tod_, bs = 'cc', k = %1$d, by = sl_)"
    ),
    k_tod
  )
  hc_fs <- sprintf("s(HR_center_end, wiscland_end, bs = 'fs', k = %d)", k_ndvi)

  # Always-fit structural models
  f1 <- move # 1 movement only
  f2 <- paste(move, "+ s(HR_edge_end)")
  f3 <- paste(move, "+ s(HR_center_end)")
  f6 <- paste(move, "+ s(HR_center_end) +", hc_fs)

  # Resource models for slots 4 & 5 — season-dependent. Non-breeding (winter)
  # deer have no valid MODIS NDVI (snow / dormancy -> all-or-mostly NA), so the
  # NDVI models (4 & 5) can't be fit; we substitute the landcover-only models
  # and record them in slots 4 & 5. Either way the object holds 6 models
  # and slots 4/5 always mean "the resource-selection model".
  if (season == "nb") {
    # Winter NDVI is unusable; the landcover-only substitutes use a penalised
    # random effect on landcover, s(wiscland_end, bs = 're'), rather than a
    # fixed factor. It shrinks sparse classes and avoids the near-separation
    # that makes the raw factor slow / non-convergent in cox.ph (validated:
    # ~17x faster, no step failures).
    f4 <- paste(move, "+ s(HR_center_end) + s(wiscland_end, bs = 're')")
    f5 <- paste(move, "+ s(wiscland_end, bs = 're')")
  } else {
    # Hierarchical NDVI block, restricted to NDVI_VEG_CLASSES (see the long
    # comment above and prepare_gam_data for is_veg / ndvi_veg / veg_class): an
    # all-class landcover random intercept, a shared veg NDVI response, and
    # per-veg-class 'fs' deviations -- the latter two switched off on non-veg
    # steps via by = is_veg.
    ndvi_block <- sprintf(
      paste0(
        "s(wiscland_end, bs = 're') + ",
        "s(ndvi_veg, by = is_veg) + ",
        "s(ndvi_veg, veg_class, bs = 'fs', k = %d, by = is_veg)"
      ),
      k_ndvi
    )
    f4 <- paste(move, "+ s(HR_center_end) +", ndvi_block)
    f5 <- paste(move, "+", ndvi_block)
  }

  stats::setNames(c(f1, f2, f3, f4, f5, f6), c("1", "2", "3", "4", "5", "6"))
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
dir.create("results", showWarnings = FALSE)

# Loop (parallel across deer) -------------------------------------------------
# fit_GAM is low-memory tabular work (no rasters), so deer are independent,
# CPU-bound units; parallelising across deer gives the best load-balancing (a
# free worker grabs the next deer rather than waiting on a deer's slowest model)
# and avoids re-shipping each deer's table to model-level workers. Each worker
# reads one deer, fits all candidate models, saves its rds, and returns its
# status + k/shrinkage diagnostics; the main process reduces them below.
process_deer <- function(i) {
  key <- keys[i]
  out_path <- sprintf("results/results_gam_%s.rds", key)

  # Skip: output already exists and we're not overwriting
  if (!overwrite && file.exists(out_path)) {
    return(list(status = "skip", key = key, audit = NULL))
  }

  # Wrapped so a per-deer error (data load, missing columns) doesn't halt the
  # run; fit_gam_mod itself already returns failure-status objects per model.
  tryCatch(
    {
      one_deer <- readRDS(track_files[i])
      gam_data <- prepare_gam_data(one_deer[[gam_input]][[1]])

      # Per-deer cyclic-tod basis cap: k must stay below the number of distinct
      # times of day this deer was fixed at (cc identifiability). Normally a
      # no-op at K_TOD = 8 (distinct tod ~12; min ~10 across the dataset).
      n_tod <- length(unique(gam_data$tod_))
      k_tod <- max(3L, min(K_TOD, n_tod - 1L))

      # Season picks the slot-4/5 resource models (NDVI vs landcover; see
      # make_formulas). results_gam is named by model number ("1".."6").
      season <- strsplit(key, "_")[[1]][2]
      formulas <- make_formulas(k_tod, K_NDVI, season)

      results_gam <- purrr::map(
        formulas,
        function(f) fit_gam_mod(gam_data, f, select = SELECT)
      )

      saveRDS(results_gam, out_path)

      # Per-smooth k / shrinkage diagnostics for the end-of-run audit. imap's
      # index is the model-number name ("1".."6").
      audit <- purrr::imap_dfr(results_gam, function(r, model_no) {
        sd <- r$smooth_diag
        if (is.null(sd)) {
          return(NULL)
        }
        sd$key <- key
        sd$model <- model_no
        sd$k_tod <- k_tod
        sd
      })

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
  saveRDS(audit, "results/k_audit_gam.rds")

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
      "\nselect = TRUE simplified these smooths (collapsed to linear / removed):\n"
    )
    print(table(smooth = simp$smooth, status = simp$status))
  }
  cat("\nFull per-fit audit -> results/k_audit_gam.rds\n")
}
