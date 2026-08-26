#' @description
#' One-step-ahead log score for a single deer, computed for every numbered amt
#' model in results/amt/results_amt_<key>.rds AND for the GAM NULL model in
#' results/gam/results_gam_null_<key>.rds.
#'
#' There is no amt null. The amt models are measured against the GAM null, which
#' is scored here in the SAME run so that delta_logp is a like-for-like
#' difference over an identical step set. Only the numbered file is written:
#'   filters/amt/logscore_amt_<key>.rds   numbered models (+ delta_logp)
#'
#' ORDERING: fit_GAM.R must have run first, since the null comes from the GAM
#' side. The script stops with a clear message if the null file is missing.
#'
#' CAVEAT worth remembering: the GAM null carries a cyclic time-of-day smooth,
#' a more flexible movement block than the day/night factor the amt models use.
#' delta_logp therefore mixes "does habitat help" with "my movement block is
#' less flexible than the reference's", and the amt models are handicapped by
#' roughly the movement-block gain (~0.017-0.066 log units per step).
#'
#' Usage: Rscript run_logscore_amt.R <id> <season> <year>
#'   id     — deer ID
#'   season — season string (e.g. "fa", "nb")
#'   year   — year (integer)

# Parse command line arguments -------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste0(
      "Usage: Rscript run_logscore.R <id> <season> <year>\n",
      "Example: Rscript run_logscore.R 5000 fa 2017"
    )
  )
}

id <- args[1]
season <- args[2]
year <- as.integer(args[3])

key <- sprintf("%s_%s_%d", id, season, year)
cat(sprintf("Running log score for deer %s\n", key))

# Load packages ----------------------------------------------------------------
library(amt)
library(terra)
library(tidyverse)
library(foreach)
library(sf)
library(furrr)
library(parallel)

# helper functions
source("scripts/helper_functions.R")

# Load data --------------------------------------------------------------------
start_time <- Sys.time()

# Deer movement data — single-row per-deer file
deer_mvt <- readRDS(sprintf("data/tracks/data_%s.rds", key))

# landscape data — season-specific annual landcover (categorical band + per-
# class binary indicator layers). env_old terrain covariates are dropped.
env_raster <- load_landcover(year, season)

# NDVI data
ndvi_year <- load_ndvi(year)

# amt models
results_amt <- readRDS(sprintf("results/amt/results_amt_%s.rds", key))

# The reference model comes from the GAM side: there is no amt null.
null_path <- sprintf("results/gam/results_gam_null_%s.rds", key)
if (!file.exists(null_path)) {
  stop(
    "GAM null not found: ", null_path, "\n",
    "amt models are scored against the GAM null, so fit_GAM.R must run first."
  )
}
results_gam_null <- readRDS(null_path)
null_gam <- if (is.character(results_gam_null$gam)) {
  NULL
} else {
  results_gam_null$gam
}

# The GAM kernel needs the tentative movement distributions; they ride along on
# stp.var as attributes.
sl_distr <- attr(deer_mvt$stp.var[[1]], "sl_")
ta_distr <- attr(deer_mvt$stp.var[[1]], "ta_")

# Score every fitted model — null/failed models filtered out automatically
# inside the precompute step (returns NULL for those).
n_models <- length(results_amt)
models_to_run <- seq_len(n_models)

# Pre-crop rasters for this single deer ---------------------------------------
crop_extent <- sf::st_buffer(
  sf::st_as_sf(
    deer_mvt$stp[[1]],
    coords = c('x1_', 'y1_'),
    crs = 6610
  ),
  CROP_BUFFER_M
)

env_cropped <- terra::crop(env_raster, crop_extent)

# HR rasters — no HR_bin (no model uses it). HR_edge stays in metres; for
# HR_center we also build the log1p transformed layer the formulas reference.
env_cropped$HR_edge <- load_hr_edge_raster(id, season, year, env_cropped)
env_cropped$HR_center <- load_hr_center_raster(id, season, year, env_cropped)
env_cropped$HR_center_log <- log1p(env_cropped$HR_center)

deer_input <- list(
  crop_env = terra::wrap(env_cropped),
  crop_ndvi = terra::wrap(terra::crop(ndvi_year, crop_extent)),
  stp = deer_mvt$stp[[1]]
)

# Precompute simulation models ------------------------------------------------
cat("Precomputing simulation models...\n")

model_sims <- purrr::map(models_to_run, function(m) {
  iss_i <- results_amt[[m]]$iss

  if (is.character(iss_i)) {
    return(NULL)
  }

  train_i <- deer_mvt$stp.var[[1]]

  coefs <- iss_i$model$coefficients
  names(coefs) <- rename_landcover_coefs(names(coefs))
  coefs <- coefs[!is.na(coefs)]

  dummy_sim <- amt::make_issf_model(coefs = coefs)
  mm_names <- colnames(model.matrix(
    amt:::ssf_formula(dummy_sim$model$formula),
    data = train_i
  ))

  for (idx in seq_along(coefs)) {
    if (grepl(":", names(coefs)[idx]) && !(names(coefs)[idx] %in% mm_names)) {
      parts <- strsplit(names(coefs)[idx], ":")[[1]]
      perms <- combinat::permn(parts)
      for (p in perms) {
        candidate <- paste(p, collapse = ":")
        if (candidate %in% mm_names) {
          names(coefs)[idx] <- candidate
          break
        }
      }
    }
  }

  amt::make_issf_model(
    coefs = coefs,
    sl = iss_i$sl_,
    ta = iss_i$ta_
  )
})
names(model_sims) <- as.character(models_to_run)

# Free large objects
rm(env_raster, ndvi_year, results_amt, results_gam_null, deer_mvt)
gc()

# Run log score across the null and every numbered model ---------------------
# The null is carried as a named unit ("null") rather than a model number so it
# can never be mistaken for one. It is a GAM, so it is scored with
# onestep_logscore_gam(); the numbered models are amt fits scored with
# onestep_logscore(). Both are evaluated over the same observed steps.
cat("Computing one-step log scores...\n")

score_units <- c(list(null = null_gam), model_sims)
units_to_run <- names(score_units)

future::plan(
  multisession,
  workers = min(length(units_to_run), parallel::detectCores() - 1)
)

results_logscore <- furrr::future_map(
  units_to_run,
  function(m) {
    cat("  Model:", m, "\n")

    model_obj <- score_units[[m]]

    if (is.null(model_obj)) {
      return(data.frame(
        burst_ = NA, step_index = NA_integer_,
        t1_ = as.POSIXct(NA), logp = NA_real_,
        status = "model_failed_to_fit", model = m
      ))
    }

    env_local <- terra::unwrap(deer_input$crop_env)
    ndvi_local <- terra::unwrap(deer_input$crop_ndvi)

    ll_df <- if (identical(m, "null")) {
      onestep_logscore_gam(
        stp_data = deer_input$stp,
        env_test = env_local,
        ndvi_test = ndvi_local,
        gam_train = model_obj,
        sl_distr = sl_distr,
        ta_distr = ta_distr
      )
    } else {
      onestep_logscore(
        stp_data = deer_input$stp,
        env_test = env_local,
        ndvi_test = ndvi_local,
        amt_train = model_obj
      )
    }

    # Return the PER-STEP scores. Totals are formed in the main process once the
    # step set common to every model is known -- summing here would let each
    # model total over its own step set.
    ll_df$model <- m
    ll_df
  },
  .options = furrr::furrr_options(
    packages = c(
      "mgcv", "amt", "terra", "sf", "dplyr", "lubridate", "foreach"
    ),
    stdout = TRUE,
    seed = TRUE
  )
)

future::plan(sequential)
gc()

# Total over the COMMON step set ---------------------------------------------
# Every model is totalled over exactly the steps that ALL of them scored -- the
# null included -- so n_steps is identical across models by construction and
# delta_logp is a like-for-like difference rather than one that happens to line
# up.
#
# Two quite different things produce an unscored step, and they are counted
# separately rather than both vanishing into na.rm = TRUE:
#   skipped_no_heading  first step of a burst; nothing precedes it, so no
#                       incoming heading exists. By design, ~19% of steps.
#   failed_*            a real failure, overwhelmingly failed_outside_disc: the
#                       observed step was longer than the candidate disc's
#                       radius. Routine, and previously invisible.
per_step <- dplyr::bind_rows(results_logscore)

scored <- per_step %>% dplyr::filter(!is.na(logp))
n_units <- dplyr::n_distinct(per_step$model)
common <- scored %>%
  dplyr::count(burst_, step_index) %>%
  dplyr::filter(n == n_units) %>%
  dplyr::mutate(k = paste(burst_, step_index)) %>%
  dplyr::pull(k)

totals <- per_step %>%
  dplyr::mutate(k = paste(burst_, step_index)) %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    total_logp = if (all(is.na(logp))) {
      NA_real_
    } else {
      sum(logp[k %in% common])
    },
    n_steps = sum(k %in% common & !is.na(logp)),
    n_scored_alone = sum(!is.na(logp)),
    n_skipped_no_heading = sum(status == "skipped_no_heading"),
    n_failed = sum(grepl("^failed_", status)),
    n_dropped_for_common = sum(!is.na(logp) & !(k %in% common)),
    .groups = "drop"
  )

cat(sprintf(
  "Steps: %d observed | %d scored by every model (common set)\n",
  dplyr::n_distinct(paste(per_step$burst_, per_step$step_index)),
  length(common)
))
if (any(totals$n_failed > 0)) {
  cat("Failure breakdown (non-heading):\n")
  print(per_step %>%
    dplyr::filter(grepl("^failed_", status)) %>%
    dplyr::count(model, status))
}

# delta_logp against the GAM null. The null's own scores are NOT written -- they
# belong to the GAM path, which saves them itself.
null_results <- totals %>% dplyr::filter(model == "null")
null_logp <- null_results$total_logp

results <- totals %>%
  dplyr::filter(model != "null") %>%
  dplyr::mutate(model = as.integer(model)) %>%
  dplyr::arrange(model) %>%
  dplyr::mutate(delta_logp = total_logp - null_logp)

cat("GAM null (reference, not saved here):\n")
print(null_results)
cat("Numbered amt models:\n")
print(results)

# Save -------------------------------------------------------------------------
dir.create("filters/amt", showWarnings = FALSE, recursive = TRUE)
saveRDS(results, sprintf("filters/amt/logscore_amt_%s.rds", key))

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("Deer %s completed in %.1f minutes\n", key, elapsed))
