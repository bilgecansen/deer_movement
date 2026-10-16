#' @description
#' Score every deer-year's simulated paths against its real path, for each
#' final model: HR overlap, svf agreement and energy score. SET picks the
#' test year or the training years (sim_paths_xgb.R simulates both).
#'
#' The simulated paths step every 4 hours through the gaps in the real
#' fixes. Before scoring they are cut to the real fix times, so the real and
#' the simulated tracks are sampled alike. The first fix is where every path
#' starts, so the energy score leaves it out.
#'
#' The scores are the shared ones in helpers_metrics.R:
#'   * overlap_ud(): the Bhattacharyya overlap of the real deer's AKDE with
#'     the mean AKDE of the simulated paths (hr_overlap); 1 is best. It also
#'     returns the agreement of the ctmm models fitted to the real and the
#'     simulated paths, kept as svf_agree_ctmm for comparison only: the mean
#'     of the simulated paths' fits often does not describe those paths
#'   * variogram_agreement(): the agreement of the real path's exact
#'     variogram with the mean of the simulated paths' exact variograms, on
#'     log lag, up to half the real path's span (svf_agree); 1 is best
#'   * calc_energy_score(): at each real fix, the energy score of the real
#'     position against the simulated positions at that time, averaged over
#'     the fixes; in metres, lower is better
#'   * calc_p_step(): the share of the real steps longer than the energy
#'     score at the fix they start from (p_step); higher is better
#' Of these, HR overlap and svf agreement are path filters
#' (filter_paths_xgb.R, with the log score of logscore_xgb.R). The energy
#' score and p_step are kept in the table for reference: on paths that run
#' free from the first fix, they mostly reflect the home range after a day
#' or two.
#'
#' Deer-years are scored in parallel, one at a time per worker, the longest
#' first, and each is saved as it finishes, so an interrupted run resumes
#' where it stopped. A deer-year that fails is reported and left unscored.
#'
#' Inputs:  xgb_sims_dir(SET, REACH): sims/xgb_paths/ (test) or
#'          sims/xgb_paths_train/ (training), with a "_gamma99" suffix for
#'          that reach; sims_xgb_<key>.rds from sim_paths_xgb.R;
#'          data/tracks/data_<key>.rds
#' Output:  filters/xgb/scores_<tag>/scores_xgb_<key>.rds, one per
#'          deer-year, and filters/xgb/scores_xgb_<tag>.rds, all of them in
#'          one table; <tag> is xgb_scores_tag(SET, REACH), e.g. "train" or
#'          "train_gamma99"
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
# The deer-years to score: "test" or "train"
SET <- "train"
# The reach rule the paths were simulated with: "longest" or "gamma99"
REACH <- "gamma99"
# TRUE: score only the deer types
DEER_ONLY <- TRUE
# Deer-years to score: NULL for every simulated one not yet scored
KEYS <- NULL
N_PAR <- parallel::detectCores() - 1L

# Load packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(furrr)
})

# helper functions
source("scripts/helper_functions.R")

sims_dir <- xgb_sims_dir(SET, REACH)
tag <- xgb_scores_tag(SET, REACH)
out_dir <- sprintf("filters/xgb/scores_%s", tag)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (is.null(KEYS)) {
  files <- list.files(sims_dir, pattern = "^sims_xgb_.*[.]rds$",
                      full.names = TRUE)
  # Largest files first: the longest paths, which take longest to score
  files <- files[order(-file.size(files))]
  KEYS <- sub("^sims_xgb_(.*)[.]rds$", "\\1", basename(files))
  KEYS <- KEYS[!file.exists(file.path(out_dir,
                                      sprintf("scores_xgb_%s.rds", KEYS)))]
}
cat(sprintf("scoring %d %s deer-years on %d workers (%s)\n", length(KEYS),
            SET, N_PAR, format(Sys.time(), "%H:%M")))

score_one <- function(key) {
  source("scripts/helper_functions.R")
  sims <- readRDS(file.path(sims_dir, sprintf("sims_xgb_%s.rds", key)))
  obs <- readRDS(sprintf("data/tracks/data_%s.rds", key))$stp[[1]]
  # The real fixes: every step's start, and each burst's last end
  last_of_burst <- !duplicated(obs$burst_, fromLast = TRUE)
  fix_hours <- lubridate::round_date(
    sort(unique(c(obs$t1_, obs$t2_[last_of_burst]))), "hour"
  )
  obs_es <- obs[lubridate::round_date(obs$t1_, "hour") != min(fix_hours), ]
  configs <- names(sims$paths)
  if (DEER_ONLY) {
    configs <- configs[grepl("_deer$", configs)]
  }

  out <- purrr::map_dfr(configs, function(cf) {
    t0 <- Sys.time()
    sim <- sims$paths[[cf]]
    sim <- sim[lubridate::round_date(sim$t_, "hour") %in% fix_hours &
                 !is.na(sim$x_), ]
    ud <- tryCatch(
      suppressMessages(suppressWarnings(
        overlap_ud(obs, sim, n_sim = max(sim$nsim))
      )),
      error = function(e) {
        list(bat_uds = NA_real_, svf_agree = NA_real_,
             error = conditionMessage(e))
      }
    )
    svf <- tryCatch(
      suppressMessages(suppressWarnings(variogram_agreement(obs, sim))),
      error = function(e) structure(NA_real_, error = conditionMessage(e))
    )
    pick <- sims$picks[sims$picks$config == cf, ]
    es <- calc_energy_score(obs = obs_es, sim = sim)
    p_step <- calc_p_step(obs_es, sim)
    tibble::tibble(
      set = SET, reach_rule = REACH, reach = sims$reach,
      key = key, season = sims$season, config = cf,
      deer = pick$deer, best = pick$best, simplest = pick$simplest,
      hr_overlap = ud$bat_uds, svf_agree = as.numeric(svf),
      svf_agree_ctmm = ud$svf_agree,
      energy_score = es,
      p_step = p_step,
      ud_error = if (is.null(ud$error)) NA_character_ else ud$error,
      svf_error = if (is.null(attr(svf, "error"))) {
        NA_character_
      } else {
        attr(svf, "error")
      },
      n_fixes = length(fix_hours),
      seconds = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    )
  })
  saveRDS(out, file.path(out_dir, sprintf("scores_xgb_%s.rds", key)))
  "ok"
}

if (length(KEYS)) {
  future::plan(future::multisession, workers = min(N_PAR, length(KEYS)))
  status <- furrr::future_map_chr(
    KEYS,
    function(key) {
      tryCatch(score_one(key),
               error = function(e) paste("failed:", conditionMessage(e)))
    },
    .options = furrr::furrr_options(seed = TRUE, stdout = FALSE,
                                    chunk_size = 1L)
  )
  future::plan(future::sequential)
  failed <- status != "ok"
  if (any(failed)) {
    cat(sprintf("%d deer-years failed:\n", sum(failed)))
    cat(sprintf("  %s: %s\n", KEYS[failed], status[failed]), sep = "")
  }
}

scores <- purrr::map_dfr(
  list.files(out_dir, pattern = "^scores_xgb_.*[.]rds$", full.names = TRUE),
  readRDS
)
out_path <- sprintf("filters/xgb/scores_xgb_%s.rds", tag)
saveRDS(scores, out_path)
cat(sprintf("-> %s: %d deer-years, %d rows (%s)\n", out_path,
            dplyr::n_distinct(scores$key), nrow(scores),
            format(Sys.time(), "%H:%M")))
