#' @description
#' Score every deer-year's simulated paths against its real path, for each
#' model: HR overlap, svf agreement and energy score. SET picks the test year
#' or the training years. MODELS picks the paths: "xgb", the final xgb models
#' (sim_paths_xgb.R), or "gam", the GAM null and the Gilbertson model
#' (sim_paths_gam_xgb.R), which carry no xgb pick, so their deer, best and
#' simplest columns are NA.
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
#' BATCH scores only that many of the deer-years left, the longest first,
#' then writes the table and stops; run_scores_xgb.sh runs batch after batch
#' in a fresh R each time, so a ctmm fit that hangs costs one batch.
#'
#' Inputs:  xgb_sims_dir(SET, REACH): sims/xgb_paths/ (test) or
#'          sims/xgb_paths_train/ (training), with a "_gamma99" suffix for
#'          that reach; sims_<MODELS>_<key>.rds; data/tracks/data_<key>.rds
#' Output:  per deer-year, filters/xgb/scores_<tag>/scores_xgb_<key>.rds
#'          (xgb) or filters/xgb/scores_gam_<tag>/scores_gam_<key>.rds
#'          (gam); and all of them in one table,
#'          filters/xgb/scores_<MODELS>_<tag>.rds; <tag> is
#'          xgb_scores_tag(SET, REACH), e.g. "train" or "train_gamma99"
#'
#' Configuration: edit the block below before running. The set, the models
#' and the batch can also be given on the command line:
#'   Rscript scripts/xgb/score_sims_xgb.R train gam 50

# Configuration ---------------------------------------------------------------
# The deer-years to score: "test" or "train"
SET <- "train"
# The paths to score: "xgb" or "gam"
MODELS <- "xgb"
# Deer-years to score in this run, the longest left first: NULL for all
BATCH <- NULL
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1) {
  SET <- .args[1]
}
if (length(.args) >= 2) {
  MODELS <- .args[2]
}
if (length(.args) >= 3) {
  BATCH <- as.integer(.args[3])
}
stopifnot(SET %in% c("test", "train"), MODELS %in% c("xgb", "gam"))
# The reach rule the paths were simulated with: "longest" or "gamma99"
REACH <- "gamma99"
# TRUE: score only the xgb deer types
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
out_dir <- if (MODELS == "xgb") {
  sprintf("filters/xgb/scores_%s", tag)
} else {
  sprintf("filters/xgb/scores_gam_%s", tag)
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
sims_file <- function(key) {
  file.path(sims_dir, sprintf("sims_%s_%s.rds", MODELS, key))
}
score_file <- function(key) {
  file.path(out_dir, sprintf("scores_%s_%s.rds", MODELS, key))
}
if (is.null(KEYS)) {
  files <- list.files(sims_dir, pattern = sprintf("^sims_%s_.*[.]rds$",
                                                  MODELS),
                      full.names = TRUE)
  # Largest files first: the longest paths, which take longest to score
  files <- files[order(-file.size(files))]
  KEYS <- sub(sprintf("^sims_%s_(.*)[.]rds$", MODELS), "\\1",
              basename(files))
  KEYS <- KEYS[!file.exists(score_file(KEYS))]
}
n_left <- length(KEYS)
if (!is.null(BATCH)) {
  KEYS <- utils::head(KEYS, BATCH)
}
cat(sprintf("scoring %d of %d %s %s deer-years left on %d workers (%s)\n",
            length(KEYS), n_left, SET, MODELS, N_PAR,
            format(Sys.time(), "%H:%M")))

score_one <- function(key) {
  source("scripts/helper_functions.R")
  sims <- readRDS(sims_file(key))
  obs <- readRDS(sprintf("data/tracks/data_%s.rds", key))$stp[[1]]
  # The real fixes: every step's start, and each burst's last end
  last_of_burst <- !duplicated(obs$burst_, fromLast = TRUE)
  fix_hours <- lubridate::round_date(
    sort(unique(c(obs$t1_, obs$t2_[last_of_burst]))), "hour"
  )
  obs_es <- obs[lubridate::round_date(obs$t1_, "hour") != min(fix_hours), ]
  configs <- names(sims$paths)
  if (MODELS == "xgb" && DEER_ONLY) {
    configs <- configs[grepl("_deer$", configs)]
  }
  # The GAM paths carry no xgb pick
  picks <- if (MODELS == "xgb") {
    sims$picks
  } else {
    tibble::tibble(config = configs, deer = NA, best = NA, simplest = NA)
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
    pick <- picks[picks$config == cf, ]
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
  saveRDS(out, score_file(key))
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
  list.files(out_dir, pattern = sprintf("^scores_%s_.*[.]rds$", MODELS),
             full.names = TRUE),
  readRDS
)
out_path <- sprintf("filters/xgb/scores_%s_%s.rds", MODELS, tag)
saveRDS(scores, out_path)
cat(sprintf("-> %s: %d deer-years, %d rows (%s)\n", out_path,
            dplyr::n_distinct(scores$key), nrow(scores),
            format(Sys.time(), "%H:%M")))
