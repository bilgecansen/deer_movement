#' @description
#' Apply the xgb path filters to every deer-year and deer model of one set:
#'   HR overlap   hr_overlap >= HR_MIN (score_sims_xgb.R)
#'   svf          svf_agree >= SVF_MIN, the variogram agreement
#'                (score_sims_xgb.R)
#'   log score    the one-step log score beats the null by at least
#'                LOGSCORE_MIN over the deer-year. The null is the
#'                deer-year's own GAM null (LOGSCORE_NULL "null_gam",
#'                logscore_gam_xgb.R) or the xgb deer null ("null_deer",
#'                logscore_xgb.R)
#'   habitat      the one-step log score beats the deer-year's own
#'                Gilbertson et al. (2025) habitat model by at least
#'                HABITAT_MIN (logscore_gam_xgb.R, GAM "gilbertson_gam")
#' A model passes for a deer-year when it passes all four; a missing score
#' (a failed ctmm fit, or a habitat model that did not fit) fails its
#' filter. The energy score and p_step in the score table are not filters.
#'
#' Inputs:  filters/xgb/scores_xgb_<tag>.rds (score_sims_xgb.R),
#'          filters/xgb/logscore_null_gam_<tag>.rds or
#'          filters/xgb/logscore_xgb_<tag>.rds, by LOGSCORE_NULL, and
#'          filters/xgb/logscore_gilbertson_gam_<tag>.rds; <tag> is
#'          xgb_scores_tag(SET, REACH)
#' Output:  filters/xgb/passed_xgb_<tag>.rds — one row per deer-year and
#'          deer model: the four scores, a pass flag for each filter and
#'          for all four; and a count printed per season
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
SET <- "train"
REACH <- "gamma99"
HR_MIN <- 0.8
SVF_MIN <- 0.8
LOGSCORE_MIN <- 3
# The null the log score must beat: "null_gam" or "null_deer"
LOGSCORE_NULL <- "null_gam"
HABITAT_MIN <- 3

# Load packages ---------------------------------------------------------------
library(tidyverse)

# helper functions
source("scripts/helper_functions.R")

tag <- xgb_scores_tag(SET, REACH)
scores <- readRDS(sprintf("filters/xgb/scores_xgb_%s.rds", tag)) |>
  filter(deer)
logscore_path <- switch(
  LOGSCORE_NULL,
  null_gam = sprintf("filters/xgb/logscore_null_gam_%s.rds", tag),
  null_deer = sprintf("filters/xgb/logscore_xgb_%s.rds", tag),
  stop("LOGSCORE_NULL must be 'null_gam' or 'null_deer'")
)
logscore <- readRDS(logscore_path)
habitat <- readRDS(sprintf("filters/xgb/logscore_gilbertson_gam_%s.rds",
                           tag))
stopifnot(setequal(paste(scores$key, scores$config),
                   paste(logscore$key, logscore$config)))
# A deer-year whose habitat model did not fit has no habitat score
stopifnot(all(paste(habitat$key, habitat$config) %in%
                paste(scores$key, scores$config)))

passed <- scores |>
  select(key, season, config, best, simplest, hr_overlap, svf_agree) |>
  left_join(select(logscore, key, config, n_scored, delta_logp),
            by = c("key", "config")) |>
  left_join(select(habitat, key, config, n_scored_habitat = n_scored,
                   delta_logp_habitat = delta_logp),
            by = c("key", "config")) |>
  mutate(
    logscore_null = LOGSCORE_NULL,
    pass_hr = !is.na(hr_overlap) & hr_overlap >= HR_MIN,
    pass_svf = !is.na(svf_agree) & svf_agree >= SVF_MIN,
    pass_logscore = !is.na(delta_logp) & delta_logp >= LOGSCORE_MIN,
    pass_habitat = !is.na(delta_logp_habitat) &
      delta_logp_habitat >= HABITAT_MIN,
    pass = pass_hr & pass_svf & pass_logscore & pass_habitat
  ) |>
  arrange(season, key, config)
out_path <- sprintf("filters/xgb/passed_xgb_%s.rds", tag)
saveRDS(passed, out_path)

counts <- passed |>
  group_by(season) |>
  summarise(deer_years = n_distinct(key),
            hr = n_distinct(key[pass_hr]),
            svf = n_distinct(key[pass_svf]),
            logscore = n_distinct(key[pass_logscore]),
            habitat = n_distinct(key[pass_habitat]),
            simplest = sum(pass & simplest), best = sum(pass & best),
            either = n_distinct(key[pass]), .groups = "drop")
counts <- bind_rows(counts, summarise(counts, season = "all",
                                      across(-season, sum)))
cat(sprintf(paste0("=== %s, log score vs %s: deer-years passing (hr, svf, ",
                   "logscore, habitat: with either deer model; simplest, ",
                   "best, either: all four) ===\n"), tag, LOGSCORE_NULL))
print(as.data.frame(counts), row.names = FALSE)
cat(sprintf("-> %s\n", out_path))
