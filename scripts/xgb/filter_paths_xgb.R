#' @description
#' Apply the xgb path filters to every deer-year and deer model of one set:
#'   HR overlap   hr_overlap >= HR_MIN (score_sims_xgb.R)
#'   svf          svf_agree >= SVF_MIN, the variogram agreement
#'                (score_sims_xgb.R)
#'   log score    the one-step log score beats the deer null by at least
#'                LOGSCORE_MIN over the deer-year (logscore_xgb.R)
#' A model passes for a deer-year when it passes all three; a missing score
#' (a failed ctmm fit, say) fails its filter. The energy score and p_step in
#' the score table are not filters.
#'
#' Inputs:  filters/xgb/scores_xgb_<tag>.rds (score_sims_xgb.R) and
#'          filters/xgb/logscore_xgb_<tag>.rds (logscore_xgb.R); <tag> is
#'          xgb_scores_tag(SET, REACH)
#' Output:  filters/xgb/passed_xgb_<tag>.rds — one row per deer-year and
#'          deer model: the three scores, a pass flag for each filter and
#'          for all three; and a count printed per season
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
SET <- "train"
REACH <- "gamma99"
HR_MIN <- 0.8
SVF_MIN <- 0.8
LOGSCORE_MIN <- 3

# Load packages ---------------------------------------------------------------
library(tidyverse)

# helper functions
source("scripts/helper_functions.R")

tag <- xgb_scores_tag(SET, REACH)
scores <- readRDS(sprintf("filters/xgb/scores_xgb_%s.rds", tag)) |>
  filter(deer)
logscore <- readRDS(sprintf("filters/xgb/logscore_xgb_%s.rds", tag))
stopifnot(setequal(paste(scores$key, scores$config),
                   paste(logscore$key, logscore$config)))

passed <- scores |>
  select(key, season, config, best, simplest, hr_overlap, svf_agree) |>
  left_join(select(logscore, key, config, n_scored, delta_logp),
            by = c("key", "config")) |>
  mutate(
    pass_hr = !is.na(hr_overlap) & hr_overlap >= HR_MIN,
    pass_svf = !is.na(svf_agree) & svf_agree >= SVF_MIN,
    pass_logscore = !is.na(delta_logp) & delta_logp >= LOGSCORE_MIN,
    pass = pass_hr & pass_svf & pass_logscore
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
            simplest = sum(pass & simplest), best = sum(pass & best),
            either = n_distinct(key[pass]), .groups = "drop")
counts <- bind_rows(counts, summarise(counts, season = "all",
                                      across(-season, sum)))
cat(sprintf(paste0("=== %s: deer-years passing (hr, svf, logscore: with ",
                   "either deer model; simplest, best, either: all three) ",
                   "===\n"), tag))
print(as.data.frame(counts), row.names = FALSE)
cat(sprintf("-> %s\n", out_path))
