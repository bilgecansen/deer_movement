#' @description
#' GAM analogue of scripts/amt/filter_models.R. Combine per-deer GAM filter
#' outputs (udoverlap_gam, logscore_gam, energy score es_gam) plus the observed
#' step-length distribution (from the track), into one long (deer x model) data
#' frame, and apply four sequential model-selection gates.
#'
#' Inputs:
#'   filters/gam/udoverlap_gam_<key>.rds — list keyed by model: list(bat_uds,
#'     svf_agree)
#'   filters/gam/logscore_gam_<key>.rds — df: model, total_logp, n_steps,
#'     delta_logp
#'   filters/gam/es_gam.rds — df: id, season, year, model, energy_score
#'   data/tracks/data_<key>.rds       — observed steps (stp) -> observed sl_
#'
#' Outputs:
#'   filters/gam/filter_combined_gam.rds  — every (deer, model) with all metrics
#'   filters/gam/filter_selected_gam.rds  — rows surviving all four gates
#'   plots/filter_violins_pre_gam.png, plots/filter_violins_post_gam.png
#'
#' ONLY THE NUMBERED MODELS ARE FILTERED. The null model (movement +
#' s(HR_center_end)) is fit and log-scored separately and does not go through
#' the gates: it is the reference the numbered models are measured against, not
#' a candidate competing with them. Concretely —
#'   * it is never simulated, so it has no bat_uds / svf_agree / energy_score;
#'   * its log score lives in logscore_gam_null_<key>.rds, which this script
#'     never reads (delta_logp, already differenced against it, is enough);
#'   * there is consequently no "null model as fallback" branch at gate 3 —
#'     a deer with no passing numbered model simply contributes nothing.
#'
#' GAM model set is numbered 1..4 (see fit_GAM.R): 1 movement only, 2 & 3
#' resource (HR-center + NDVI x landcover / NDVI x landcover alone;
#' landcover-only substitutes in winter), 4 HR-center x landcover. The
#' environmental (resource-selection) models are 2, 3, 4; model 1
#' is movement-only.
#'
#' Gates (applied sequentially — only survivors flow into the next step):
#'   1. bat_uds    >= 0.8
#'   2. svf_agree  >= 0.8
#'   3. delta_logp >= 3 (must BEAT the null by this margin — see delta_logp_min)
#'   4. p_excd     >= 0.5  P(observed sl > energy_score); equivalent to ES below
#'                         the deer's median step length.

library(tidyverse)

# Pipeline mode ---------------------------------------------------------------
# FALSE: use in-sample filter outputs (filters/gam/udoverlap_gam_*.rds,
#        filters/gam/logscore_gam_*.rds, filters/gam/es_gam.rds) and write
#        filter_combined_gam.rds / filter_selected_gam.rds + non-suffixed plots.
# TRUE: use out-of-sample (test) filter outputs
#   (filters/gam/udoverlap_gam_test_*.rds,
#        filters/gam/logscore_gam_test_*.rds, filters/gam/es_gam_test.rds) and
#        write filter_combined_gam_test.rds / filter_selected_gam_test.rds +
#        _test plots.
test_mode <- F

# delta_logp_min is a *positive* margin: a numbered model must beat the null by
# at least this much. Under the previous HR-edge null the threshold was -3
# ("don't lose to the null by more than 3"), which is the wrong stance now that
# the null is the HR-center model — the same reference used in the earlier
# null-2 vs null-3 comparison, where null 3 was paired with +3.
delta_logp_min <- 3
ud_min <- 0.8
svf_min <- 0.8
p_excd_min <- 0.5

# File-naming derived from test_mode -----------------------------------------
suffix <- if (test_mode) "_test" else ""
udov_prefix <- sprintf("udoverlap_gam%s_", suffix)
logs_prefix <- sprintf("logscore_gam%s_", suffix)
es_file <- sprintf("filters/gam/es_gam%s.rds", suffix)
combined_out <- sprintf("filters/gam/filter_combined_gam%s.rds", suffix)
selected_out <- sprintf("filters/gam/filter_selected_gam%s.rds", suffix)
plot_pre <- sprintf("plots/filter_violins_pre_gam%s.png", suffix)
plot_post <- sprintf("plots/filter_violins_post_gam%s.png", suffix)

# Discover keys ---------------------------------------------------------------
# Broad listing matches both test and non-test files (the regex
# `^udoverlap_gam_.*` happily eats `udoverlap_gam_test_…`); narrow to the
# requested mode explicitly. The `_gam_` infix keeps amt outputs
# (udoverlap_<key>, logscore_<key>) out of this GAM pipeline entirely.
udov_files <- list.files(
  "filters/gam",
  "^udoverlap_gam_.*\\.rds$",
  full.names = TRUE
)
logs_files <- list.files(
  "filters/gam",
  "^logscore_gam_.*\\.rds$",
  full.names = TRUE
)

# Drop the null-model log scores in BOTH modes. They match the broad regex
# (`logscore_gam_null_<key>` and `logscore_gam_null_test_<key>`) but are not
# candidates: the null does not go through the gates, and the comparison against
# it already rides on the numbered rows as delta_logp. Without this they would
# be read as deer with keys like "null_5004_fa_2017".
logs_files <- logs_files[!grepl("^logscore_gam_null_", basename(logs_files))]

if (test_mode) {
  udov_files <- udov_files[grepl("^udoverlap_gam_test_", basename(udov_files))]
  logs_files <- logs_files[grepl("^logscore_gam_test_", basename(logs_files))]
} else {
  udov_files <- udov_files[!grepl("^udoverlap_gam_test_", basename(udov_files))]
  logs_files <- logs_files[!grepl("^logscore_gam_test_", basename(logs_files))]
}

keys_udov <- gsub(
  sprintf("^%s(.*)\\.rds$", udov_prefix),
  "\\1",
  basename(udov_files)
)
keys_logs <- gsub(
  sprintf("^%s(.*)\\.rds$", logs_prefix),
  "\\1",
  basename(logs_files)
)
keys <- intersect(keys_udov, keys_logs)

cat(sprintf(
  "Found %d deer with udoverlap and logscore (udov: %d, logs: %d)\n",
  length(keys),
  length(keys_udov),
  length(keys_logs)
))

parse_key <- function(k) {
  parts <- strsplit(k, "_")[[1]]
  tibble::tibble(
    id = parts[1],
    season = parts[2],
    year = as.integer(parts[3])
  )
}

# Load udoverlap + logscore per deer ------------------------------------------
per_deer_df <- purrr::map_dfr(keys, function(k) {
  ud <- readRDS(file.path("filters/gam", sprintf("%s%s.rds", udov_prefix, k)))
  ls <- readRDS(file.path("filters/gam", sprintf("%s%s.rds", logs_prefix, k)))

  ud_df <- purrr::imap_dfr(ud, function(x, m) {
    if (length(x) == 1 && is.na(x)) {
      tibble::tibble(
        model = as.integer(m),
        bat_uds = NA_real_,
        svf_agree = NA_real_
      )
    } else {
      # udoverlap files written before the rename carry this as `svf_score`.
      # Accept either name so existing outputs do not need a (multi-hour) rerun.
      tibble::tibble(
        model = as.integer(m),
        bat_uds = x$bat_uds,
        svf_agree = if (is.null(x[["svf_agree"]])) {
          x[["svf_score"]]
        } else {
          x[["svf_agree"]]
        }
      )
    }
  })

  ls_df <- tibble::as_tibble(ls) |>
    dplyr::select(model, total_logp, n_steps, delta_logp) |>
    dplyr::mutate(perplexity = exp(-total_logp / n_steps))

  dplyr::full_join(ud_df, ls_df, by = "model") |>
    dplyr::bind_cols(parse_key(k)) |>
    dplyr::mutate(key = k, .before = 1)
})

# Energy scores (one combined file for all deer) ------------------------------
es_df <- readRDS(es_file) |>
  dplyr::mutate(
    key = paste(id, season, year, sep = "_"),
    model = as.integer(model)
  ) |>
  dplyr::select(key, model, energy_score)

# Observed step lengths per deer, from the track ------------------------------
# Model-independent (unlike the amt filter, which recovers sl_ from the clogit
# model frame): the observed steps live on data/tracks/data_<key>.rds$stp. We
# keep the full vector (list-column) so p_excd can be computed per (deer, model)
# against its own energy_score.
sl_df <- purrr::map_dfr(unique(per_deer_df$key), function(k) {
  tfile <- sprintf("data/tracks/data_%s.rds", k)
  if (!file.exists(tfile)) {
    return(tibble::tibble(key = k, sl_obs = list(NA_real_)))
  }
  stp <- readRDS(tfile)$stp[[1]]
  sl_vec <- stp$sl_[!is.na(stp$sl_)]
  if (length(sl_vec) == 0) {
    sl_vec <- NA_real_
  }
  tibble::tibble(key = k, sl_obs = list(sl_vec))
})

# One long data frame ---------------------------------------------------------
# p_excd = P(observed sl > energy_score); >= 0.5 iff energy_score below the
# deer's median step length.
all_df <- per_deer_df |>
  dplyr::left_join(es_df, by = c("key", "model")) |>
  dplyr::left_join(sl_df, by = "key") |>
  dplyr::mutate(
    p_excd = purrr::map2_dbl(sl_obs, energy_score, function(sl, es) {
      if (is.null(sl) || all(is.na(sl)) || is.na(es)) {
        return(NA_real_)
      }
      mean(sl > es, na.rm = TRUE)
    })
  ) |>
  dplyr::select(-sl_obs) |>
  dplyr::arrange(key, model)

# Sequential gates ------------------------------------------------------------
step1 <- all_df |>
  dplyr::filter(!is.na(bat_uds), bat_uds >= ud_min)

step2 <- step1 |>
  dplyr::filter(!is.na(svf_agree), svf_agree >= svf_min)

# Step 3: a model passes iff it beats the null by at least delta_logp_min. This
# is now a plain row-wise filter with no per-deer grouping: the null is no
# longer one of the candidates, so there is no null row to fall back to when a
# deer has no passing model. Such a deer simply drops out here.
step3 <- step2 |>
  dplyr::filter(!is.na(delta_logp), delta_logp >= delta_logp_min)

step4 <- step3 |>
  dplyr::filter(!is.na(p_excd), p_excd >= p_excd_min)

# Annotate each row with the step that eliminated it (NA if it survived all
# four gates). Useful for downstream "what was dropped where" diagnostics.
key_id <- function(d) paste(d$key, d$model)
in_s1 <- key_id(all_df) %in% key_id(step1)
in_s2 <- key_id(all_df) %in% key_id(step2)
in_s3 <- key_id(all_df) %in% key_id(step3)
in_s4 <- key_id(all_df) %in% key_id(step4)
all_df <- all_df |>
  dplyr::mutate(
    dropped_at = dplyr::case_when(
      !in_s1 ~ 1L,
      !in_s2 ~ 2L,
      !in_s3 ~ 3L,
      !in_s4 ~ 4L,
      TRUE ~ NA_integer_
    )
  )

# Save ------------------------------------------------------------------------
saveRDS(all_df, combined_out)
saveRDS(step4, selected_out)

cat(sprintf(
  "\nDeer-model rows: %d -> step1 %d -> step2 %d -> step3 %d -> step4 %d\n",
  nrow(all_df),
  nrow(step1),
  nrow(step2),
  nrow(step3),
  nrow(step4)
))
cat(sprintf(
  "Unique deer:     %d total -> step1 %d -> step2 %d -> step3 %d -> step4 %d\n",
  dplyr::n_distinct(all_df$key),
  dplyr::n_distinct(step1$key),
  dplyr::n_distinct(step2$key),
  dplyr::n_distinct(step3$key),
  dplyr::n_distinct(step4$key)
))

# Plots -----------------------------------------------------------------------
# Two figures of four violin panels (one per metric, one panel per model).
#   *_pre.png  — every (deer, model) row in all_df, before any filtering.
#   *_post.png — each panel shows the gate's input (previous step) split into
#                kept vs. dropped at that gate, so you can see what
#                got eliminated:
#                   bat_uds:    all_df vs. step1
#                   svf_agree:  step1  vs. step2
#                   delta_logp: step2  vs. step3
#                   p_excd:     step3  vs. step4
library(patchwork)

violin_panel <- function(df, y, hline, fill, title = y) {
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = as.factor(model), y = .data[[y]])
  ) +
    ggplot2::geom_violin(fill = fill, alpha = 0.5, trim = TRUE) +
    ggplot2::geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
    ggplot2::geom_hline(
      yintercept = hline,
      colour = "red",
      linetype = 2,
      linewidth = 0.6
    ) +
    ggplot2::labs(x = "Model", y = y, title = title) +
    ggplot2::theme_minimal()
}

# Panel for the post figure: one violin per model over the gate's input set
# (previous step's survivors), with jittered points colored by whether each
# row was kept or dropped at that gate. The threshold line plus the point
# colors are enough to see what was eliminated.
violin_panel_diff <- function(input_df, kept_df, y, hline, fill, title = y) {
  dropped_df <- dplyr::anti_join(
    input_df,
    kept_df,
    by = c("key", "model")
  )
  combined <- dplyr::bind_rows(
    kept_df |> dplyr::mutate(status = "kept"),
    dropped_df |> dplyr::mutate(status = "dropped")
  ) |>
    dplyr::mutate(status = factor(status, levels = c("dropped", "kept")))

  ggplot2::ggplot(
    combined,
    ggplot2::aes(x = as.factor(model), y = .data[[y]])
  ) +
    ggplot2::geom_violin(fill = fill, alpha = 0.4, trim = TRUE) +
    ggplot2::geom_jitter(
      ggplot2::aes(colour = status),
      width = 0.15,
      size = 0.7,
      alpha = 0.7
    ) +
    ggplot2::geom_hline(
      yintercept = hline,
      colour = "red",
      linetype = 2,
      linewidth = 0.6
    ) +
    ggplot2::scale_colour_manual(
      values = c(dropped = "grey50", kept = "#1F77B4"),
      drop = FALSE
    ) +
    ggplot2::labs(x = "Model", y = y, title = title, colour = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

panels_pre <- (violin_panel(all_df, "bat_uds", ud_min, "#FF644E") |
  violin_panel(all_df, "svf_agree", svf_min, "#16E7CF")) /
  (violin_panel(all_df, "delta_logp", delta_logp_min, "#BF5AF2") |
    violin_panel(all_df, "p_excd", p_excd_min, "#61D836")) +
  patchwork::plot_annotation(title = "Pre-filter (all deer x model rows)")

panels_post <- (violin_panel_diff(all_df, step1, "bat_uds", ud_min, "#FF644E") |
  violin_panel_diff(step1, step2, "svf_agree", svf_min, "#16E7CF")) /
  (violin_panel_diff(step2, step3, "delta_logp", delta_logp_min, "#BF5AF2") |
    violin_panel_diff(step3, step4, "p_excd", p_excd_min, "#61D836")) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom")

panels_post <- panels_post +
  patchwork::plot_annotation(
    title = "Per-step filtering (each gate's input split into kept vs. dropped)"
  )

dir.create("plots", showWarnings = FALSE)
ggplot2::ggsave(
  plot_pre,
  panels_pre,
  width = 10,
  height = 7,
  dpi = 300
)
ggplot2::ggsave(
  plot_post,
  panels_post,
  width = 10,
  height = 7,
  dpi = 300
)

# Per-deer pass/fail vs. covariates ------------------------------------------
# A deer "passed" iff at least one environmental-variable model (models 2..4)
# survived all four gates. Deer whose only survivor is the movement-only model
# (1), or with nothing surviving, are counted as failed.
env_models <- c(2L, 3L, 4L)

# The non-environmental candidate. Not the null model — the null never reaches
# this table — but the numbered model with no habitat covariate at all, which
# plays the analogous role *within* the filtered set.
nonenv_models <- 1L

passed_keys <- step4 |>
  dplyr::filter(model %in% env_models) |>
  dplyr::pull(key) |>
  unique()

deer_status <- purrr::map_dfr(keys, parse_key) |>
  dplyr::mutate(
    key = keys,
    id = as.character(id),
    status = factor(
      ifelse(key %in% passed_keys, "passed", "failed"),
      levels = c("failed", "passed")
    )
  )

# Join sex / age covariates from the curated deer track table.
sw_meta <- readRDS("library/SW_filtered_deer.RData") |>
  tibble::as_tibble() |>
  dplyr::transmute(
    id = as.character(id),
    season = as.character(season),
    year = as.integer(year),
    sex = sex,
    age = `age.at.col1`
  ) |>
  dplyr::distinct(id, season, year, .keep_all = TRUE)

deer_status <- deer_status |>
  dplyr::left_join(sw_meta, by = c("id", "season", "year"))

cat(sprintf(
  "\nPer-deer pass/fail (env models %d..%d): passed %d / failed %d (of %d)\n",
  min(env_models),
  max(env_models),
  sum(deer_status$status == "passed"),
  sum(deer_status$status == "failed"),
  nrow(deer_status)
))

# Proportion bar charts: stacked-to-100% so passed/failed shares are
# directly comparable across categories with very different sample sizes.
# Per-category total n is printed above each bar so the reader can still see
# how much data drives each proportion.
# (Sex column stays in deer_status above; just not plotted here.)
status_fill <- c(failed = "grey60", passed = "#1F77B4")

make_status_prop <- function(df, var, x_label, fills = status_fill) {
  d <- df |> dplyr::filter(!is.na(.data[[var]]))

  summ <- d |>
    dplyr::count(.data[[var]], status, name = "n", .drop = FALSE) |>
    dplyr::group_by(.data[[var]]) |>
    dplyr::mutate(prop = n / sum(n), total_n = sum(n)) |>
    dplyr::ungroup()

  n_per <- summ |>
    dplyr::distinct(.data[[var]], total_n)

  ggplot2::ggplot(
    summ,
    ggplot2::aes(x = .data[[var]], y = prop, fill = status)
  ) +
    ggplot2::geom_col() +
    ggplot2::geom_text(
      ggplot2::aes(
        label = ifelse(prop > 0, scales::percent(prop, accuracy = 1), "")
      ),
      position = ggplot2::position_stack(vjust = 0.5),
      colour = "white",
      size = 3
    ) +
    ggplot2::geom_text(
      data = n_per,
      ggplot2::aes(
        x = .data[[var]],
        y = 1.02,
        label = sprintf("n = %d", total_n)
      ),
      inherit.aes = FALSE,
      size = 3,
      vjust = 0
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent,
      limits = c(0, 1.12),
      breaks = seq(0, 1, 0.25),
      expand = ggplot2::expansion(mult = c(0.01, 0))
    ) +
    ggplot2::scale_fill_manual(values = fills, drop = FALSE) +
    ggplot2::labs(x = x_label, y = "Proportion of deer", fill = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

p_season <- make_status_prop(deer_status, "season", "Season")
p_age <- make_status_prop(deer_status, "age", "Age")
# Cast year to factor so it discretizes on the x-axis (otherwise geom_col with
# an integer x renders ungrouped bars on a continuous scale).
p_year <- make_status_prop(
  deer_status |> dplyr::mutate(year = factor(year)),
  "year",
  "Year"
)

panels_cov <- (p_season | p_age | p_year) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom")

panels_cov <- panels_cov +
  patchwork::plot_annotation(
    title = "Per-deer pass/fail by covariate",
    subtitle = sprintf(
      paste("Passed = >=1 surviving model in formulas %d-%d;",
            "failed = none or only movement-only (%s)."),
      min(env_models),
      max(env_models),
      paste(nonenv_models, collapse = ", ")
    )
  )

ggplot2::ggsave(
  sprintf("plots/filter_covariates_gam%s.png", suffix),
  panels_cov,
  width = 13,
  height = 5,
  dpi = 300
)

# Movement-only vs env among deer that did pass the filters -------------------
# "Did pass" here = at least one model survived all four gates (i.e. the deer
# has rows in step4). Deer with nothing in step4 are dropped entirely. Among the
# rest, classify the deer by whether any surviving model is an env model (2-4);
# otherwise its only survivor is the movement-only model (1).
#
# This used to be "null vs env" and compared against models 2/3 while they were
# the null-ish structural models. Those are resource models now and the real
# null sits outside this table entirely, so the contrast is movement-only
# vs env.
deer_compare <- step4 |>
  dplyr::group_by(key) |>
  dplyr::summarize(
    has_env = any(model %in% env_models),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    status = factor(
      ifelse(has_env, "env", "movement"),
      levels = c("movement", "env")
    )
  ) |>
  dplyr::select(key, status) |>
  dplyr::left_join(
    deer_status |> dplyr::select(key, season, year, sex, age),
    by = "key"
  )

cat(sprintf(
  "Movement vs env among passing deer: env %d / movement %d (of %d)\n",
  sum(deer_compare$status == "env"),
  sum(deer_compare$status == "movement"),
  nrow(deer_compare)
))

null_env_fill <- c(movement = "grey60", env = "#1F77B4")

p_season_c <- make_status_prop(
  deer_compare,
  "season",
  "Season",
  fills = null_env_fill
)
p_age_c <- make_status_prop(
  deer_compare,
  "age",
  "Age",
  fills = null_env_fill
)
p_year_c <- make_status_prop(
  deer_compare |> dplyr::mutate(year = factor(year)),
  "year",
  "Year",
  fills = null_env_fill
)

panels_compare <- (p_season_c | p_age_c | p_year_c) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom")

panels_compare <- panels_compare +
  patchwork::plot_annotation(
    title = "Movement-only vs env among deer with surviving models",
    subtitle = sprintf(
      paste("Deer with nothing in step4 excluded.",
            "env = >=1 surviving model in %d-%d;",
            "movement = only model %s."),
      min(env_models),
      max(env_models),
      paste(nonenv_models, collapse = ", ")
    )
  )

ggplot2::ggsave(
  sprintf("plots/filter_movement_vs_env_gam%s.png", suffix),
  panels_compare,
  width = 13,
  height = 5,
  dpi = 300
)
