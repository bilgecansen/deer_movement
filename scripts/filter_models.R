#' @description
#' Combine per-deer filter outputs (udoverlap, logscore, energy score) plus the
#' observed median step length from the fitted iSSF results, into one long
#' (deer x model) data frame, and apply four sequential model-selection gates.
#'
#' Inputs:
#'   filters/udoverlap_<key>.rds  — list keyed by model: list(bat_uds, svf_score)
#'   filters/logscore_<key>.rds   — df: model, total_logp, n_steps, delta_logp
#'   filters/es.rds               — df: id, season, year, model, energy_score
#'   results/results_issf_<key>.rds — used to recover observed sl_
#'
#' Outputs:
#'   filters/filter_combined.rds  — every (deer, model) with all metrics
#'   filters/filter_selected.rds  — rows surviving all four gates
#'   plots/filter_violins_pre.png, plots/filter_violins_post.png
#'
#' Gates (applied sequentially — only survivors flow into the next step):
#'   1. bat_uds    >= 0.8
#'   2. svf_score  >= 0.8
#'   3. delta_logp >= -3  (delta vs. model 2; null model itself passes only as a
#'                         fallback when no other model in that deer passes)
#'   4. p_excd     >= 0.5  P(observed sl > energy_score); equivalent to ES below
#'                         the deer's median step length.

library(tidyverse)

null_model <- 2L
delta_logp_min <- -3
ud_min <- 0.8
svf_min <- 0.8
p_excd_min <- 0.5

# Discover keys ---------------------------------------------------------------
udov_files <- list.files("filters", "^udoverlap_.*\\.rds$", full.names = TRUE)
logs_files <- list.files("filters", "^logscore_.*\\.rds$", full.names = TRUE)

keys_udov <- gsub("^udoverlap_(.*)\\.rds$", "\\1", basename(udov_files))
keys_logs <- gsub("^logscore_(.*)\\.rds$", "\\1", basename(logs_files))
keys <- intersect(keys_udov, keys_logs)

cat(sprintf(
  "Found %d deer with both udoverlap and logscore outputs (udov: %d, logs: %d)\n",
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
  ud <- readRDS(file.path("filters", sprintf("udoverlap_%s.rds", k)))
  ls <- readRDS(file.path("filters", sprintf("logscore_%s.rds", k)))

  ud_df <- purrr::imap_dfr(ud, function(x, m) {
    if (length(x) == 1 && is.na(x)) {
      tibble::tibble(
        model = as.integer(m),
        bat_uds = NA_real_,
        svf_score = NA_real_
      )
    } else {
      tibble::tibble(
        model = as.integer(m),
        bat_uds = x$bat_uds,
        svf_score = x$svf_score
      )
    }
  })

  ls_df <- tibble::as_tibble(ls) |>
    dplyr::select(model, total_logp, n_steps, delta_logp)

  dplyr::full_join(ud_df, ls_df, by = "model") |>
    dplyr::bind_cols(parse_key(k)) |>
    dplyr::mutate(key = k, .before = 1)
})

# Energy scores (one combined file for all deer) ------------------------------
es_df <- readRDS("filters/es.rds") |>
  dplyr::mutate(
    key = paste(id, season, year, sep = "_"),
    model = as.integer(model)
  ) |>
  dplyr::select(key, model, energy_score)

# Observed step lengths per deer, from the fitted iSSF model frame -----------
# The clogit model frame in results_issf carries `log(sl_)` as a column and the
# response Surv object (which holds case_ as its "status" column). The observed
# steps are case_ == 1. Same vector for every formula on a deer, so we use the
# first non-erroring fit. We keep the full vector (list-column) so that p_excd
# can be computed per (deer, model) against its own energy_score.
sl_df <- purrr::map_dfr(unique(per_deer_df$key), function(k) {
  rfile <- sprintf("results/results_issf_%s.rds", k)
  if (!file.exists(rfile)) {
    return(tibble::tibble(key = k, sl_obs = list(NA_real_)))
  }
  rr <- readRDS(rfile)
  sl_vec <- NA_real_
  for (mm in rr) {
    if (is.character(mm$iss)) {
      next
    }
    frame <- mm$iss$model$model
    if (is.null(frame) || !"log(sl_)" %in% names(frame)) {
      next
    }
    surv_resp <- frame[[1]]
    case_vec <- surv_resp[, "status"]
    log_sl <- frame[["log(sl_)"]]
    sl_vec <- exp(log_sl[case_vec == 1])
    break
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
  dplyr::filter(!is.na(svf_score), svf_score >= svf_min)

# Step 3: every non-null model passes iff its delta_logp (vs. model 2) >=
# threshold. The null model (id = 2) is kept only as a fallback for deer where
# no non-null model passed. Operates on step2 survivors, so a null model that
# itself failed steps 1 or 2 is not available as a fallback.
step3 <- step2 |>
  dplyr::group_by(key) |>
  dplyr::group_modify(function(g, .y) {
    passing_alt <- g |>
      dplyr::filter(
        model != null_model,
        !is.na(delta_logp),
        delta_logp >= delta_logp_min
      )
    if (nrow(passing_alt) > 0) {
      passing_alt
    } else {
      g |> dplyr::filter(model == null_model)
    }
  }) |>
  dplyr::ungroup()

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
      TRUE   ~ NA_integer_
    )
  )

# Save ------------------------------------------------------------------------
saveRDS(all_df, "filters/filter_combined.rds")
saveRDS(step4, "filters/filter_selected.rds")

cat(sprintf(
  "\nDeer-model rows: %d total -> step1 %d -> step2 %d -> step3 %d -> step4 %d\n",
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
#                kept vs. dropped at that gate, so you can see what got
#                eliminated:
#                   bat_uds:    all_df vs. step1
#                   svf_score:  step1  vs. step2
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
  violin_panel(all_df, "svf_score", svf_min, "#16E7CF")) /
  (violin_panel(all_df, "delta_logp", delta_logp_min, "#BF5AF2") |
    violin_panel(all_df, "p_excd", p_excd_min, "#61D836")) +
  patchwork::plot_annotation(title = "Pre-filter (all deer x model rows)")

panels_post <- (
  violin_panel_diff(all_df, step1, "bat_uds",    ud_min,         "#FF644E") |
  violin_panel_diff(step1,  step2, "svf_score",  svf_min,        "#16E7CF")
) / (
  violin_panel_diff(step2,  step3, "delta_logp", delta_logp_min, "#BF5AF2") |
  violin_panel_diff(step3,  step4, "p_excd",     p_excd_min,     "#61D836")
) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom")

panels_post <- panels_post +
  patchwork::plot_annotation(
    title = "Per-step filtering (each gate's input split into kept vs. dropped)"
  )

dir.create("plots", showWarnings = FALSE)
ggplot2::ggsave(
  "plots/filter_violins_pre.png",
  panels_pre,
  width = 10,
  height = 7,
  dpi = 300
)
ggplot2::ggsave(
  "plots/filter_violins_post.png",
  panels_post,
  width = 10,
  height = 7,
  dpi = 300
)
