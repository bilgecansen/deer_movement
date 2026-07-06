#' @description
#' Per-deer plots of observed track vs. simulated paths, faceted by model.
#' Reads `filters/issf/filter_combined_issf.rds` (which carries a `dropped_at` column
#' marking the gate that eliminated each row; NA = survived all four).
#'
#' Outputs five PDFs in plots/, all with the same aesthetic:
#'   deer_paths_selected.pdf            — models that passed every gate
#'   deer_paths_step1_eliminated.pdf    — failed bat_uds    gate
#'   deer_paths_step2_eliminated.pdf    — failed svf_score  gate
#'   deer_paths_step3_eliminated.pdf    — failed delta_logp gate (incl. null
#'                                         model when an alt passed)
#'   deer_paths_step4_eliminated.pdf    — failed p_excd     gate

library(tidyverse)

annotated <- readRDS("filters/issf/filter_combined_issf_test.rds")

dir.create("plots", showWarnings = FALSE)

# Plot one deer: observed path (orange) + simulated paths from the chosen
# models (black, low alpha), faceted by model. Returns NULL if track or sim
# file is missing, or if no requested model has usable simulations.
plot_deer_models <- function(key, model_info, title, sim_alpha = 0.25) {
  track_file <- sprintf("data/tracks/data_%s.rds", key)
  sim_file <- sprintf("sims/issf/sims_issf_%s.rds", key)

  if (!file.exists(track_file) || !file.exists(sim_file)) {
    return(NULL)
  }

  one_deer <- readRDS(track_file)
  obs <- one_deer$stp[[1]] |>
    dplyr::select(x_ = x1_, y_ = y1_, t_ = t1_, burst_)

  results_sim <- readRDS(sim_file)
  keep_models <- as.character(model_info$model)

  sim_df <- purrr::imap_dfr(results_sim, function(sim_m, m) {
    if (!m %in% keep_models) {
      return(NULL)
    }
    if (is.null(sim_m) || (length(sim_m) == 1 && is.na(sim_m))) {
      return(NULL)
    }
    sim_m |>
      dplyr::as_tibble() |>
      dplyr::mutate(model = as.integer(m))
  })

  if (nrow(sim_df) == 0) {
    return(NULL)
  }

  sim_df <- sim_df |>
    dplyr::mutate(path_id = paste(model, nsim, sep = "_"))

  # Facet strip: model number + the four filter metrics so the reader can see
  # which one failed. NA prints as "NA" — no special casing needed.
  fmt <- function(x, d) ifelse(is.na(x), "NA", sprintf(sprintf("%%.%df", d), x))
  strip_lbl <- stats::setNames(
    sprintf(
      "model %d\nUD %s | SVF %s | dlogp %s | P(sl>es) %s",
      model_info$model,
      fmt(model_info$bat_uds, 2),
      fmt(model_info$svf_score, 2),
      fmt(model_info$delta_logp, 1),
      fmt(model_info$p_excd, 2)
    ),
    as.character(model_info$model)
  )

  ggplot() +
    geom_path(
      data = obs,
      aes(x = x_, y = y_, group = burst_),
      colour = "orange",
      linewidth = 0.7
    ) +
    geom_path(
      data = sim_df,
      aes(x = x_, y = y_, group = path_id),
      colour = "black",
      alpha = sim_alpha,
      linewidth = 0.3
    ) +
    coord_equal() +
    facet_wrap(~model, labeller = ggplot2::as_labeller(strip_lbl)) +
    labs(title = title) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      strip.text = element_text(size = 8)
    )
}

# Write one PDF, one page per deer, faceted across the rows of `subset_df`.
write_pdf <- function(subset_df, outfile, title_prefix) {
  keys <- sort(unique(subset_df$key))
  cat(sprintf(
    "%s: %d deer, %d (deer, model) pairs -> %s\n",
    title_prefix,
    length(keys),
    nrow(subset_df),
    outfile
  ))

  if (length(keys) == 0) {
    cat("  [skip] nothing to plot\n")
    return(invisible(NULL))
  }

  pdf(outfile, width = 10, height = 8, onefile = TRUE)
  for (k in keys) {
    model_info <- subset_df |>
      dplyr::filter(key == k) |>
      dplyr::arrange(model)
    p <- plot_deer_models(
      k,
      model_info,
      title = sprintf(
        "%s — %s (%d model%s)",
        title_prefix,
        k,
        nrow(model_info),
        if (nrow(model_info) == 1) "" else "s"
      )
    )
    if (!is.null(p)) {
      print(p)
    } else {
      cat(sprintf("  [skip] %s: no usable sims for requested models\n", k))
    }
  }
  dev.off()
}

# Five outputs ----------------------------------------------------------------
write_pdf(
  annotated |> dplyr::filter(is.na(dropped_at)),
  "plots/deer_paths_selected.pdf",
  "Selected"
)

for (s in 1:4) {
  write_pdf(
    annotated |> dplyr::filter(dropped_at == s),
    sprintf("plots/deer_paths_step%d_eliminated.pdf", s),
    sprintf("Step %d eliminated", s)
  )
}
