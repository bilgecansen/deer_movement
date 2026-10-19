#' @description
#' Path-level used-habitat calibration (UHC; Fieberg et al. 2018) for NDVI:
#' do the simulated paths spend their time at the same NDVI as the real
#' deer? The NDVI counterpart of uhc_landcover_xgb.R, on the same paths:
#'   the season's two xgb deer picks (sim_paths_xgb.R),
#'   the GAM null and the Gilbertson et al. (2025) model
#'   (sim_paths_gam_xgb.R).
#' Every model's paths start at the deer's first fix and run free on the
#' same clock, so the first fix is left out of both. Each simulated path is
#' read at the clock tick of every later real fix (xgb_path_clock() puts a
#' tick within TOLERANCE_MINS of each), and every value, real or simulated,
#' from the monthly NDVI layer of that real fix's time
#' (xgb_path_ndvi_layer(), as the kernels read it); a position with no NDVI
#' layer in reach, or off the window, has no value. nb has no NDVI and is
#' left out.
#'
#' Per deer-year and model, the Wasserstein-1 distance (wasserstein1())
#' between the NDVI at the real fixes and the NDVI over the model's 10 paths
#' pooled: how far, in NDVI units, the simulated values sit from the real
#' ones on average. Each model is set against the GAM null deer-year by
#' deer-year.
#'
#' Inputs:  data/tracks/data_<key>.rds, sims_xgb_<key>.rds and
#'          sims_gam_<key>.rds in xgb_sims_dir(SET, REACH), the rasters the
#'          wrangle read
#' Output:  filters/xgb/uhc_ndvi_<tag>.rds — the NDVI values per deer-year,
#'          source (real or model) and path; and
#'          plots/uhc_ndvi_w1_<tag>.png; <tag> is xgb_scores_tag(SET, REACH)
#'
#' Configuration: edit the block below before running. The set can also be
#' given on the command line: Rscript scripts/xgb/uhc_ndvi_xgb.R train

# Configuration ---------------------------------------------------------------
SET <- "test"
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1) {
  SET <- .args[1]
}
stopifnot(SET %in% c("test", "train"))
REACH <- "gamma99"
TOLERANCE_MINS <- 30
# FALSE: reuse the saved values if they exist, and only summarise
RECOUNT <- FALSE
N_PAR <- parallel::detectCores() - 1L

# Load packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(furrr)
})

# helper functions
source("scripts/helper_functions.R")

tag <- xgb_scores_tag(SET, REACH)
sims_dir <- xgb_sims_dir(SET, REACH)
deer <- canonical_deer(readRDS("library/SW_filtered_deer.RData"))
deer <- deer[deer$set == SET & deer$season != "nb", ]
keys <- sprintf("%s_%s_%d", deer$id, deer$season, as.integer(deer$year))
have <- file.exists(file.path(sims_dir, sprintf("sims_xgb_%s.rds", keys))) &
  file.exists(file.path(sims_dir, sprintf("sims_gam_%s.rds", keys)))
if (any(!have)) {
  cat(sprintf("%d deer-years lack xgb or GAM paths and are left out\n",
              sum(!have)))
}
keys <- keys[have]

values_one <- function(key) {
  source("scripts/helper_functions.R")
  track <- readRDS(sprintf("data/tracks/data_%s.rds", key))
  stp <- track$stp[[1]]
  xs <- readRDS(file.path(sims_dir, sprintf("sims_xgb_%s.rds", key)))
  gs <- readRDS(file.path(sims_dir, sprintf("sims_gam_%s.rds", key)))
  stopifnot(identical(xs$times, gs$times))
  times <- xs$times

  # The real fixes after the first, and the clock tick of each
  last_of_burst <- !duplicated(stp$burst_, fromLast = TRUE)
  real <- dplyr::bind_rows(
    tibble::tibble(t_ = stp$t1_, x_ = stp$x1_, y_ = stp$y1_),
    tibble::tibble(t_ = stp$t2_[last_of_burst], x_ = stp$x2_[last_of_burst],
                   y_ = stp$y2_[last_of_burst])
  ) |>
    dplyr::distinct(t_, .keep_all = TRUE) |>
    dplyr::arrange(t_)
  real <- real[-1, ]
  gap <- abs(outer(as.numeric(real$t_), as.numeric(times), "-"))
  tick <- max.col(-gap, ties.method = "first")
  stopifnot(gap[cbind(seq_along(tick), tick)] <= TOLERANCE_MINS * 60)

  land <- xgb_path_landscape(track, load_landcover(track$year, track$season),
                             load_ndvi(track$year), load_topo())
  layer <- xgb_path_ndvi_layer(real$t_, land$ndvi_times)
  ndvi_at <- function(x, y) {
    land$ndvi[cbind(xgb_path_cell(land, x, y), layer)]
  }

  picks <- xs$picks
  pick_label <- ifelse(picks$best, "xgb best", "xgb simplest")
  names(pick_label) <- picks$config
  model_paths <- c(
    stats::setNames(xs$paths, pick_label[names(xs$paths)]),
    list(`GAM null` = gs$paths$null_gam,
         `Gilbertson` = gs$paths$gilbertson_gam)
  )
  sim <- purrr::imap_dfr(model_paths, function(p, model) {
    purrr::map_dfr(sort(unique(p$nsim)), function(r) {
      pr <- p[p$nsim == r, ]
      stopifnot(identical(pr$t_, times))
      pos <- pr[tick, ]
      tibble::tibble(source = model, path = r, ndvi = ndvi_at(pos$x_, pos$y_))
    })
  })
  dplyr::bind_rows(
    tibble::tibble(source = "real", path = NA_integer_,
                   ndvi = ndvi_at(real$x_, real$y_)),
    sim
  ) |>
    dplyr::filter(!is.na(ndvi)) |>
    dplyr::mutate(key = key, season = track$season)
}

out_path <- sprintf("filters/xgb/uhc_ndvi_%s.rds", tag)
if (RECOUNT || !file.exists(out_path)) {
  cat(sprintf("reading NDVI for %d %s deer-years (%s)\n", length(keys), tag,
              format(Sys.time(), "%H:%M")))
  future::plan(future::multisession, workers = min(N_PAR, length(keys)))
  vals <- furrr::future_map_dfr(keys, values_one,
                                .options = furrr::furrr_options(
                                  seed = TRUE, chunk_size = 1L
                                ))
  future::plan(future::sequential)
  saveRDS(vals, out_path)
  cat(sprintf("-> %s (%s)\n", out_path, format(Sys.time(), "%H:%M")))
} else {
  vals <- readRDS(out_path)
  stopifnot(setequal(unique(vals$key), keys))
  cat(sprintf("reusing %s\n", out_path))
}

# Wasserstein-1 per deer-year and model --------------------------------------
models <- c("xgb best", "xgb simplest", "GAM null", "Gilbertson")
w1 <- vals |>
  dplyr::filter(source != "real") |>
  dplyr::group_by(key, season, source) |>
  dplyr::summarise(sim = list(ndvi), .groups = "drop") |>
  dplyr::left_join(vals |>
                     dplyr::filter(source == "real") |>
                     dplyr::group_by(key) |>
                     dplyr::summarise(real = list(ndvi), .groups = "drop"),
                   by = "key") |>
  dplyr::mutate(w1 = purrr::map2_dbl(real, sim, wasserstein1),
                source = factor(source, models)) |>
  dplyr::select(key, season, source, w1)
cat("\n=== Wasserstein distance between real and simulated NDVI, median",
    "over deer-years (NDVI units) ===\n")
print(as.data.frame(w1 |>
                      dplyr::group_by(season, source) |>
                      dplyr::summarise(w1 = median(w1), .groups = "drop") |>
                      tidyr::pivot_wider(names_from = source,
                                         values_from = w1)),
      digits = 3, row.names = FALSE)
w1_wide <- w1 |> tidyr::pivot_wider(names_from = source, values_from = w1)
cat("\n=== deer-years where a model comes closer than the GAM null ===\n")
for (m in c("xgb best", "xgb simplest", "Gilbertson")) {
  d <- w1_wide[[m]] - w1_wide[["GAM null"]]
  cat(sprintf("%-13s closer in %d of %d (median difference %.4f)\n", m,
              sum(d < 0), length(d), median(d)))
}

# Plot -----------------------------------------------------------------------
cols <- c("xgb best" = "#1B9E77", "xgb simplest" = "#66A61E",
          "Gilbertson" = "#D95F02")
season_lab <- c(br = "Breeding", fa = "Fawning", pf = "Post-fawning")
d_w1 <- w1_wide |>
  tidyr::pivot_longer(c("xgb best", "xgb simplest", "Gilbertson"),
                      names_to = "model", values_to = "w1") |>
  dplyr::mutate(diff = w1 - .data[["GAM null"]],
                model = factor(model, names(cols)))
n_closer <- d_w1 |>
  dplyr::group_by(season, model) |>
  dplyr::summarise(label = sprintf("%d of %d", sum(diff < 0), dplyr::n()),
                   .groups = "drop")
ylim_w1 <- stats::quantile(d_w1$diff, c(0.01, 0.99))
g <- ggplot(d_w1, aes(model, diff, fill = model)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey30") +
  geom_boxplot(width = 0.6, linewidth = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.15, height = 0, size = 0.8, alpha = 0.45) +
  geom_text(data = n_closer, aes(model, ylim_w1[1], label = label),
            inherit.aes = FALSE, vjust = 0, size = 3.3) +
  scale_fill_manual(NULL, values = cols, guide = "none") +
  facet_wrap(~season, ncol = 3, labeller = labeller(season = season_lab)) +
  coord_cartesian(ylim = ylim_w1) +
  labs(x = NULL,
       y = "Wasserstein distance minus the GAM null's, per deer-year",
       title = sprintf(paste("NDVI per deer-year: distance from the real",
                             "deer against the GAM null, %s"), SET),
       subtitle = paste("Below 0: the model's paths sit at NDVI closer to",
                        "the real deer's than the null's do. Numbers:",
                        "deer-years below 0."),
       caption = paste("Wasserstein-1 distance (NDVI units) between NDVI at",
                       "the real fixes and at the model's 10 paths pooled;",
                       "the 1% most extreme differences fall outside the",
                       "axis.")) +
  theme_bw(base_size = 11)
ggsave(sprintf("plots/uhc_ndvi_w1_%s.png", tag), g, width = 12, height = 6,
       dpi = 140)
cat(sprintf("-> plots/uhc_ndvi_w1_%s.png\n", tag))
