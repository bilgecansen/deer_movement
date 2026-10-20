#' @description
#' Paths of the deer-years residency_tau.R flags as not range resident
#' (initial tau over TAU_MAX_DAYS), longest tau first. Each panel is
#' one deer-year: its fixes joined within bursts and coloured by day since
#' the first fix, on axes in km from the home-range centre the models use
#' (the cell nearest ctmm's mu in data/HR/HRcenter_<key>.tif, marked +).
#'
#' Inputs:  data/residency_tau.rds, library/SW_filtered_deer.RData,
#'          data/HR/HRcenter_<key>.tif
#' Output:  plots/residency_tau_flagged_<page>.png, PER_PAGE panels a page
#'
#' Usage: Rscript scripts/plot_residency_tau.R

# Configuration ---------------------------------------------------------------
PER_PAGE <- 12
N_COL <- 4

# Load packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(terra)
  library(patchwork)
})

# helper functions
source("scripts/helper_functions.R")

tau <- readRDS("data/residency_tau.rds")
flagged <- tau |>
  dplyr::filter(flagged) |>
  dplyr::arrange(dplyr::desc(tau_days))

tracks <- readRDS("library/SW_filtered_deer.RData")
tracks$key <- sprintf("%s_%s_%d", tracks$id, tracks$season,
                      as.integer(tracks$year))
stopifnot(all(flagged$key %in% tracks$key))

# Every fix of a deer-year, the step starts plus each burst's last fix, with
# its burst and day since the first fix
deer_fixes <- function(stp) {
  ends <- stp |>
    dplyr::group_by(burst_) |>
    dplyr::slice_tail(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(burst_, x = x2_, y = y2_, t = t2_)
  dplyr::bind_rows(dplyr::select(stp, burst_, x = x1_, y = y1_, t = t1_),
                   ends) |>
    dplyr::arrange(t) |>
    dplyr::distinct(t, .keep_all = TRUE) |>
    dplyr::mutate(day = as.numeric(t - min(t), units = "days"))
}

max_day <- max(flagged$span_days)
panel <- function(k) {
  row <- flagged[flagged$key == k, ]
  fixes <- deer_fixes(tracks$stp[[which(tracks$key == k)]])
  centre <- terra::rast(sprintf("data/HR/HRcenter_%s.tif", k))
  xy <- terra::xyFromCell(centre, terra::where.min(centre)[1, "cell"])
  fixes$east <- (fixes$x - xy[1, "x"]) / 1000
  fixes$north <- (fixes$y - xy[1, "y"]) / 1000
  # A square window around the fixes and the centre, so every panel has the
  # same shape whatever the path's
  east_range <- range(c(0, fixes$east))
  north_range <- range(c(0, fixes$north))
  half <- 0.52 * max(diff(east_range), diff(north_range))
  lims_east <- mean(east_range) + c(-half, half)
  lims_north <- mean(north_range) + c(-half, half)
  ggplot(fixes, aes(east, north, colour = day)) +
    geom_path(aes(group = burst_), linewidth = 0.3, alpha = 0.7) +
    geom_point(size = 0.5) +
    annotate("point", x = 0, y = 0, shape = 3, size = 4, stroke = 1.2,
             colour = "red") +
    scale_colour_viridis_c("Day", limits = c(0, max_day)) +
    coord_equal(xlim = lims_east, ylim = lims_north, expand = FALSE) +
    labs(x = "East (km)", y = "North (km)",
         title = sprintf("%s (%s)\ntau = %.1f days", k, row$set,
                         row$tau_days)) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 11))
}

pages <- split(flagged$key, ceiling(seq_along(flagged$key) / PER_PAGE))
for (p in seq_along(pages)) {
  g <- patchwork::wrap_plots(lapply(pages[[p]], panel), ncol = N_COL) +
    patchwork::plot_layout(guides = "collect")
  n_row <- ceiling(length(pages[[p]]) / N_COL)
  out <- sprintf("plots/residency_tau_flagged_%d.png", p)
  ggsave(out, g, width = 14, height = 3.4 * n_row + 0.5, dpi = 130)
  cat(sprintf("-> %s\n", out))
}
