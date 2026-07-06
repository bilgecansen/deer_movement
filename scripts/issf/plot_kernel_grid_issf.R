#' @description
#' Plot a 3x3 grid of redistribution kernels for one deer + one fitted model,
#' at 9 evenly-spaced observed steps from the deer's track. Each panel shows
#' the kernel raster cropped tight around the observed start (blue point) and
#' observed end (red point).
#'
#' Built with tidyterra::geom_spatraster + patchwork because terra::plot
#' rewrites the device layout per call, which mangles base-R multi-panel
#' grids when writing to a file device.
#'
#' Usage:
#'   Rscript scripts/issf/plot_kernel_grid_issf.R <id> <season> <year> [model_num]
#'     id        — deer ID
#'     season    — season string (e.g. "fa", "nb")
#'     year      — year (integer)
#'     model_num — fitted-model index from results_issf_<key>.rds (default 2)
#'
#' Output: plots/kernel_all_<key>_m<model_num>.png

# Parse CLI args ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3 || length(args) > 4) {
  stop(
    "Usage: Rscript scripts/issf/plot_kernel_grid_issf.R <id> <season> <year> [model_num]"
  )
}

id <- args[1]
season <- args[2]
year <- as.integer(args[3])
model_num <- if (length(args) == 4) as.integer(args[4]) else 2L

key <- sprintf("%s_%s_%d", id, season, year)
cat(sprintf("Deer: %s   Model: %d\n", key, model_num))

# Load packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(amt)
  library(terra)
  library(tidyverse)
  library(sf)
  library(ggplot2)
  library(patchwork)
  library(tidyterra)
})

source("scripts/helper_functions.R")

# Load data --------------------------------------------------------------------
deer_mvt <- readRDS(sprintf("data/tracks/data_%s.rds", key))

# Season-specific annual landcover (categorical band + per-class binary
# indicator layers). env_old terrain covariates are dropped.
env_raster <- load_landcover(year, season)

ndvi_year <- load_ndvi(year)

results_issf <- readRDS(sprintf("results/issf/results_issf_%s.rds", key))

# Crop env around the deer's track + add HR layers ---------------------------
stp_data <- deer_mvt$stp[[1]]

crop_extent <- sf::st_buffer(
  sf::st_as_sf(stp_data, coords = c('x1_', 'y1_'), crs = 6610),
  CROP_BUFFER_M
)

env_cropped <- terra::crop(env_raster, crop_extent)
env_cropped$HR_edge <- load_hr_edge_raster(id, season, year, env_cropped)
env_cropped$HR_center <- load_hr_center_raster(id, season, year, env_cropped)
env_cropped$HR_center_log <- log1p(env_cropped$HR_center)

ndvi_cropped <- terra::crop(ndvi_year, crop_extent)

# Build the issf model object (with coefficient-name renaming) ---------------
# Same coefficient-rename + interaction-permutation trick used by run_sims.R
# and run_logscore.R so the model matches model.matrix columns at runtime.
iss_i <- results_issf[[model_num]]$iss
if (is.character(iss_i)) {
  stop(sprintf("Model %d for %s did not fit successfully.", model_num, key))
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

issf_model <- amt::make_issf_model(
  coefs = coefs,
  sl = iss_i$sl_,
  ta = iss_i$ta_
)

# Pick 9 evenly-spaced steps across the track --------------------------------
step_ids <- round(seq(1, nrow(stp_data), length.out = 9))
cat(sprintf("Step rows: %s\n", paste(step_ids, collapse = ", ")))

# Kernel builder for a single observed step ----------------------------------
build_kernel <- function(i) {
  row_i <- stp_data[i, ]
  mo <- lubridate::month(row_i$t1_)
  env_local <- env_cropped
  env_local$ndvi <- terra::resample(
    ndvi_cropped[[mo]],
    env_local,
    method = "near"
  )

  start_pt <- row_i[, c("x1_", "y1_", "t1_")] |>
    amt::make_track(
      .x = x1_,
      .y = y1_,
      .t = t1_,
      crs = terra::crs(env_local)
    ) |>
    amt::make_start() |>
    amt::mutate(dt = lubridate::hours(4))

  amt::redistribution_kernel(
    x = issf_model,
    map = env_local,
    fun = function(xy, map) {
      xy |>
        amt::extract_covariates(map, where = "both") |>
        amt::time_of_day(include.crepuscule = FALSE, where = "both") |>
        dplyr::mutate(
          tod_start_day = as.integer(tod_start_ == "day"),
          tod_start_night = as.integer(tod_start_ == "night"),
          days = lubridate::yday(t2_) - min(lubridate::yday(t2_)) + 1
        )
    },
    start = start_pt,
    landscape = "discrete",
    as.rast = TRUE
  )
}

# Build the 9 panels ---------------------------------------------------------
cat(sprintf("Building %d kernels...\n", length(step_ids)))

panels <- purrr::imap(step_ids, function(i, panel_idx) {
  row_i <- stp_data[i, ]
  cat(sprintf("  Panel %d (step %d)\n", panel_idx, i))

  k <- tryCatch(build_kernel(i), error = function(e) {
    message(sprintf("    failed: %s", conditionMessage(e)))
    NULL
  })

  if (is.null(k)) {
    return(
      ggplot2::ggplot() +
        ggplot2::labs(title = sprintf("Step %d — failed", i)) +
        ggplot2::theme_void()
    )
  }

  # Crop tight around the two points so the kernel fills the panel.
  pad <- 500 # metres, CRS 6610 units
  ext <- terra::ext(
    min(row_i$x1_, row_i$x2_) - pad,
    max(row_i$x1_, row_i$x2_) + pad,
    min(row_i$y1_, row_i$y2_) - pad,
    max(row_i$y1_, row_i$y2_) + pad
  )
  r <- terra::crop(k$redistribution.kernel, ext)

  ggplot2::ggplot() +
    tidyterra::geom_spatraster(data = r) +
    ggplot2::scale_fill_viridis_c(
      na.value = "transparent",
      guide = "none"
    ) +
    ggplot2::annotate(
      "point",
      x = row_i$x1_,
      y = row_i$y1_,
      colour = "white",
      fill = "blue",
      shape = 21,
      size = 3
    ) +
    ggplot2::annotate(
      "point",
      x = row_i$x2_,
      y = row_i$y2_,
      colour = "white",
      fill = "red",
      shape = 21,
      size = 3
    ) +
    ggplot2::coord_sf(datum = terra::crs(r), expand = FALSE) +
    ggplot2::labs(title = sprintf("Step %d", i)) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 9, hjust = 0.5)
    )
})

combined <- patchwork::wrap_plots(panels, nrow = 3, ncol = 3) +
  patchwork::plot_annotation(
    title = sprintf(
      "Redistribution kernels — deer %s, model %d",
      key,
      model_num
    )
  )

# Save -----------------------------------------------------------------------
dir.create("plots", showWarnings = FALSE)
out_path <- sprintf("plots/kernel_all_%s_m%d.png", key, model_num)

ggplot2::ggsave(
  out_path,
  combined,
  width = 12,
  height = 12,
  dpi = 200
)

cat(sprintf("Saved %s\n", out_path))
