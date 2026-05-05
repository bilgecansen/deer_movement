#' @description
#' Fit ctmm to all training+test locations for each (deer_id, season, year)
#' combination, calculate AKDE, and save three per-deer rasters:
#'   * HRbin_<id>_<season>_<year>.tif    — binary 95% HR (upper CI)
#'   * HRedge_<id>_<season>_<year>.tif   — signed distance to HR edge (+inside)
#'   * HRcenter_<id>_<season>_<year>.tif — distance to ctmm μ (HR center)
#'
#' Usage: Rscript create_hr_rasters.R

# Load packages ----------------------------------------------------------------
library(ctmm)
library(terra)
library(tidyverse)
library(sf)
library(furrr)
library(parallel)

# Load data --------------------------------------------------------------------
cat("Loading data...\n")

# Deer movement data
sw_deer_tracks <- readRDS('data/SW_filtered_deer.RData')

# Filter to the rows we want to model
deer_mvt <- sw_deer_tracks %>%
  dplyr::filter(keep == T & excursion == F & unstable_hr_center == F) %>%
  dplyr::filter(year %in% 2017:2021)

# Pre-extract per-deer inputs as a flat list — avoids serializing the full
# nested deer_mvt object (with all nested columns) into every worker process.
# Each element carries the metadata needed for output filenames plus the FULL
# track (train + test) used to estimate the HR boundary.
step_list <- lapply(seq_len(nrow(deer_mvt)), function(i) {
  list(
    id = deer_mvt$id[i],
    season = deer_mvt$season[i],
    year = deer_mvt$year[i],
    stp = deer_mvt$stp[[i]]
  )
})
rm(sw_deer_tracks, deer_mvt)
gc()

# Load env_raster — used only as spatial template (CRS + resolution)
env_raster <- terra::rast("data/env/wiscland/wiscland2_binary.tif")
env_old <- terra::rast("Example_code/Env_2017.tif")
env_raster <- terra::crop(env_raster, env_old) %>%
  terra::resample(env_old)

# Wrap a single layer for sending to parallel workers
env_template <- terra::wrap(env_raster[[1]])
rm(env_raster, env_old)
gc()

# Create output directory
dir.create("data/HR", recursive = TRUE, showWarnings = FALSE)

# Process all deer in parallel -------------------------------------------------
cat(sprintf("Processing %d deer on 5 workers...\n", length(step_list)))

future::plan(multisession, workers = 5)

furrr::future_map(
  step_list,
  function(deer_input) {
    id <- deer_input$id
    season <- deer_input$season
    year <- deer_input$year
    deer_key <- sprintf("%s_%s_%d", id, season, year)

    tryCatch(
      {
        cat(sprintf("Deer %s: starting\n", deer_key))

        # Step 1: Reconstruct full track (training + test) ---------------------
        stp_all <- deer_input$stp

        # All step start locations (every GPS fix except the last per burst)
        locs_start <- stp_all %>%
          dplyr::select(x = x1_, y = y1_, timestamp = t1_)

        # Last endpoint of each burst (the final GPS fix)
        locs_end <- stp_all %>%
          dplyr::group_by(burst_) %>%
          dplyr::slice_tail(n = 1) %>%
          dplyr::ungroup() %>%
          dplyr::select(x = x2_, y = y2_, timestamp = t2_)

        track_df <- dplyr::bind_rows(locs_start, locs_end) %>%
          dplyr::arrange(timestamp) %>%
          dplyr::distinct(timestamp, .keep_all = TRUE)

        rm(locs_start, locs_end)

        # Step 2: Convert to ctmm telemetry ------------------------------------
        # ctmm expects longitude/latitude — same pattern as overlap_ud()
        tele <- track_df %>%
          sf::st_as_sf(coords = c("x", "y"), crs = 6610) %>%
          sf::st_transform(4326) %>%
          dplyr::mutate(
            longitude = sf::st_coordinates(geometry)[, 1],
            latitude = sf::st_coordinates(geometry)[, 2]
          ) %>%
          sf::st_drop_geometry() %>%
          as.data.frame() %>%
          ctmm::as.telemetry()

        rm(track_df)
        gc()

        # Step 3: Fit ctmm model -----------------------------------------------
        invisible(capture.output(
          fit <- ctmm::ctmm.select(
            tele,
            ctmm::ctmm.guess(tele, interactive = FALSE),
            verbose = FALSE,
            cores = 1
          )
        ))

        # Step 3b: Capture HR center (ctmm μ) in CRS 6610 ----------------------
        # μ is in ctmm's internal AEQD projection (set by `tele`); reproject it
        # to the working CRS so the runtime distance raster lives on the same
        # grid as HRbin/HRedge.
        mu_proj <- ctmm::projection(tele)
        mu_pt <- sf::st_sfc(
          sf::st_point(as.numeric(fit$mu)),
          crs = mu_proj
        ) %>%
          sf::st_transform(6610)

        # Step 4: Calculate AKDE -----------------------------------------------
        invisible(capture.output(
          ud <- ctmm::akde(tele, fit)
        ))

        rm(tele, fit)
        gc()

        # Step 5: Extract 95% HR polygon at upper CI bound --------------------
        # Returns SpatialPolygonsDataFrame with 3 rows: low, est, high CI
        hr_spdf <- ctmm:::SpatialPolygonsDataFrame.UD(
          ud,
          level.UD = 0.95,
          level = 0.95
        )
        hr_high <- sf::st_as_sf(hr_spdf[3, ]) %>% # row 3 = high CI
          sf::st_transform(6610)

        rm(ud, hr_spdf)
        gc()

        # Step 6: Create binary raster -----------------------------------------
        # Crop extent: 5000m buffer around all step start locations (full
        # track), so the resulting HR rasters cover both training and test
        # buffers used downstream.
        crop_extent <- sf::st_buffer(
          sf::st_as_sf(stp_all, coords = c("x1_", "y1_"), crs = 6610),
          5000
        )

        env_local <- terra::unwrap(env_template) %>%
          terra::crop(crop_extent)

        hr_binary <- terra::rasterize(
          terra::vect(hr_high),
          env_local,
          field = 1,
          background = 0
        )

        # Step 7: Signed distance to HR edge -----------------------------------
        # Positive inside HR, negative outside — gradient always points inward,
        # so the deer is pulled back toward the HR even if it slips outside.
        hr_boundary <- sf::st_boundary(hr_high)
        edge_dist <- terra::distance(env_local, terra::vect(hr_boundary))
        hr_edge <- terra::ifel(hr_binary == 1, edge_dist, -edge_dist)
        names(hr_edge) <- "HR_edge"

        # Step 8: Distance to HR center (ctmm μ) -------------------------------
        # Positive everywhere; the further the cell from the home-range center,
        # the larger the value. Lives on the same grid as HRbin/HRedge.
        hr_center <- terra::distance(env_local, terra::vect(mu_pt))
        names(hr_center) <- "HR_center"

        rm(
          hr_high,
          hr_boundary,
          edge_dist,
          mu_pt,
          env_local,
          crop_extent,
          stp_all
        )
        gc()

        # Save
        out_bin <- sprintf("data/HR/HRbin_%s.tif", deer_key)
        out_edge <- sprintf("data/HR/HRedge_%s.tif", deer_key)
        out_center <- sprintf("data/HR/HRcenter_%s.tif", deer_key)

        terra::writeRaster(hr_binary, out_bin, overwrite = TRUE)
        terra::writeRaster(hr_edge, out_edge, overwrite = TRUE)
        terra::writeRaster(hr_center, out_center, overwrite = TRUE)

        cat(sprintf(
          "Deer %s: saved %s, %s, %s\n",
          deer_key,
          out_bin,
          out_edge,
          out_center
        ))
      },
      error = function(e) {
        warning(sprintf(
          "Deer %s failed: %s",
          deer_key,
          conditionMessage(e)
        ))
      }
    )
  },
  .options = furrr::furrr_options(
    packages = c("ctmm", "terra", "sf", "dplyr"),
    stdout = TRUE,
    seed = TRUE
  )
)

future::plan(sequential)
cat("Done.\n")
