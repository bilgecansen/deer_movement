#' @description
#' Fit ctmm to training data for all deer, calculate AKDE, and save a binary
#' home range raster for each deer.
#'
#' Binary raster: inside 95% HR (upper CI bound) = 1, outside = 0.
#' Output: data/HR/HR_<row_no>.tif  (one file per deer)
#'
#' Usage: Rscript create_hr_rasters.R

# Load packages ----------------------------------------------------------------
library(ctmm)
library(terra)
library(tidyverse)
library(sf)
library(furrr)
library(parallel)

# Row range — keep in sync with wrangle_deer_mvt.R so HR file indices match
# the deer row numbers in data_deer_1_119.rds
row_start <- 1
row_end <- 119

# Load data --------------------------------------------------------------------
cat("Loading data...\n")

# Deer movement data
raw_data <- readRDS('Example_code/SW_filtered_deer.RData')

# Prepare sample data
deer_mvt <- raw_data %>%
  dplyr::filter(keep == T & excursion == F & unstable_hr_center == F) %>%
  dplyr::filter(year %in% c(2017, 2018)) %>%
  dplyr::slice(row_start:row_end)

# Pre-extract the FULL step track per deer (training + test). The HR boundary
# is estimated more accurately when all available locations are used; the
# resulting distance-to-edge raster is then read into training data via
# extract_step_variables().
step_list <- lapply(seq_len(nrow(deer_mvt)), function(i) {
  deer_mvt$stp[[i]]
})
rm(raw_data, deer_mvt)
gc()

# Load env_raster — used only as spatial template (CRS + resolution)
env_raster <- terra::rast("env/wiscland/wiscland2_binary.tif")
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
  seq_len(length(step_list)),
  function(i) {
    tryCatch(
      {
        cat(sprintf("Deer %d: starting\n", i))

        # Step 1: Reconstruct full track (training + test) ---------------------
        stp_all <- step_list[[i]]

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
        dist_raster <- terra::distance(env_local, terra::vect(hr_boundary))
        hr_dist <- terra::ifel(hr_binary == 1, dist_raster, -dist_raster)
        names(hr_dist) <- "HR_edge"

        rm(
          hr_high,
          hr_boundary,
          dist_raster,
          env_local,
          crop_extent,
          stp_all
        )
        gc()

        # Save
        out_path <- sprintf("data/HR/HRbin_%d.tif", i)
        terra::writeRaster(hr_binary, out_path, overwrite = TRUE)

        out_path_dist <- sprintf("data/HR/HRedge_%d.tif", i)
        terra::writeRaster(hr_dist, out_path_dist, overwrite = TRUE)

        cat(sprintf(
          "Deer %d: saved to %s and %s\n",
          i,
          out_path,
          out_path_dist
        ))
      },
      error = function(e) {
        warning(sprintf("Deer %d failed: %s", i, conditionMessage(e)))
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
