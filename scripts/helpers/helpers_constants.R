# Shared constants -------------------------------------------------------------
#
# Values used across every stage of the pipeline and by the runner scripts
# directly. They live on their own so that raster, track, amt and GAM helpers
# can all depend on them without depending on each other.
#
# Part of the helper library split out of scripts/helper_functions.R, which
# now sources every file in this folder. Scripts keep sourcing that one
# aggregator, so nothing here needs to be sourced directly.

# Landcover factor levels, reference (intercept) level first. Sourced from the
# annual `library/landcover/landcover_<year>.tif` rasters, whose seasonal bands
# carry 11 categories. `open_water` is intentionally omitted: it is the water /
# step-exclusion mask (see make_water_mask), not a predictor. `forest` is the
# reference level, so it gets no indicator layer and no coefficient.
LANDCOVER_LEVELS <- c(
  "forest",
  "corn",
  "soybeans",
  "alfalfa_hay",
  "small_grains",
  "other_ag",
  "wetland_forested",
  "wetland_open",
  "grassland",
  "developed"
)

# Crop buffer (metres, CRS 6610) added around each deer's track when cropping
# rasters for random-step generation, covariate extraction, HR raster creation,
# simulation, and scoring. Sized to clear the longest single observed/simulated
# step (~2.3 km max) plus room for simulated-path drift before the HR_center /
# HR_edge terms pull the path back. Keep create_hr_rasters' buffer >= this so HR
# rasters fully cover the runtime crop window.
CROP_BUFFER_M <- 3000

