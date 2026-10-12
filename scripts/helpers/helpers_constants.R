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

# Landcover classes counted as forest for the forest-edge distance (see
# forest_edge_distance). Every other class, open_water included, is non-forest.
FOREST_CLASSES <- c("forest", "wetland_forested")

# Crop buffer (metres, CRS 6610) added around each deer's track when cropping
# rasters for random-step generation, covariate extraction, HR raster creation,
# simulation, and scoring: room for simulated-path drift before the HR_center /
# HR_edge terms pull the path back. Random-step generation, covariate
# extraction and HR raster creation widen it per deer, to the longest step
# plus CROP_MARGIN_M where that is larger (see deer_crop_buffer), so every
# random endpoint lies inside the window and inside the HR rasters.
CROP_BUFFER_M <- 3000
CROP_MARGIN_M <- 250

# Non-parametric ("nonp") random points, the design the xgb path reads: per
# deer-year, NONP_INNER_SHARE of each step's points are drawn uniformly in
# distance from 0 to the deer's NONP_INNER_QUANTILE step length, the rest
# uniformly in distance from there to its longest step (see
# make_random_pt_extraction).
NONP_INNER_QUANTILE <- 0.9
NONP_INNER_SHARE <- 0.9

# The canonical deer list (see canonical_deer). Years are the deer's `year`
# field. Training and calibration use TRAIN_YEARS, testing uses TEST_YEAR.
TRAIN_YEARS <- 2016:2020
TEST_YEAR <- 2021

# Deer collared at 20 months enter the list only from this age (months, at
# the start of the track; see deer_track_age).
MIN_AGE_MONTHS <- 24

# Deer-seasons (<id>_<season>) left out by hand after inspecting their paths.
# Both the training and the test years go.
EXCLUDED_DEER_SEASONS <- c("7938_nb")

