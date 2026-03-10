#' @description
#' Create binary indicator layers from the categorical wiscland raster.
#' One layer per non-reference category (forest is the reference level).
#' Saves the result as a multi-layer raster for use in iSSF simulations.

library(terra)

# Load categorical raster
env_raster <- rast("env/wiscland/wiscland2.tif")
names(env_raster) <- "wiscland"

# Create binary indicator layers for each non-reference level
# (forest is the reference category, so it is excluded)
agriculture <- ifel(env_raster$wiscland == "agriculture", 1, 0)
names(agriculture) <- "agriculture"

grassland <- ifel(env_raster$wiscland == "grassland", 1, 0)
names(grassland) <- "grassland"

hay <- ifel(env_raster$wiscland == "hay", 1, 0)
names(hay) <- "hay"

oak <- ifel(env_raster$wiscland == "oak", 1, 0)
names(oak) <- "oak"

central.hardwoods <- ifel(env_raster$wiscland == "central.hardwoods", 1, 0)
names(central.hardwoods) <- "central.hardwoods"

other.forest <- ifel(env_raster$wiscland == "other.forest", 1, 0)
names(other.forest) <- "other.forest"

other <- ifel(env_raster$wiscland == "other", 1, 0)
names(other) <- "other"

# Distance to forest edge (central.hardwoods + other.forest)
forest_binary <- ifel(
  env_raster$wiscland == "central.hardwoods" |
    env_raster$wiscland == "other.forest",
  1,
  0
)
dist_to_forest <- terra::distance(ifel(forest_binary == 1, 1, NA), unit = "m")
dist_to_nonforest <- terra::distance(
  ifel(forest_binary == 0, 1, NA),
  unit = "m"
)
dist_forest_edge <- dist_forest_edge <- min(c(
  dist_to_forest,
  dist_to_nonforest
))

names(dist_forest_edge) <- "dist_forest_edge"

# Stack all layers together
env_binary <- c(
  env_raster,
  agriculture,
  grassland,
  hay,
  oak,
  central.hardwoods,
  other.forest,
  other,
  dist_forest_edge
)

# Save
writeRaster(env_binary, "env/wiscland/wiscland2_binary.tif", overwrite = TRUE)

cat("Saved wiscland2_binary.tif with layers:", names(env_binary), "\n")
