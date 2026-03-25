#' @description
#' Create binary indicator layers from the categorical wiscland raster.
#' One layer per non-reference category (forest is the reference level).
#' Saves the result as a multi-layer raster for use in iSSF simulations.

library(terra)

# Load categorical raster
env_raster <- rast("env/wiscland/wiscland2.tif")

# Reclassify: merge hay -> agriculture, other.forest -> other
levs <- cats(env_raster)[[1]]
rcl <- matrix(c(8, 2, 4, 3), ncol = 2, byrow = TRUE)
env_raster <- classify(env_raster, rcl)
levs <- levs[!levs$label %in% c("hay", "other.forest"), ]
levels(env_raster) <- levs
names(env_raster) <- "wiscland"

# Create binary indicator layers for each non-reference level
# (forest is the reference category, so it is excluded)
agriculture <- ifel(env_raster$wiscland == "agriculture", 1, 0)
names(agriculture) <- "agriculture"

grassland <- ifel(env_raster$wiscland == "grassland", 1, 0)
names(grassland) <- "grassland"

oak <- ifel(env_raster$wiscland == "oak", 1, 0)
names(oak) <- "oak"

central.hardwoods <- ifel(env_raster$wiscland == "central.hardwoods", 1, 0)
names(central.hardwoods) <- "central.hardwoods"

other <- ifel(env_raster$wiscland == "other", 1, 0)
names(other) <- "other"

# Stack all layers together
env_binary <- c(
  env_raster,
  agriculture,
  grassland,
  oak,
  central.hardwoods,
  other
)

# Save
writeRaster(env_binary, "env/wiscland/wiscland2_binary.tif", overwrite = TRUE)

cat("Saved wiscland2_binary.tif with layers:", names(env_binary), "\n")
