library(terra)
library(sf)

# the wiscland2 raster saved on my research drive is the template for
# environmental rasters we will use for this project and others in WI it is in
# NAD83(2011) / Wisconsin Transverse Mercator (EPSG:6610) with origin 0,0, a 30
# m resolution (mesh), and an extent (in the given CRS) of xmin = 294840, xmax =
# 770040, ymin = 225120, ymax = 734400.
wiscland_archive_path <- "/Volumes/checastaldo/special/wi_data/env/wiscland/"
wiscland2 <- terra::rast(paste0(wiscland_archive_path, "wiscland2.tif"))
# get the extents and resolution for use later
wl2_ext <- as.vector(terra::ext(wiscland2))
wl2_res <- terra::res(wiscland2)
# here r is the environmental raster you will be resampling onto EPSG:6610, or,
# if it is already in EPSG:6610, this will adjust the origin, mesh, and extent
# to be stackable this process is called warping and you will use gdal directly
# to do this through the sf package for speed and customization

# write the raster to this file. Here we are using datatype = "INT1U" since the
# maximum integer value in this raster example is 255 However this choice is
# really important and depends on the type of numeric data you have. Here are
# some choices
# INT1U  0 to 255            small categorical (cropland data for example)
# INT2U  0 to 65,535         Wiscland2-style 4-digit codes
# INT2S  -32,768 to 32,767   signed integer data with negatives
# INT4U  0 to ~4.3 billion   pixel counts, large IDs
# INT4S  ±~2.1 billion       signed large integers
# FLT4S  32-bit float, range ±3.4e38, ~7 significant digits
# FLT8S  64-bit double, full R precision
# note rasters never contain character data - these are always treated as
# factors with integers and usually involve a RAT table that accompanies the
# raster itself as a file ending tif.aux.xml the RAT table is a data frame with
# unique pixel values and their label (category) The NAflag is the category for
# integers that is NA, as this often is a number such as 0, or 255 etc., best to
# omit this for non-integer data

# make a temporary folder to write the raster input for warping and choose where
# the warped output will live tmp_dir holds the input as an intermediate;
# outfile should be a permanent path you want to keep
tmp_dir <- tempfile()
dir.create(tmp_dir)
infile  <- file.path(tmp_dir, "in.tif")
# CHANGE THIS to where you want the warped raster saved
outfile <- "your/permanent/output/path.tif"

# save the focal raster for warping, r, to the input tif
terra::writeRaster(
  r, infile, overwrite = TRUE,
  wopt = list(datatype = "INT1U", NAflag = 0, gdal = "COMPRESS=LZW")
)

# Warp to target CRS, snapped to Wiscland2's grid using -te (explicit extent
# from Wiscland2) instead of -tap so cells align with the wiscland2 raster for
# integer data we will be using nearest neighbor resampling, "near", with the
# "-r" argument however for floating point data this should be bilinear or cubic
# -ot must match the writeRaster datatype above:
#   INT1U -> "Byte",   INT2U -> "UInt16", INT2S -> "Int16",
#   INT4U -> "UInt32", INT4S -> "Int32",  FLT4S -> "Float32", FLT8S -> "Float64"
# note also the -srcnodata argument needs to match the NAflag argument above
sf::gdal_utils(
  util = "warp",
  source = infile,
  destination = outfile,
  options = c(
    # replace outfile if it already exists
    "-overwrite",
    # target CRS: NAD83(2011) / Wisconsin Transverse Mercator
    "-t_srs", "EPSG:6610",
    # resampling method (near for categorical, bilinear/cubic for continuous)
    "-r", "near",
    # output GDAL datatype (must match writeRaster datatype)
    "-ot", "Byte",
    # target resolution: x and y pixel size, in CRS units (meters here)
    "-tr", as.character(wl2_res[1]), as.character(wl2_res[2]),
    # target extent: xmin ymin xmax ymax, in target CRS units;
    "-te", as.character(wl2_ext["xmin"]), as.character(wl2_ext["ymin"]),
           # together with -tr this forces cell-for-cell alignment with
           # wiscland2
           as.character(wl2_ext["xmax"]), as.character(wl2_ext["ymax"]),
    # value in the input that means "no data" (must match NAflag from
    # writeRaster)
    "-srcnodata", "0",
    # value to use for "no data" in the output
    "-dstnodata", "0",
    # warp option: parallelize across all available CPU cores
    "-wo", "NUM_THREADS=ALL_CPUS",
    # creation option: tiled GeoTIFF (faster random access for large rasters)
    "-co", "TILED=YES",
    # creation option: lossless LZW compression (smaller files, no quality loss)
    "-co", "COMPRESS=LZW",
    # creation option: use BigTIFF format if file might exceed 4 GB
    "-co", "BIGTIFF=IF_SAFER"
))

# If your input raster has a RAT (categorical raster with class labels), carry
# the sidecar over to the output. gdalwarp doesn't propagate the .tif.aux.xml
# file that holds factor labels. Pixel values are unchanged by reprojection, so
# the original RAT remains valid for the warped output.
src_aux <- paste0(infile,  ".aux.xml")
out_aux <- paste0(outfile, ".aux.xml")
if (file.exists(src_aux)) file.copy(src_aux, out_aux, overwrite = TRUE)

# verify the output aligns with wiscland2; if it doesn't stack cleanly,
# downstream analyses will fail
warped <- terra::rast(outfile)
if (terra::compareGeom(wiscland2, warped, stopOnError = FALSE)) {
  message("Output aligned with wiscland2 ✓")
} else {
  warning("Output does NOT match wiscland2 grid. Stacking will fail.")
}
