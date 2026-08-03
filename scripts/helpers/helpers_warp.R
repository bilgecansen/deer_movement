# Raster warping utility -------------------------------------------------------
#
# warp_to_template(): a standalone gdalwarp wrapper with datatype/range
# checking, used by the one-off raster preparation scripts (wrangle_wiscland.R,
# assign_landcover.R, combine_ndvi.R). Nothing in the modelling pipeline calls
# it; it is kept apart so it does not sit in the load path of every run.
#
# Part of the helper library split out of scripts/helper_functions.R, which
# now sources every file in this folder. Scripts keep sourcing that one
# aggregator, so nothing here needs to be sourced directly.

#' Warp a raster to match a template's CRS / extent / resolution
#'
#' Wraps gdalwarp via sf::gdal_utils. Same job as terra::project()/resample(),
#' but multi-threaded, streaming-capable, and with full control over output
#' datatype, nodata, compression, and tiling.
#'
#' Default datatype is INT2U — covers up to 65535, fits Wiscland-style
#' 4-digit codes and most WI categorical rasters in the smallest size.
#' Override with FLT4S for continuous data (NDVI, elevation, etc.).
#'
#' Before warping, the function checks the input's data range against the
#' chosen datatype's bounds; out-of-range values raise an error (rather than
#' being silently clipped by gdalwarp). It also warns when float-storage
#' input is being converted to an integer datatype, since fractional values
#' would be truncated.
#'
#' @param r           SpatRaster or path to input raster
#' @param template    SpatRaster or path defining target CRS / extent / resolution
#' @param outfile     Output file path (string)
#' @param method      GDAL resampling method ("near", "bilinear", "cubic", ...).
#'                    Default "near".
#' @param datatype    Terra-style datatype string. Options:
#'                    INT1U, INT2U, INT2S, INT4U, INT4S, FLT4S, FLT8S.
#'                    Default "INT2U".
#' @param nodata      Value treated as NA in input AND output. NA disables the
#'                    -srcnodata/-dstnodata flags. Default NA.
#' @param compress    GDAL compression scheme. Default "LZW".
#' @param threads     Worker threads for warping. Default "ALL_CPUS".
#' @param verify      If TRUE, run terra::compareGeom() against template after
#'                    warp and warn on misalignment. Default TRUE.
#' @param check_range If TRUE, scan input min/max and verify they fit the
#'                    target datatype's bounds before warping. Adds a one-time
#'                    minmax compute cost; turn off for very large rasters
#'                    you've already validated. Default TRUE.
#' @return SpatRaster wrapping the warped output file (lazy reference).
warp_to_template <- function(
  r,
  template,
  outfile,
  method = "near",
  datatype = "INT2U",
  nodata = NA,
  compress = "LZW",
  threads = "ALL_CPUS",
  verify = TRUE,
  check_range = TRUE
) {
  # Terra datatype -> GDAL -ot. The two must match or downstream stacks
  # will silently disagree on dtype.
  gdal_dt <- c(
    INT1U = "Byte",
    INT2U = "UInt16",
    INT2S = "Int16",
    INT4U = "UInt32",
    INT4S = "Int32",
    FLT4S = "Float32",
    FLT8S = "Float64"
  )
  if (!datatype %in% names(gdal_dt)) {
    stop(
      "datatype must be one of: ",
      paste(names(gdal_dt), collapse = ", ")
    )
  }
  ot <- gdal_dt[[datatype]]

  # Resolve input to a SpatRaster handle (for the range / type check below).
  if (inherits(r, "SpatRaster")) {
    r_in <- r
  } else if (is.character(r)) {
    r_in <- terra::rast(r)
  } else {
    stop("r must be a SpatRaster or a file path string")
  }

  # Datatype compatibility check — catches silent clipping that would
  # otherwise happen during gdalwarp / writeRaster.
  if (check_range) {
    dt_bounds <- list(
      INT1U = c(0, 255),
      INT2U = c(0, 65535),
      INT2S = c(-32768, 32767),
      INT4U = c(0, 4294967295),
      INT4S = c(-2147483648, 2147483647),
      FLT4S = c(-3.4e38, 3.4e38),
      FLT8S = c(-Inf, Inf)
    )
    bounds <- dt_bounds[[datatype]]

    mm <- terra::minmax(r_in, compute = TRUE)
    data_min <- suppressWarnings(min(mm[1, ], na.rm = TRUE))
    data_max <- suppressWarnings(max(mm[2, ], na.rm = TRUE))

    if (is.finite(data_min) && data_min < bounds[1]) {
      stop(sprintf(
        "Input min %g is below datatype %s lower bound %g — 
        would silently clip. Pick a wider datatype.",
        data_min,
        datatype,
        bounds[1]
      ))
    }
    if (is.finite(data_max) && data_max > bounds[2]) {
      stop(sprintf(
        "Input max %g is above datatype %s upper bound %g — 
        would silently clip. Pick a wider datatype.",
        data_max,
        datatype,
        bounds[2]
      ))
    }

    # Float-storage -> integer datatype: would truncate fractional values.
    # Storage check first (cheap); if storage is float, also confirm the
    # data actually has fractional values via the min/max we already have.
    if (datatype %in% c("INT1U", "INT2U", "INT2S", "INT4U", "INT4S")) {
      in_dtype <- terra::datatype(r_in)
      if (any(grepl("FLT", in_dtype))) {
        has_fractional <- any(mm != round(mm), na.rm = TRUE)
        if (has_fractional) {
          warning(sprintf(
            "Input has float storage (%s) with fractional values; 
            target datatype %s is integer — values will be truncated.",
            paste(unique(in_dtype), collapse = ","),
            datatype
          ))
        }
      }
    }
  }

  # Resolve template (path or SpatRaster) and pull its grid spec.
  if (is.character(template)) {
    template <- terra::rast(template)
  }
  tmpl_ext <- as.vector(terra::ext(template))
  tmpl_res <- terra::res(template)

  # Use EPSG:CODE form when available (compact, unambiguous); fall back to
  # full WKT otherwise. -t_srs accepts both.
  crs_info <- terra::crs(template, describe = TRUE)
  tmpl_crs_str <- if (!is.na(crs_info$code)) {
    paste0(crs_info$authority, ":", crs_info$code)
  } else {
    terra::crs(template)
  }

  # gdalwarp reads from disk, so a SpatRaster argument gets written to a temp
  # file first. Path inputs go through directly.
  cleanup_dir <- NULL
  if (inherits(r, "SpatRaster")) {
    cleanup_dir <- tempfile()
    dir.create(cleanup_dir)
    infile <- file.path(cleanup_dir, "in.tif")
    wopt <- list(
      datatype = datatype,
      gdal = sprintf("COMPRESS=%s", compress)
    )
    if (!is.na(nodata)) {
      wopt$NAflag <- nodata
    }
    terra::writeRaster(r, infile, overwrite = TRUE, wopt = wopt)
  } else {
    infile <- r
  }

  # -te + -tr force cell-for-cell alignment with the template (more
  # deterministic than letting gdalwarp infer the output grid).
  opts <- c(
    "-overwrite",
    "-t_srs",
    tmpl_crs_str,
    "-r",
    method,
    "-ot",
    ot,
    "-tr",
    as.character(tmpl_res[1]),
    as.character(tmpl_res[2]),
    "-te",
    as.character(tmpl_ext["xmin"]),
    as.character(tmpl_ext["ymin"]),
    as.character(tmpl_ext["xmax"]),
    as.character(tmpl_ext["ymax"]),
    "-wo",
    sprintf("NUM_THREADS=%s", threads),
    "-co",
    "TILED=YES",
    "-co",
    sprintf("COMPRESS=%s", compress),
    "-co",
    "BIGTIFF=IF_SAFER"
  )
  if (!is.na(nodata)) {
    opts <- c(
      opts,
      "-srcnodata",
      as.character(nodata),
      "-dstnodata",
      as.character(nodata)
    )
  }

  sf::gdal_utils(
    util = "warp",
    source = infile,
    destination = outfile,
    options = opts
  )

  # Copy RAT / band-metadata sidecar if present. gdalwarp doesn't move
  # .aux.xml, but pixel values are unchanged by reprojection so the original
  # RAT (and any layer-name / time metadata terra writes) stays valid.
  src_aux <- paste0(infile, ".aux.xml")
  out_aux <- paste0(outfile, ".aux.xml")
  if (file.exists(src_aux)) {
    file.copy(src_aux, out_aux, overwrite = TRUE)
  }

  # Clean up temp input
  if (!is.null(cleanup_dir)) {
    unlink(cleanup_dir, recursive = TRUE)
  }

  # Verify alignment
  warped <- terra::rast(outfile)
  if (verify && !terra::compareGeom(template, warped, stopOnError = FALSE)) {
    warning(
      "Output does NOT match template grid; downstream stacking will fail."
    )
  }

  warped
}

