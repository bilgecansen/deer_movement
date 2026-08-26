# Helper functions -------------------------------------------------------------
#
# Aggregator. The helper library used to be one ~2100-line file; it now lives as
# themed files in scripts/helpers/ and this file just sources all of them.
#
# Every script keeps its existing `source("scripts/helper_functions.R")` line
# and gets the whole library, exactly as before. Nothing needs to know which
# file a given function ended up in, and a helper can be moved between files
# without touching any caller.
#
#   helpers_constants.R  LANDCOVER_LEVELS, CROP_BUFFER_M
#   helpers_raster.R     landcover / NDVI / home-range raster loading
#   helpers_warp.R       warp_to_template (gdalwarp wrapper; prep scripts only)
#   helpers_track.R      random-point generation, step covariate extraction
#   helpers_issf.R       iSSF fitting, coefficient tidying, log score
#   helpers_gam.R GAM fitting, diagnostics, redistribution kernel, log score
#   helpers_simulate.R   simulate_movement (the iSSF / GAM seam)
#   helpers_metrics.R    energy score, SVF, UD overlap (both paths)
#   helpers_plots.R      track plotting / animation
#
# Load order is for readability only: R resolves function and constant
# references when a function is *called*, not when it is defined, so the
# cross-file references (e.g. extract_step_variables -> load_hr_edge_raster,
# gam_cov_fun -> LANDCOVER_LEVELS) work regardless of ordering.
#
# NOTE: source() resolves paths against the working directory, not against this
# file. Every script in this project runs from the repository root, which is why
# the plain relative paths below work — keep running them that way.

# Evaluate the parts into whatever environment THIS file is being sourced into.
# Without local = , a nested source() always writes to the global environment,
# so `source("scripts/helper_functions.R", local = e)` would silently leave `e`
# empty and dump the whole library into globalenv instead.
.helper_env <- environment()

for (.helper_file in c(
  "helpers_constants.R",
  "helpers_raster.R",
  "helpers_warp.R",
  "helpers_track.R",
  "helpers_issf.R",
  "helpers_gam.R",
  "helpers_simulate.R",
  "helpers_metrics.R",
  "helpers_plots.R"
)) {
  .helper_path <- file.path("scripts", "helpers", .helper_file)
  if (!file.exists(.helper_path)) {
    stop(sprintf(
      "Helper file not found: %s (are you running from the project root?)",
      .helper_path
    ))
  }
  source(.helper_path, local = .helper_env)
}
rm(.helper_file, .helper_path, .helper_env)
