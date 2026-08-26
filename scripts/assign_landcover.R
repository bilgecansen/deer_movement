library(terra)

# Build a landcover raster from a Wiscland2 source raster `src`
# from Appendix B in https://p.widencdn.net/8ghipa/Wiscland_2_User_Guide_September_2016
#
# Two mutually exclusive modes:
#   (a) Canned RAT-based level: pass L in 1:4. Uses the embedded raster
#     attribute
#       table (RAT) to map every source pixel value to its level-L class.
#       Returns a factor raster with values 1..n and labels = level-L class
#       descriptions from the RAT. At L=3/4, classes that are not actually
#       subdivided to that level can be set to NA via na_when_unsubdivided =
#       TRUE; otherwise the parent class label is carried forward (matching
#       DNR's published behavior).
#   (b) Custom reclass: pass `rc` as a data.frame with columns
#       `from`, `to`, `id`, `label`. Uses range-based reclassification with
#       half-open intervals [from, to). `id` values must be contiguous integers
#       1..n. Returns a factor raster with values 1..n and the supplied labels.
#
# Mode precedence: if `rc` is supplied, it takes precedence and L is ignored
# (warning issued if both given). If neither is supplied, L defaults to 1.
#
# NA handling: 65535 (Wiscland2's no-data sentinel) is converted to NA
# at the start of every call. Any NaN in the output is normalized to NA at the
# end so that downstream `unique(values(out))` shows a single missing category
# rather than separate NA and NaN rows.
#
# Validation (Mode B only): Mode B fails fast with informative errors on (i)
# missing required columns, (ii) NA in id/label, (iii) inverted/empty ranges
# (from >= to), (iv) non-contiguous ids, (v) inconsistent labels for the same
# id, and (vi) overlapping ranges. Coverage validation (every source value must
# fall in some [from, to) range) is controlled by `na_unmapped`: by default
# (na_unmapped = FALSE), uncovered values cause a hard error listing the
# offending codes -- typo-protection for handcrafted rc tables. Pass na_unmapped
# = TRUE to allow uncovered values to silently map to NA, useful when you want
# to factorize only a few classes of interest from a much larger
# source vocabulary.

#-------------------------------------------------------------------------------
assign_landcover <- function(
  src,
  L = NULL,
  rc = NULL,
  na_when_unsubdivided = FALSE,
  na_unmapped = FALSE
) {
  # Setup: capture RAT and normalize source NA
  # The RAT must be captured BEFORE any terra::classify() call, since classify()
  # drops categorical metadata as a side effect
  rat <- terra::cats(src)[[1]]
  # Convert Wiscland2's 65535 sentinel to NA. Done once, up front, so neither
  # mode has to think about it again
  src <- terra::classify(src, cbind(65535, NA))
  # If both L and rc are supplied, rc wins. Warn so the user knows L was ignored
  # rather than silently dropped.
  if (!is.null(rc) && !is.null(L)) {
    warning(
      "Both `L` and `rc` were supplied; `L` is ignored when `rc` is given.",
      call. = FALSE
    )
  }

  # ============================================================================
  # Mode B: custom reclassification
  # ============================================================================

  if (!is.null(rc)) {
    # Document the interval semantics every time. Half-open intervals [from, to)
    # are easy to forget and easy to misread
    message(
      "Custom reclass uses half-open intervals [from, to): `from` is included, 
      `to` is excluded. ",
      "E.g. c(5000, 5001, 1) matches only the value 5000."
    )

    # Validate input shape
    if (!is.data.frame(rc)) {
      stop("`rc` must be a data.frame with columns: from, to, id, label.")
    }
    required_cols <- c("from", "to", "id", "label")
    missing_cols <- setdiff(required_cols, names(rc))
    if (length(missing_cols) > 0) {
      stop(
        "`rc` is missing required columns: ",
        paste(missing_cols, collapse = ", ")
      )
    }
    if (anyNA(rc$id) || anyNA(rc$label)) {
      stop("`rc$id` and `rc$label` must not contain NA.")
    }

    # Validate range geometry: each row must be a non-empty interval
    # from >= to means an empty range (==) or an inverted one (>); either is
    # almost always a typo and would silently drop pixels.
    if (any(rc$from >= rc$to)) {
      bad <- which(rc$from >= rc$to)
      stop(
        "`rc` has rows where `from >= to` (empty or inverted ranges):\n",
        paste(
          sprintf("  row %d: from=%s, to=%s", bad, rc$from[bad], rc$to[bad]),
          collapse = "\n"
        )
      )
    }

    # Validate ids: must be contiguous integers 1..n Contiguity matters because
    # the output factor will have integer levels 1..n; gaps would produce a
    # sparse factor that's confusing
    unique_ids <- sort(unique(rc$id))
    expected <- seq_len(max(unique_ids))
    if (!identical(as.integer(unique_ids), as.integer(expected))) {
      stop(
        "`rc$id` must be contiguous integers starting at 1. ",
        "Got: ",
        paste(unique_ids, collapse = ", "),
        " (expected: ",
        paste(expected, collapse = ", "),
        ")."
      )
    }

    # Validate id<->label consistency: each id must have exactly one label
    # Multiple rows can share an id (that's how you collapse several ranges
    # into one class), but they must all carry the same label
    id_label <- unique(rc[, c("id", "label")])
    dup_ids <- id_label$id[duplicated(id_label$id)]
    if (length(dup_ids) > 0) {
      offenders <- id_label[id_label$id %in% dup_ids, ]
      offenders <- offenders[order(offenders$id), ]
      stop(
        "`rc` has conflicting labels for the same id:\n",
        paste(
          sprintf("  id=%s -> %s", offenders$id, shQuote(offenders$label)),
          collapse = "\n"
        )
      )
    }

    # Validate non-overlap: no two rows can claim the same source value
    # Two half-open intervals [a, b) and [c, d) overlap iff a < d AND c < b.
    # Strict < is correct because the upper endpoint is excluded, so e.g.
    # [1000, 2000) and [2000, 3000) are cleanly adjacent, not overlapping.
    n <- nrow(rc)
    if (n >= 2) {
      overlaps <- character(0)
      for (i in seq_len(n - 1)) {
        for (j in (i + 1):n) {
          if (rc$from[i] < rc$to[j] && rc$from[j] < rc$to[i]) {
            overlaps <- c(
              overlaps,
              sprintf(
                paste("  rows %d and %d: [%s, %s) (id=%s, %s)",
                      "overlaps [%s, %s) (id=%s, %s)"),
                i,
                j,
                rc$from[i],
                rc$to[i],
                rc$id[i],
                shQuote(rc$label[i]),
                rc$from[j],
                rc$to[j],
                rc$id[j],
                shQuote(rc$label[j])
              )
            )
          }
        }
      }
      if (length(overlaps) > 0) {
        stop("`rc` has overlapping ranges:\n", paste(overlaps, collapse = "\n"))
      }
    }

    # Validate coverage: every source value must fall in some [from, to), unless
    # na_unmapped = TRUE. When na_unmapped = FALSE (default), uncovered source
    # values cause a hard error -- typo-protection for handcrafted rc tables.
    # When na_unmapped = TRUE, uncovered source values silently map to NA via
    # the `others = NA` arg of terra::classify()
    if (!na_unmapped) {
      src_vals <- sort(unique(terra::values(src, mat = FALSE, na.rm = TRUE)))
      in_any_range <- vapply(
        src_vals,
        function(v) any(v >= rc$from & v < rc$to),
        logical(1)
      )
      uncovered <- src_vals[!in_any_range]
      if (length(uncovered) > 0) {
        stop(
          "`rc` does not cover all source values present in the raster. 
          Uncovered codes: ",
          paste(uncovered, collapse = ", "),
          "\n  (To allow uncovered values to silently map to NA, pass 
          na_unmapped = TRUE.)"
        )
      }
    }

    # Apply reclassification right = FALSE, include.lowest = TRUE makes
    # terra::classify() use [from, to) semantics matching the validation above
    # others = NA is defensive: any source value the validator missed becomes NA
    # rather than silently passing through
    rcl_mat <- as.matrix(rc[, c("from", "to", "id")])
    out <- terra::classify(
      src,
      rcl_mat,
      right = FALSE,
      include.lowest = TRUE,
      others = NA
    )
    # Normalize NaN -> NA so the output has a single uniform missing category as
    # terra sometimes emits NaN where NA is expected on numeric output
    out <- terra::classify(out, cbind(NaN, NA))
    # Attach factor labels: integer id -> human-readable label, ordered by
    # id ascending
    lbl <- id_label[order(id_label$id), ]
    levels(out) <- data.frame(value = lbl$id, label = lbl$label)
    names(out) <- "landcover_custom"
    return(out)
  }

  # ============================================================================
  # Mode A: canned RAT-based level
  # ============================================================================

  # Default to L = 1 if the user supplied neither rc nor L
  if (is.null(L)) {
    L <- 1L
  }
  stopifnot(L %in% 1:4)

  # Mode A depends entirely on the RAT. If the source raster has no categorical
  # metadata, Mode A can't function. Therefore fail with a clear message rather
  # than the cryptic downstream error that an empty RAT would produce
  if (is.null(rat) || nrow(rat) == 0) {
    stop(
      "`src` has no embedded categorical metadata (RAT). 
      Mode A requires the Wiscland2 RAT with cls_lvl_* / cls_desc_* columns."
    )
  }

  # Pick the level columns
  code_col <- paste0("cls_lvl_", L)
  desc_col <- paste0("cls_desc_", L)
  to_code <- rat[[code_col]]

  # Optionally NA-out classes that aren't actually subdivided to level L In the
  # RAT, classes that don't exist at finer levels are encoded by repeating the
  # parent code (e.g. "Cranberries" has cls_lvl_3 == cls_lvl_2) When the user
  # wants to *see* where the classification thins out, set those rows to NA so
  # the output raster is missing there. Only meaningful at L >= 3; at L=1 and
  # L=2 every class has a real subdivision, so the equality check would
  # never trigger
  if (na_when_unsubdivided && L >= 3) {
    parent_col <- paste0("cls_lvl_", L - 1)
    to_code[rat[[code_col]] == rat[[parent_col]]] <- NA
  }

  # Build the level-L vocabulary: distinct (code, label) pairs, with new
  # contiguous integer ids 1..n
  # The output raster will use these contiguous ids as cell values;
  # the RAT's original 4-digit codes are dropped on output
  vocab <- data.frame(code = to_code, label = rat[[desc_col]])
  vocab <- unique(vocab[!is.na(vocab$code), ])
  vocab <- vocab[order(vocab$code), ]
  vocab$id <- seq_len(nrow(vocab))
  # Build the reclass table mapping RAT$Value -> new id, then apply
  # match() preserves NAs from to_code (rows that were unsubdivided get id = NA,
  # which classify() turns into NA in the raster)
  to_id <- vocab$id[match(to_code, vocab$code)]
  rcl <- cbind(from = rat$Value, to = to_id)
  out <- terra::classify(src, rcl = rcl, others = NA)
  # Normalize NaN -> NA so the output has a single uniform missing category
  out <- terra::classify(out, cbind(NaN, NA))
  # Attach factor labels: integer id -> level-L description
  levels(out) <- data.frame(value = vocab$id, label = vocab$label)
  names(out) <- sprintf("level%d", L)
  out
}

# Usage Load wiscland2 data which eoncodes pixels for all levels with
# accompanying RAT table this has already been resampled 'EPSG:6610' wisconsin
# transverse mercator projection wiscland_archive_path <-
# "/Volumes/checastaldo/special/wi_data/env/wiscland/" wiscland2 <-
# terra::rast(paste0(wiscland_archive_path, "wiscland2.tif"))

# Create a level 3 wiscland2 raster
# r <- assign_landcover(src = wiscland2, L = 3)

# Create a level 3 wiscland2 raster where pixels without level 3 data are set to
# NA r <- assign_landcover(src = wiscland2, L = 3, na_when_unsubdivided = TRUE)

# Create a custom wiscland2 raster that combines certain categories
# Here the rc input must cover all possible pixel vales
# rc <- tibble::tribble(
#   ~from,  ~to,  ~id,  ~label,
#    5000, 5001,    1,  "water",
#    1000, 2000,    2,  "other",
#    6000, 8001,    2,  "other",
#    2000, 3000,    3,  "agriculture",
#    3110, 3111,    4,  "hay",
#    3000, 3110,    5,  "grassland",
#    3119, 4000,    5,  "grassland",
#    4230, 4240,    6,  "oak",
#    4240, 4241,    7,  "central.hardwoods",
#    4100, 4200,    8,  "other.forest",
#    4210, 4230,    8,  "other.forest",
#    4250, 5000,    8,  "other.forest"
# )
# r <- assign_landcover(src = wiscland2, rc = rc)

# Create a custom wiscland2 raster that combines certain categories
# and excludes all others as NA
# rc <- tibble::tribble(
#   ~from,  ~to,  ~id,  ~label,
#    4000, 4301,    1,  "forest",
#    6400, 6451,    2,  "forested wetlands",
#    6300, 6331,    3,  "lowland scrub/scrub")
# Permissive: rc only needs to cover what you care about; everything else is NA
# r <- assign_landcover(src = wiscland2, rc = rc, na_unmapped = TRUE)
