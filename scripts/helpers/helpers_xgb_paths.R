# Path simulation from fitted xgb models -------------------------------------
#
# The xgb path's redistribution kernel and the path simulator built on it. A
# model from fit_final_xgb.R scores every 30 m cell within reach of the
# deer. A cell's weight is exp(score) divided by its distance from the deer,
# and the deer's next position is drawn from those weights.
#
# Dividing by distance matches the kernel to how the model was fitted. The
# fit's random points were spread evenly in distance (nonp_random_points()),
# so the trees score a cell against points that crowd near the start. A
# grid of cells is spread evenly in area instead, with more cells in each
# ring the further out the ring lies. The deer's own cell gets no weight,
# so every step leaves it, as in the GAM kernel.
#
# The kernel reads every column the way the wrangle extracted it for the
# fit: the same rasters, cropped to the same window, read at the cell under
# the point. xgb_path_check() confirms that on a deer's own observed and
# random steps before any path is drawn.
#
# The hr and hab blocks read only the cell under the endpoint, so their
# score is computed once per cell (xgb_path_static()). Only the move and
# modifier blocks, which read the step itself, are scored at every step.
#
# Part of the helper library split out of scripts/helper_functions.R, which
# now sources every file in this folder. Scripts keep sourcing that one
# aggregator, so nothing here needs to be sourced directly.

#' Load a model saved by fit_final_xgb.R, ready to predict
#'
#' @param path File written by fit_final_xgb.R
#' @param nthread Threads per prediction
#' @return The saved list, with `boosters` loaded and `boosters` and `specs`
#'   named by block
xgb_path_model <- function(path, nthread = 1L) {
  m <- readRDS(path)
  m$boosters <- lapply(m$boosters, function(raw) {
    bst <- xgboost::xgb.load.raw(raw)
    xgboost::xgb.model.parameters(bst) <- list(nthread = nthread)
    bst
  })
  names(m$specs) <- vapply(m$specs, `[[`, character(1), "name")
  names(m$boosters) <- names(m$specs)
  m
}

#' Score of one block
#'
#' @param model Model from xgb_path_model()
#' @param block Block name: "move", "modifier", "hr" or "hab"
#' @param d Data frame carrying the block's columns
#' @return Score per row; NA where any of the block's columns is missing,
#'   because the fit never saw a missing value
xgb_path_block <- function(model, block, d) {
  spec <- model$specs[[block]]
  X <- as.matrix(d[, spec$feats, drop = FALSE])
  out <- rep(NA_real_, nrow(X))
  ok <- stats::complete.cases(X)
  if (any(ok)) {
    bst <- model$boosters[[block]]
    dm <- xgb_matrix(X[ok, , drop = FALSE], spec$categorical)
    out[ok] <- xgb_raw(bst, dm, 1, xgboost::xgb.get.num.boosted.rounds(bst))
  }
  out
}

#' A deer-year's landscape, as the wrangle read it
#'
#' Crops the rasters to the window extract_step_variables() used for this
#' deer-year (deer_crop_buffer() around the start of every step with random
#' points) and keeps each column the models read as one value per cell.
#' Landcover becomes the 0-based codes of LANDCOVER_LEVELS, NA on open
#' water. Cells are numbered as terra numbers them, row by row from the top
#' left.
#'
#' @param track One row of data/tracks/ (id, season, year, stp,
#'   stp.var.nonp)
#' @param landcover load_landcover() for the deer's year and season
#' @param ndvi load_ndvi() for the deer's year
#' @param topo load_topo()
#' @return list(ext, nrow, ncol, res, cols, ndvi, ndvi_times): `cols` holds
#'   one vector per column, `ndvi` one matrix column per monthly layer
xgb_path_landscape <- function(track, landcover, ndvi, topo) {
  window <- step_crop_window(track$stp.var.nonp[[1]],
                             deer_crop_buffer(track$stp[[1]]$sl_))
  lc <- terra::crop(landcover, window)
  nd <- terra::crop(ndvi, window)
  tp <- terra::crop(topo, window)
  hr <- load_hr_center_raster(track$id, track$season, track$year, lc)
  for (r in list(nd, tp, hr)) {
    stopifnot(terra::compareGeom(lc, r, stopOnError = FALSE))
  }
  cats <- terra::levels(lc[["wiscland"]])[[1]]
  value <- terra::values(lc[["wiscland"]], mat = FALSE)
  label <- cats[[2]][match(value, cats[[1]])]
  list(
    ext = as.vector(terra::ext(lc)),
    nrow = terra::nrow(lc),
    ncol = terra::ncol(lc),
    res = terra::res(lc)[1],
    cols = list(
      landcover = match(label, LANDCOVER_LEVELS) - 1L,
      forest_edge = terra::values(lc[["forest_edge"]], mat = FALSE),
      elevation = terra::values(tp[["elevation"]], mat = FALSE),
      northness = terra::values(tp[["northness"]], mat = FALSE),
      eastness = terra::values(tp[["eastness"]], mat = FALSE),
      HR_center = terra::values(hr, mat = FALSE)
    ),
    ndvi = terra::values(nd, mat = TRUE),
    ndvi_times = terra::time(nd)
  )
}

#' Cell under each point
#'
#' @param land Landscape from xgb_path_landscape()
#' @param x,y Coordinates
#' @return Cell number, NA off the window
xgb_path_cell <- function(land, x, y) {
  col <- floor((x - land$ext[["xmin"]]) / land$res) + 1
  row <- floor((land$ext[["ymax"]] - y) / land$res) + 1
  inside <- col >= 1 & col <= land$ncol & row >= 1 & row <= land$nrow
  ifelse(inside, (row - 1) * land$ncol + col, NA_real_)
}

#' NDVI layer for each time, as amt::extract_covariates_var_time() picks it
#'
#' The layer whose time stamp is nearest, NA if that is more than
#' `max_days` away; the wrangle allowed 31.
#'
#' @param t Times
#' @param layer_times Time stamps of the NDVI layers
#' @param max_days Largest gap allowed
#' @return Layer index per time
xgb_path_ndvi_layer <- function(t, layer_times, max_days = 31) {
  gap <- abs(outer(as.numeric(t), as.numeric(layer_times), "-"))
  i <- max.col(-gap, ties.method = "first")
  i[gap[cbind(seq_along(i), i)] > max_days * 86400] <- NA
  i
}

#' Clock hour of a step's start, as extract_step_variables() records tod_
#'
#' @param t Time of the step's start, in the tracks' local time zone
#' @return Decimal hour, 0-24
xgb_path_hour <- function(t) {
  lubridate::hour(t) + lubridate::minute(t) / 60
}

#' Endpoint columns of some cells
#'
#' @param land Landscape from xgb_path_landscape()
#' @param cell Cell numbers
#' @param layer NDVI layer for each cell
#' @return Data frame with the endpoint columns of the pooled tables
xgb_path_end_columns <- function(land, cell, layer) {
  data.frame(
    landcover = land$cols$landcover[cell],
    HR_center_end = land$cols$HR_center[cell],
    ndvi_end = land$ndvi[cbind(cell, layer)],
    forest_edge_end = land$cols$forest_edge[cell],
    elevation_end = land$cols$elevation[cell],
    northness_end = land$cols$northness[cell],
    eastness_end = land$cols$eastness[cell]
  )
}

#' Start-of-step columns of some cells
#'
#' @param land Landscape from xgb_path_landscape()
#' @param cell Cell numbers of the steps' starts
#' @param layer NDVI layer for each start
#' @return Data frame with the start columns of the pooled tables
xgb_path_start_columns <- function(land, cell, layer) {
  data.frame(
    wiscland_start = land$cols$landcover[cell],
    ndvi_start = land$ndvi[cbind(cell, layer)],
    forest_edge_start = land$cols$forest_edge[cell]
  )
}

#' Score of the hr and hab blocks in every cell of the window
#'
#' Those blocks read only the cell and the deer, so they are scored once
#' per NDVI layer rather than at every step. A cell with any of their
#' columns missing, open water among them, scores NA.
#'
#' @param model Model from xgb_path_model()
#' @param land Landscape from xgb_path_landscape()
#' @param layers NDVI layers to score
#' @param deer_id The deer's id code, for the deer types
#' @return Matrix, one row per cell and one column per NDVI layer, NA in
#'   the layers not asked for. Attribute "reads_ndvi": FALSE when no block
#'   reads NDVI, and then every column holds the same scores.
xgb_path_static <- function(model, land, layers, deer_id = NA_integer_) {
  blocks <- intersect(c("hr", "hab"), names(model$specs))
  reads_ndvi <- any(vapply(model$specs[blocks],
                           function(s) "ndvi_end" %in% s$feats, logical(1)))
  n_cell <- land$nrow * land$ncol
  out <- matrix(NA_real_, n_cell, ncol(land$ndvi))
  for (layer in if (reads_ndvi) layers else 1L) {
    d <- xgb_path_end_columns(land, seq_len(n_cell), layer)
    d$deer_id <- deer_id
    out[, layer] <- Reduce(`+`, lapply(blocks, xgb_path_block,
                                       model = model, d = d))
  }
  if (!reads_ndvi) {
    out[] <- out[, 1]
  }
  attr(out, "reads_ndvi") <- reads_ndvi
  out
}

#' Score of the move and modifier blocks
#'
#' @param model Model from xgb_path_model()
#' @param d Data frame with sl_, tod_hour, cos_ta, the start columns and
#'   deer_id
#' @return Score per row
xgb_path_dynamic <- function(model, d) {
  blocks <- intersect(c("move", "modifier"), names(model$specs))
  Reduce(`+`, lapply(blocks, xgb_path_block, model = model, d = d))
}

#' Check the kernel reads a deer-year's own steps as the fit read them
#'
#' On every observed and random endpoint of the deer's nonp steps, the
#' columns the kernel builds must equal the ones the wrangle stored, and
#' each model's kernel score (xgb_path_static() + xgb_path_dynamic()) must
#' equal its score on the stored columns. Stops on any difference.
#'
#' @param track One row of data/tracks/
#' @param land Landscape from xgb_path_landscape()
#' @param models Named list of models from xgb_path_model()
#' @param statics Matching list of xgb_path_static() results
#' @param deer_ids Matching list of deer id codes (NA for no-deer types)
#' @return Invisibly, the largest score difference per model
xgb_path_check <- function(track, land, models, statics, deer_ids) {
  v <- as.data.frame(track$stp.var.nonp[[1]])
  stp <- track$stp[[1]]
  # Heading of the step before: every step with random points has one in
  # its own burst, ending where this step starts
  r <- match(as.numeric(v$t1_), as.numeric(stp$t1_)) - 1
  stopifnot(!anyNA(r), all(stp$burst_[r] == v$burst_),
            all(stp$x2_[r] == v$x1_), all(stp$y2_[r] == v$y1_))
  heading <- atan2(stp$y2_[r] - stp$y1_[r], stp$x2_[r] - stp$x1_[r])

  end_cell <- xgb_path_cell(land, v$x2_, v$y2_)
  start_cell <- xgb_path_cell(land, v$x1_, v$y1_)
  layer_end <- xgb_path_ndvi_layer(v$t2_, land$ndvi_times)
  layer_start <- xgb_path_ndvi_layer(v$t1_, land$ndvi_times)
  kernel <- cbind(
    xgb_path_end_columns(land, end_cell, layer_end),
    xgb_path_start_columns(land, start_cell, layer_start),
    sl_ = sqrt((v$x2_ - v$x1_)^2 + (v$y2_ - v$y1_)^2),
    cos_ta = cos(atan2(v$y2_ - v$y1_, v$x2_ - v$x1_) - heading),
    tod_hour = xgb_path_hour(v$t1_)
  )
  code <- function(x) as.integer(factor(x, levels = LANDCOVER_LEVELS)) - 1L
  stored <- data.frame(
    landcover = code(v$wiscland_end),
    HR_center_end = v$HR_center_end,
    ndvi_end = v$ndvi_end,
    forest_edge_end = v$forest_edge_end,
    elevation_end = v$elevation_end,
    northness_end = v$northness_end,
    eastness_end = v$eastness_end,
    wiscland_start = code(v$wiscland_start),
    ndvi_start = v$ndvi_start,
    forest_edge_start = v$forest_edge_start,
    sl_ = v$sl_,
    cos_ta = cos(v$ta_),
    tod_hour = v$tod_
  )
  same <- function(a, b) {
    (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= 1e-6)
  }
  differ <- vapply(names(stored), function(nm) {
    sum(!same(kernel[[nm]], stored[[nm]]))
  }, numeric(1))
  if (any(differ > 0)) {
    stop(sprintf("The kernel reads columns differently from the wrangle: %s",
                 paste(sprintf("%s (%d rows)", names(differ)[differ > 0],
                               differ[differ > 0]), collapse = ", ")))
  }

  gaps <- vapply(names(models), function(nm) {
    m <- models[[nm]]
    feats <- unique(unlist(lapply(m$specs, `[[`, "feats")))
    kernel$deer_id <- deer_ids[[nm]]
    stored$deer_id <- deer_ids[[nm]]
    rows <- stats::complete.cases(stored[, feats])
    # A tibble, as in the fit: one-column blocks keep their column name
    reference <- Reduce(`+`, xgb_contributions(
      m$boosters, tibble::as_tibble(stored[rows, ]), m$specs
    ))
    static <- statics[[nm]]
    column <- if (attr(static, "reads_ndvi")) layer_end else 1L
    score <- static[cbind(end_cell, column)][rows] +
      xgb_path_dynamic(m, kernel[rows, ])
    max(abs(score - reference))
  }, numeric(1))
  if (any(!is.finite(gaps) | gaps > 1e-4)) {
    stop(sprintf("Kernel scores differ from the fitted scores: %s",
                 paste(sprintf("%s %.3g", names(gaps), gaps),
                       collapse = ", ")))
  }
  invisible(gaps)
}

#' Offsets to every cell a step can reach
#'
#' Every cell whose centre can lie within `reach` of a point anywhere in the
#' deer's own cell, less that cell itself; the exact distance is taken per
#' step.
#'
#' @param reach Longest step allowed (m)
#' @param res Cell size (m)
#' @return Data frame of column (dc) and row (dr) offsets
xgb_path_offsets <- function(reach, res) {
  k <- ceiling(reach / res) + 1
  off <- expand.grid(dc = -k:k, dr = -k:k)
  off[sqrt(off$dc^2 + off$dr^2) * res <= reach + res &
        (off$dc != 0 | off$dr != 0), ]
}

#' Simulate paths from one model
#'
#' Every path starts at `start` with a random heading and takes one step
#' per interval of `times`, through any gaps in the real fixes. At each
#' step every cell whose centre lies within `reach` of the deer is a
#' candidate, except the deer's own cell. Cells off the window, or with a
#' column missing (open water among them), get no weight. The others weigh
#' exp(score) / distance, and one is drawn. A path whose candidates all
#' have no weight stops there, its later positions NA.
#'
#' @param model Model from xgb_path_model()
#' @param land Landscape from xgb_path_landscape()
#' @param static Cell scores from xgb_path_static()
#' @param start x and y of the first fix
#' @param times Fix times, from the first fix; one step per interval
#' @param reach Longest step allowed (m)
#' @param n_paths Paths to draw
#' @param deer_id The deer's id code, for the deer types
#' @return tibble(nsim, x_, y_, t_), every path at every time. Attributes:
#'   "off_window", per step and path the share of the cells within reach
#'   that lie off the window; "stopped", the paths that stopped early.
xgb_simulate_paths <- function(model, land, static, start, times, reach,
                               n_paths, deer_id = NA_integer_) {
  res <- land$res
  off <- xgb_path_offsets(reach, res)
  n_off <- nrow(off)
  n_steps <- length(times) - 1
  x <- matrix(NA_real_, n_steps + 1, n_paths)
  y <- x
  x[1, ] <- start[1]
  y[1, ] <- start[2]
  heading <- stats::runif(n_paths, -pi, pi)
  off_window <- matrix(NA_real_, n_steps, n_paths)
  alive <- rep(TRUE, n_paths)
  end_column <- function(s) {
    if (attr(static, "reads_ndvi")) {
      xgb_path_ndvi_layer(times[s + 1], land$ndvi_times)
    } else {
      1L
    }
  }

  for (s in seq_len(n_steps)) {
    live <- which(alive)
    if (!length(live)) {
      break
    }
    here <- xgb_path_cell(land, x[s, live], y[s, live])
    stopifnot(!anyNA(here))
    row <- (here - 1) %/% land$ncol + 1
    col <- (here - 1) %% land$ncol + 1

    # Candidate cells for every live path
    p <- rep(live, each = n_off)
    cr <- rep(row, each = n_off) + off$dr
    cc <- rep(col, each = n_off) + off$dc
    cx <- land$ext[["xmin"]] + (cc - 0.5) * res
    cy <- land$ext[["ymax"]] - (cr - 0.5) * res
    dx <- cx - x[s, p]
    dy <- cy - y[s, p]
    sl <- sqrt(dx^2 + dy^2)
    within <- sl <= reach
    inside <- cr >= 1 & cr <= land$nrow & cc >= 1 & cc <= land$ncol
    off_window[s, live] <- 1 - tapply(inside[within],
                                      factor(p[within], levels = live), mean)
    keep <- within & inside
    cell <- (cr[keep] - 1) * land$ncol + cc[keep]
    p <- p[keep]
    cx <- cx[keep]
    cy <- cy[keep]
    dx <- dx[keep]
    dy <- dy[keep]
    sl <- sl[keep]

    score <- static[cbind(cell, end_column(s))]
    j <- match(p, live)
    starts <- xgb_path_start_columns(
      land, here, xgb_path_ndvi_layer(times[s], land$ndvi_times)
    )
    d <- cbind(
      data.frame(sl_ = sl, tod_hour = xgb_path_hour(times[s]),
                 cos_ta = cos(atan2(dy, dx) - heading[p])),
      starts[j, , drop = FALSE],
      deer_id = deer_id
    )
    ok <- !is.na(score)
    score[ok] <- score[ok] + xgb_path_dynamic(model, d[ok, , drop = FALSE])

    for (q in live) {
      i <- which(p == q & !is.na(score))
      if (!length(i)) {
        alive[q] <- FALSE
        next
      }
      w <- exp(score[i] - max(score[i])) / sl[i]
      pick <- i[sample.int(length(i), 1, prob = w)]
      x[s + 1, q] <- cx[pick]
      y[s + 1, q] <- cy[pick]
      heading[q] <- atan2(dy[pick], dx[pick])
    }
  }

  out <- tibble::tibble(
    nsim = rep(seq_len(n_paths), each = n_steps + 1),
    x_ = as.vector(x),
    y_ = as.vector(y),
    t_ = rep(times, n_paths)
  )
  attr(out, "off_window") <- off_window
  attr(out, "stopped") <- which(!alive)
  out
}

#' One-step log score of a deer-year's real steps under one model
#'
#' At every real step the kernel of xgb_simulate_paths() is built at the
#' real start: every cell whose centre lies within `reach`, except the
#' start's own cell, weighs exp(score) / distance, the turning angle taken
#' from the step before. The step's log score is the log of the weight of
#' the cell where the deer ended, over the summed weight of all of them.
#'
#' A step is scored only if it has a turning angle (the first step of a
#' burst has none, and the fit left it out too). A scored step can still
#' come out without a value:
#'   "beyond_reach"  its end cell is not among the candidates: further than
#'                   `reach`, or off the window
#'   "own_cell"      it ends in the cell it started in, which gets no weight
#'   "no_weight"     its end cell has a column missing (open water among
#'                   them) and so no weight
#'
#' @param model Model from xgb_path_model()
#' @param land Landscape from xgb_path_landscape()
#' @param static Cell scores from xgb_path_static()
#' @param stp The deer-year's real steps (x1_, y1_, x2_, y2_, t1_, t2_,
#'   ta_, burst_, sl_)
#' @param reach Longest step allowed (m)
#' @param deer_id The deer's id code, for the deer types
#' @return tibble(step, t1_, sl_, status, logp): one row per step with a
#'   turning angle; logp is NA unless status is "ok"
xgb_onestep_logscore <- function(model, land, static, stp, reach,
                                 deer_id = NA_integer_) {
  res <- land$res
  off <- xgb_path_offsets(reach, res)
  steps <- which(!is.na(stp$ta_))
  # The step before ends where this one starts, in the same burst
  prev <- steps - 1L
  stopifnot(all(prev >= 1), all(stp$burst_[prev] == stp$burst_[steps]),
            all(stp$x2_[prev] == stp$x1_[steps]),
            all(stp$y2_[prev] == stp$y1_[steps]))
  heading <- atan2(stp$y2_[prev] - stp$y1_[prev],
                   stp$x2_[prev] - stp$x1_[prev])
  layer_start <- xgb_path_ndvi_layer(stp$t1_[steps], land$ndvi_times)
  layer_end <- if (attr(static, "reads_ndvi")) {
    xgb_path_ndvi_layer(stp$t2_[steps], land$ndvi_times)
  } else {
    rep(1L, length(steps))
  }

  out <- lapply(seq_along(steps), function(j) {
    s <- steps[j]
    x0 <- stp$x1_[s]
    y0 <- stp$y1_[s]
    here <- xgb_path_cell(land, x0, y0)
    end <- xgb_path_cell(land, stp$x2_[s], stp$y2_[s])
    row <- (here - 1) %/% land$ncol + 1
    col <- (here - 1) %% land$ncol + 1
    cr <- row + off$dr
    cc <- col + off$dc
    cx <- land$ext[["xmin"]] + (cc - 0.5) * res
    cy <- land$ext[["ymax"]] - (cr - 0.5) * res
    dx <- cx - x0
    dy <- cy - y0
    sl <- sqrt(dx^2 + dy^2)
    keep <- sl <= reach & cr >= 1 & cr <= land$nrow & cc >= 1 &
      cc <= land$ncol
    cell <- (cr[keep] - 1) * land$ncol + cc[keep]
    status <- if (is.na(end) || !end %in% c(here, cell)) {
      "beyond_reach"
    } else if (end == here) {
      "own_cell"
    } else {
      "ok"
    }
    logp <- NA_real_
    if (status == "ok") {
      d <- cbind(
        data.frame(sl_ = sl[keep], tod_hour = xgb_path_hour(stp$t1_[s]),
                   cos_ta = cos(atan2(dy[keep], dx[keep]) - heading[j])),
        xgb_path_start_columns(land, here, layer_start[j])[
          rep(1L, sum(keep)), , drop = FALSE
        ],
        deer_id = deer_id
      )
      score <- static[cbind(cell, layer_end[j])]
      ok <- !is.na(score)
      score[ok] <- score[ok] + xgb_path_dynamic(model, d[ok, , drop = FALSE])
      if (is.na(score[cell == end])) {
        status <- "no_weight"
      } else {
        lw <- score[ok] - log(sl[keep][ok])
        top <- max(lw)
        logp <- lw[cell[ok] == end] - (top + log(sum(exp(lw - top))))
      }
    }
    tibble::tibble(step = s, t1_ = stp$t1_[s], sl_ = stp$sl_[s],
                   status = status, logp = logp)
  })
  dplyr::bind_rows(out)
}

#' The simulation clock for one deer-year
#'
#' A step every `step_hours` from the first real fix, through the gaps in
#' the fixes. Some collars shift their schedule partway, so a later fix
#' falls between two ticks of the clock. Where a fix is more than
#' `tolerance_mins` off the clock, the step that crosses the gap before it
#' is stretched to end exactly at that fix (a 6-hour gap is one 6-hour
#' step, a 10-hour gap a 4-hour then a 6-hour step), and the clock carries
#' on from there. A deer-year whose fixes all sit on the clock gets the
#' plain grid, the first fix plus whole steps.
#'
#' @param fixes Real fix times, sorted
#' @param step_hours Hours per step
#' @param tolerance_mins How far a fix may sit from the clock and still
#'   count as on it
#' @return The clock's times, from the first fix to the last
xgb_path_clock <- function(fixes, step_hours = 4, tolerance_mins = 30) {
  step <- step_hours * 3600
  t <- as.numeric(fixes)
  out <- t[1]
  now <- t[1]
  for (f in t[-1]) {
    gap <- f - now
    k <- round(gap / step)
    if (k >= 1 && abs(gap - k * step) <= tolerance_mins * 60) {
      out <- c(out, now + step * seq_len(k))
      now <- now + step * k
    } else {
      n <- max(1, floor(gap / step))
      out <- c(out, now + step * seq_len(n - 1), f)
      now <- f
    }
  }
  as.POSIXct(out, origin = "1970-01-01", tz = lubridate::tz(fixes))
}

#' The canonical training deer-years, read once per session
#'
#' canonical_deer() of the full deer list, the training set only, with the
#' columns the path helpers read (id, season, year and the observed steps).
#' The list is read on the first call and kept for the rest of the R
#' session.
#'
#' @return Data frame, one row per training deer-year
xgb_train_deer <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      d <- canonical_deer(readRDS("library/SW_filtered_deer.RData"))
      cache <<- d[d$set == "train", c("id", "season", "year", "stp")]
    }
    cache
  }
})

#' Reach of a step for one animal in one season
#'
#' The animal's longest step with a turning angle over its training years
#' of the season: the furthest any random point reached for it in the fit.
#'
#' @param id Animal id
#' @param season Season code
#' @return Reach in metres
xgb_path_reach <- function(id, season) {
  train <- xgb_train_deer()
  train <- train[train$id == id & train$season == season, ]
  stopifnot(nrow(train) > 0)
  max(vapply(train$stp, function(s) max(s$sl_[!is.na(s$ta_)]), numeric(1)))
}

#' Folder of the simulated paths for a set and reach rule
#'
#' @param set "test" or "train"
#' @param reach "longest" or "gamma99" (see sim_paths_xgb.R)
#' @return Folder path
xgb_sims_dir <- function(set, reach = "longest") {
  stopifnot(set %in% c("test", "train"),
            reach %in% c("longest", "gamma99"))
  base <- if (set == "test") "sims/xgb_paths" else "sims/xgb_paths_train"
  if (reach == "longest") base else paste0(base, "_", reach)
}

#' Tag of the scores for a set and reach rule: the set, plus the rule
#' unless it is "longest"
#'
#' @inheritParams xgb_sims_dir
#' @return e.g. "train" or "train_gamma99"
xgb_scores_tag <- function(set, reach = "longest") {
  stopifnot(set %in% c("test", "train"),
            reach %in% c("longest", "gamma99"))
  if (reach == "longest") set else paste0(set, "_", reach)
}

#' amt's reach for a deer-year: a quantile of its tentative step length
#'
#' The gamma that amt::random_steps() fitted to the deer-year's step lengths
#' when the wrangle drew the amt random points, stored on the track as
#' attr(stp.var, "sl_"). Its 0.99 quantile, rounded up to the metre, is the
#' max.dist the amt and GAM kernels stop at. Used here only to set how far a
#' step can reach; the xgb kernel itself has no step-length distribution.
#'
#' @param track One deer-year from data/tracks/
#' @param p Quantile of the tentative gamma
#' @return Reach in metres
xgb_path_reach_gamma <- function(track, p = 0.99) {
  d <- attr(track$stp.var[[1]], "sl_")
  stopifnot(identical(d$name, "gamma"))
  ceiling(stats::qgamma(p, shape = d$params$shape, scale = d$params$scale))
}

#' The training deer-year closest in time, same animal and season
#'
#' For a training deer-year, itself. For a test deer-year, the training
#' year of the same animal and season nearest to it; a tie goes to the later
#' year.
#'
#' @param id Animal id
#' @param season Season code
#' @param year Year of the deer-year
#' @return Key (<id>_<season>_<year>) of that training deer-year
xgb_closest_train_key <- function(id, season, year) {
  train <- xgb_train_deer()
  train <- train[train$id == id & train$season == season, ]
  stopifnot(nrow(train) > 0)
  years <- as.integer(train$year)
  closest <- years[order(abs(years - year), -years)][1]
  sprintf("%s_%s_%d", id, season, closest)
}

#' Reach of a deer-year's steps under a reach rule
#'
#' "longest" is the animal's longest training step in the season
#' (xgb_path_reach()). "gamma99" is the smaller of that and the 0.99
#' quantile of the tentative gamma of the training deer-year closest in time
#' (xgb_closest_train_key(), xgb_path_reach_gamma()): the deer-year itself in
#' training, and for a test deer-year a training year of the same animal and
#' season, as the reach of a year not yet seen has to come from the ones
#' that were.
#'
#' @param track One deer-year from data/tracks/
#' @param rule "longest" or "gamma99"
#' @return Reach in metres
xgb_path_reach_rule <- function(track, rule = c("longest", "gamma99")) {
  rule <- match.arg(rule)
  longest <- xgb_path_reach(track$id, track$season)
  if (rule == "longest") {
    return(longest)
  }
  key <- xgb_closest_train_key(track$id, track$season,
                               as.integer(track$year))
  source_track <- readRDS(sprintf("data/tracks/data_%s.rds", key))
  min(xgb_path_reach_gamma(source_track), longest)
}

# Per-deer-year GAMs scored against xgb ---------------------------------------
#
# The xgb log-score filters set each selected model against GAMs fitted to
# the deer-year alone: the GAM null (fit_null_gam_xgb.R) and the Gilbertson
# et al. (2025) habitat model (fit_gilbertson_gam_xgb.R). Both are fitted on
# the wrangle's amt random points and scored on the cells, reach and
# statuses of xgb_onestep_logscore().

# Endpoint columns a per-deer-year GAM may read, as xgb_path_end_columns()
# names them; landcover enters the GAMs as the factor wiscland_end
XGB_GAM_END_VARS <- c("HR_center_end", "ndvi_end", "forest_edge_end",
                      "elevation_end", "northness_end", "eastness_end",
                      "wiscland_end")

#' Endpoint columns a fitted per-deer-year GAM reads
#'
#' A GAM that reads NDVI centred on the deer-year's mean (ndvi_c, with the
#' centre stored on the fit as `ndvi_centre`) reads ndvi_end.
#'
#' @param fit Fitted mgcv gam
#' @return The entries of XGB_GAM_END_VARS it reads
xgb_gam_end_vars <- function(fit) {
  vars <- all.vars(stats::formula(fit))
  if ("ndvi_c" %in% vars) {
    stopifnot(is.finite(xgb_gam_ndvi_centre(fit)))
    vars <- c(vars, "ndvi_end")
  }
  intersect(XGB_GAM_END_VARS, vars)
}

#' The NDVI centre of a per-deer-year GAM
#'
#' @param fit Fitted mgcv gam
#' @return fit$ndvi_centre, the mean NDVI its ndvi_c is centred on; NA for
#'   a GAM that does not read ndvi_c
xgb_gam_ndvi_centre <- function(fit) {
  if (is.null(fit$ndvi_centre)) NA_real_ else fit$ndvi_centre
}

#' A deer-year's steps for a per-deer-year GAM, in mgcv cox.ph form
#'
#' The wrangle's amt random points (stp.var: 25 per step, drawn from the
#' deer-year's own gamma and von Mises), `n_pts` of them kept per step at
#' random. The points are picked before any step is dropped, so two models
#' given the same seed are fitted on the same points. They are drawn with
#' the L'Ecuyer-CMRG generator whatever the caller's, so a seed gives the
#' same points in a furrr worker (seed = TRUE) and in a plain session; the
#' caller's generator and its state are left as they were. A step with any of
#' sl_, ta_, tod_ or `end_vars` missing on any of its rows is dropped whole,
#' as the pooled tables drop them; open water has no landcover class, so a
#' step ending on it drops out when wiscland_end is read.
#'
#' @param track One deer-year from data/tracks/
#' @param n_pts Random points kept per step
#' @param seed Seed for the draw of the points kept
#' @param end_vars Endpoint columns the model reads (XGB_GAM_END_VARS)
#' @return prepare_gam_data() output (times, stratum, obs added), with the
#'   end coordinates and times kept for xgb_gam_check(); wiscland_end, if
#'   read, as a factor with levels LANDCOVER_LEVELS
xgb_gam_data <- function(track, n_pts, seed, end_vars = "HR_center_end") {
  stopifnot(all(end_vars %in% XGB_GAM_END_VARS))
  v <- as.data.frame(track$stp.var[[1]])
  kept <- withr::with_seed(
    seed,
    v |>
      dplyr::filter(!case_) |>
      dplyr::slice_sample(n = n_pts, by = step_id_),
    .rng_kind = "L'Ecuyer-CMRG",
    .rng_normal_kind = "Inversion",
    .rng_sample_kind = "Rejection"
  )
  v <- dplyr::bind_rows(dplyr::filter(v, case_), kept) |>
    dplyr::select(step_id_, case_, x1_, y1_, x2_, y2_, t1_, t2_, sl_, ta_,
                  tod_, dplyr::all_of(end_vars))
  if ("wiscland_end" %in% end_vars) {
    v$wiscland_end <- factor(v$wiscland_end, levels = LANDCOVER_LEVELS)
  }
  v <- v |>
    dplyr::group_by(step_id_) |>
    dplyr::filter(
      !any(dplyr::if_any(dplyr::all_of(c("sl_", "ta_", "tod_", end_vars)),
                         is.na)),
      sum(case_) == 1,
      dplyr::n() >= 2
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(step_id_, dplyr::desc(case_))
  prepare_gam_data(as.data.frame(v))
}

#' Endpoint columns of some cells, as a per-deer-year GAM reads them
#'
#' @param land Landscape from xgb_path_landscape()
#' @param cell Cell numbers
#' @param layer NDVI layer for each cell
#' @param end_vars Endpoint columns the model reads
#' @param ndvi_centre The model's NDVI centre (xgb_gam_ndvi_centre()); NA
#'   for a model that does not read ndvi_c
#' @return Data frame of `end_vars`, and ndvi_c (ndvi_end less the centre)
#'   when the centre is given; wiscland_end a factor with levels
#'   LANDCOVER_LEVELS, NA on open water
xgb_gam_end_columns <- function(land, cell, layer, end_vars,
                                ndvi_centre = NA_real_) {
  d <- xgb_path_end_columns(land, cell, layer)
  d$wiscland_end <- factor(LANDCOVER_LEVELS[d$landcover + 1L],
                           levels = LANDCOVER_LEVELS)
  out <- d[, end_vars, drop = FALSE]
  if (!is.na(ndvi_centre)) {
    out$ndvi_c <- d$ndvi_end - ndvi_centre
  }
  out
}

#' Check a per-deer-year GAM's kernel reads the deer-year as its fit did
#'
#' On every endpoint the model was fitted on that lies in the window, each
#' column the kernel builds must equal the stored one, and none may be on a
#' cell the kernel gives no weight (open water, or a column missing): the
#' fit's random points were drawn off water. Stops on any difference.
#'
#' @param data The fitted rows, from xgb_gam_data()
#' @param land Landscape from xgb_path_landscape()
#' @param end_vars Endpoint columns the model reads
#' @return Invisibly, the largest difference over the numeric columns
xgb_gam_check <- function(data, land, end_vars) {
  cell <- xgb_path_cell(land, data$x2_, data$y2_)
  inside <- !is.na(cell)
  layer <- if ("ndvi_end" %in% end_vars) {
    xgb_path_ndvi_layer(data$t2_[inside], land$ndvi_times)
  } else {
    rep(1L, sum(inside))
  }
  k <- xgb_gam_end_columns(land, cell[inside], layer, end_vars)
  wet <- is.na(land$cols$landcover[cell[inside]]) |
    !stats::complete.cases(k)
  if (any(wet)) {
    stop(sprintf("%d fitted endpoints lie on cells with no weight",
                 sum(wet)))
  }
  gaps <- vapply(end_vars, function(v) {
    if (v == "wiscland_end") {
      as.numeric(any(as.character(k[[v]]) !=
                       as.character(data[[v]][inside])))
    } else {
      max(abs(k[[v]] - data[[v]][inside]))
    }
  }, numeric(1))
  if (any(!is.finite(gaps) | gaps > 1e-6)) {
    stop(sprintf("Kernel columns differ from the fitted: %s",
                 paste(sprintf("%s %.3g", end_vars, gaps), collapse = ", ")))
  }
  invisible(max(gaps))
}

#' One-step log score of a deer-year's real steps under a per-deer-year GAM
#'
#' The cells, reach and statuses of xgb_onestep_logscore(): every cell whose
#' centre lies within `reach` of the real start, except the start's own
#' cell, the turning angle taken from the step before. The GAM was fitted
#' on random points drawn from the deer-year's gamma and von Mises, so its
#' score corrects those: a cell weighs
#'   exp(eta + log tentative density) / distance,
#' as in gam_kernel_weights(). eta is the GAM's linear predictor and the
#' tentative density comes from gam_movement_kernel(). The endpoint columns
#' are read as xgb_path_end_columns() reads them, NDVI from the layer of the
#' real step's end time. Open water, where landcover is NA, gets no weight,
#' nor does a cell with any column the GAM reads missing.
#'
#' @param fit The deer-year's fitted GAM (mgcv cox.ph)
#' @param land Landscape from xgb_path_landscape()
#' @param stp The deer-year's real steps
#' @param reach Longest step allowed (m)
#' @param sl_distr,ta_distr The tentative distributions, attr(stp.var, "sl_")
#'   and attr(stp.var, "ta_")
#' @return tibble(step, t1_, sl_, status, logp), as xgb_onestep_logscore()
xgb_gam_onestep_logscore <- function(fit, land, stp, reach, sl_distr,
                                     ta_distr) {
  res <- land$res
  off <- xgb_path_offsets(reach, res)
  end_vars <- xgb_gam_end_vars(fit)
  steps <- which(!is.na(stp$ta_))
  prev <- steps - 1L
  stopifnot(all(prev >= 1), all(stp$burst_[prev] == stp$burst_[steps]),
            all(stp$x2_[prev] == stp$x1_[steps]),
            all(stp$y2_[prev] == stp$y1_[steps]))
  heading <- atan2(stp$y2_[prev] - stp$y1_[prev],
                   stp$x2_[prev] - stp$x1_[prev])
  layer_end <- if ("ndvi_end" %in% end_vars) {
    xgb_path_ndvi_layer(stp$t2_[steps], land$ndvi_times)
  } else {
    rep(1L, length(steps))
  }

  out <- lapply(seq_along(steps), function(j) {
    s <- steps[j]
    x0 <- stp$x1_[s]
    y0 <- stp$y1_[s]
    here <- xgb_path_cell(land, x0, y0)
    end <- xgb_path_cell(land, stp$x2_[s], stp$y2_[s])
    row <- (here - 1) %/% land$ncol + 1
    col <- (here - 1) %% land$ncol + 1
    cr <- row + off$dr
    cc <- col + off$dc
    cx <- land$ext[["xmin"]] + (cc - 0.5) * res
    cy <- land$ext[["ymax"]] - (cr - 0.5) * res
    dx <- cx - x0
    dy <- cy - y0
    sl <- sqrt(dx^2 + dy^2)
    keep <- sl <= reach & cr >= 1 & cr <= land$nrow & cc >= 1 &
      cc <= land$ncol
    cell <- (cr[keep] - 1) * land$ncol + cc[keep]
    status <- if (is.na(end) || !end %in% c(here, cell)) {
      "beyond_reach"
    } else if (end == here) {
      "own_cell"
    } else {
      "ok"
    }
    logp <- NA_real_
    if (status == "ok") {
      k <- xgb_gam_end_columns(land, cell, rep(layer_end[j], length(cell)),
                               end_vars, xgb_gam_ndvi_centre(fit))
      dry <- !is.na(land$cols$landcover[cell]) & stats::complete.cases(k)
      if (!dry[cell == end]) {
        status <- "no_weight"
      } else {
        ta <- atan2(dy[keep][dry], dx[keep][dry]) - heading[j]
        nd <- data.frame(
          sl_ = sl[keep][dry],
          ta_ = atan2(sin(ta), cos(ta)),
          tod_ = xgb_path_hour(stp$t1_[s]),
          k[dry, , drop = FALSE]
        )
        eta <- as.numeric(stats::predict(fit, newdata = nd, type = "link"))
        lw <- eta + gam_movement_kernel(nd, sl_distr, ta_distr) -
          log(nd$sl_)
        top <- max(lw)
        logp <- lw[cell[dry] == end] - (top + log(sum(exp(lw - top))))
      }
    }
    tibble::tibble(step = s, t1_ = stp$t1_[s], sl_ = stp$sl_[s],
                   status = status, logp = logp)
  })
  dplyr::bind_rows(out)
}

#' Simulate paths from a per-deer-year GAM
#'
#' xgb_simulate_paths() with the GAM's kernel. Every path starts at `start`
#' with a random heading and takes one step per interval of `times`, through
#' any gaps in the real fixes. At each step every cell whose centre lies
#' within `reach` of the deer is a candidate, except the deer's own cell,
#' and weighs exp(eta + log tentative density) / distance, as in
#' xgb_gam_onestep_logscore(). Cells off the window, on open water, or with
#' a column the GAM reads missing get no weight; NDVI is read from the layer
#' of the step's end time. A path whose candidates all have no weight stops
#' there, its later positions NA.
#'
#' @param fit The fitted GAM (mgcv cox.ph)
#' @param land Landscape from xgb_path_landscape()
#' @param start x and y of the first fix
#' @param times Fix times, from the first fix; one step per interval
#' @param reach Longest step allowed (m)
#' @param n_paths Paths to draw
#' @param sl_distr,ta_distr The GAM's tentative distributions
#' @return tibble(nsim, x_, y_, t_), every path at every time. Attribute
#'   "stopped": the paths that stopped early.
xgb_gam_simulate_paths <- function(fit, land, start, times, reach, n_paths,
                                   sl_distr, ta_distr) {
  res <- land$res
  off <- xgb_path_offsets(reach, res)
  n_off <- nrow(off)
  end_vars <- xgb_gam_end_vars(fit)
  n_steps <- length(times) - 1
  layer_end <- if ("ndvi_end" %in% end_vars) {
    xgb_path_ndvi_layer(times[-1], land$ndvi_times)
  } else {
    rep(1L, n_steps)
  }
  x <- matrix(NA_real_, n_steps + 1, n_paths)
  y <- x
  x[1, ] <- start[1]
  y[1, ] <- start[2]
  heading <- stats::runif(n_paths, -pi, pi)
  alive <- rep(TRUE, n_paths)

  for (s in seq_len(n_steps)) {
    live <- which(alive)
    if (!length(live)) {
      break
    }
    here <- xgb_path_cell(land, x[s, live], y[s, live])
    stopifnot(!anyNA(here))
    row <- (here - 1) %/% land$ncol + 1
    col <- (here - 1) %% land$ncol + 1

    # Candidate cells for every live path
    p <- rep(live, each = n_off)
    cr <- rep(row, each = n_off) + off$dr
    cc <- rep(col, each = n_off) + off$dc
    cx <- land$ext[["xmin"]] + (cc - 0.5) * res
    cy <- land$ext[["ymax"]] - (cr - 0.5) * res
    dx <- cx - x[s, p]
    dy <- cy - y[s, p]
    sl <- sqrt(dx^2 + dy^2)
    keep <- sl <= reach & cr >= 1 & cr <= land$nrow & cc >= 1 &
      cc <= land$ncol
    cell <- (cr[keep] - 1) * land$ncol + cc[keep]
    p <- p[keep]
    cx <- cx[keep]
    cy <- cy[keep]
    dx <- dx[keep]
    dy <- dy[keep]
    sl <- sl[keep]

    k <- xgb_gam_end_columns(land, cell, rep(layer_end[s], length(cell)),
                             end_vars, xgb_gam_ndvi_centre(fit))
    dry <- !is.na(land$cols$landcover[cell]) & stats::complete.cases(k)
    lw <- rep(NA_real_, length(cell))
    if (any(dry)) {
      ta <- atan2(dy[dry], dx[dry]) - heading[p[dry]]
      nd <- data.frame(sl_ = sl[dry], ta_ = atan2(sin(ta), cos(ta)),
                       tod_ = xgb_path_hour(times[s]),
                       k[dry, , drop = FALSE])
      eta <- as.numeric(stats::predict(fit, newdata = nd, type = "link"))
      lw[dry] <- eta + gam_movement_kernel(nd, sl_distr, ta_distr) -
        log(nd$sl_)
    }

    for (q in live) {
      i <- which(p == q & !is.na(lw))
      if (!length(i)) {
        alive[q] <- FALSE
        next
      }
      w <- exp(lw[i] - max(lw[i]))
      pick <- i[sample.int(length(i), 1, prob = w)]
      x[s + 1, q] <- cx[pick]
      y[s + 1, q] <- cy[pick]
      heading[q] <- atan2(dy[pick], dx[pick])
    }
  }

  out <- tibble::tibble(
    nsim = rep(seq_len(n_paths), each = n_steps + 1),
    x_ = as.vector(x),
    y_ = as.vector(y),
    t_ = rep(times, n_paths)
  )
  attr(out, "stopped") <- which(!alive)
  out
}
