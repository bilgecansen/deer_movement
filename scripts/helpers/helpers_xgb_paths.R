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

#' Reach of a step for one animal in one season
#'
#' The animal's longest step with a turning angle over its training years
#' of the season: the furthest any random point reached for it in the fit.
#'
#' @param id Animal id
#' @param season Season code
#' @return Reach in metres
xgb_path_reach <- function(id, season) {
  train <- canonical_deer(readRDS("library/SW_filtered_deer.RData"))
  train <- train[train$set == "train" & train$id == id &
                   train$season == season, ]
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
