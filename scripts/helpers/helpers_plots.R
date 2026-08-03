# Track plotting and animation -------------------------------------------------
#
# Interactive/ad-hoc visualisation of observed and simulated deer paths. No
# pipeline script currently calls these (scripts/plot_deer_paths.R builds its
# own figures); they are kept here so they stay available without adding to
# what every run loads.
#
# Part of the helper library split out of scripts/helper_functions.R, which
# now sources every file in this folder. Scripts keep sourcing that one
# aggregator, so nothing here needs to be sourced directly.

#' @description
#' Plot simulated paths vs. actual (test) path for a single deer.
#' Usage: source this file, then call plot_deer_paths(row_no).

plot_deer_paths <- function(
  row_no,
  data_path = "data_deer_1_119.rds",
  sim_dir = "results",
  model = NULL, # NULL = facet across all non-NA models; or integer index
  sim_alpha = 0.25
) {
  # Load deer data and pick the row
  deer_mvt <- readRDS(data_path)
  if (row_no < 1 || row_no > nrow(deer_mvt)) {
    stop(sprintf("row_no %d out of range (1:%d)", row_no, nrow(deer_mvt)))
  }
  deer_row <- deer_mvt[row_no, ]

  # Actual (test) path
  stp_test <- deer_row$stp_test[[1]]
  obs <- stp_test %>%
    dplyr::select(x_ = x1_, y_ = y1_, t_ = t1_, burst_) %>%
    dplyr::mutate(type = "observed")

  # Simulated paths
  sim_file <- file.path(sim_dir, sprintf("results_sim_%d.rds", row_no))
  if (!file.exists(sim_file)) {
    stop(sprintf("Simulation file not found: %s", sim_file))
  }
  results_sim <- readRDS(sim_file)

  # Build a long data frame of simulated paths across models
  sim_df <- purrr::imap_dfr(results_sim, function(sim_m, m) {
    if (is.null(sim_m) || (length(sim_m) == 1 && is.na(sim_m))) {
      return(NULL)
    }
    sim_m %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(model = as.integer(m))
  })

  if (nrow(sim_df) == 0) {
    stop("No non-NA simulations found for this deer.")
  }

  if (!is.null(model)) {
    sim_df <- sim_df %>% dplyr::filter(model == !!model)
    if (nrow(sim_df) == 0) {
      stop(sprintf("No simulations for model %d", model))
    }
  }

  # Group id so each (model, nsim) draws as its own line
  sim_df <- sim_df %>%
    dplyr::mutate(path_id = paste(model, nsim, sep = "_"))

  # Plot
  p <- ggplot() +
    geom_path(
      data = obs,
      aes(x = x_, y = y_, group = burst_),
      colour = "orange",
      linewidth = 0.7
    ) +
    geom_path(
      data = sim_df,
      aes(x = x_, y = y_, group = path_id),
      colour = "black",
      alpha = sim_alpha,
      linewidth = 0.3
    ) +
    coord_equal() +
    labs(
      title = sprintf(
        "Deer row %d (id: %s) - observed vs simulated (test)",
        row_no,
        if ("id" %in% names(deer_row)) as.character(deer_row$id) else ""
      ),
    ) +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.title = element_blank())

  if (is.null(model)) {
    p <- p + facet_wrap(~model)
  }

  p
}


#' Normalize an stp_test-shaped path to (x, y, burst_, path_id, frame).
#' Assigns sequential frame indices so that each burst starts `burst_gap`
#' frames after the previous burst ended, creating a visual pause during
#' which the meteorite tail can fade out before the next burst head appears.
.normalize_stp <- function(path, id, burst_gap = 5) {
  df <- tibble::tibble(
    x = path$x1_,
    y = path$y1_,
    burst_ = path$burst_
  )

  # Preserve original order while grouping rows within their burst
  df$.order <- seq_len(nrow(df))
  burst_ids <- unique(df$burst_)

  out <- purrr::imap_dfr(burst_ids, function(b, i) {
    sub <- df[df$burst_ == b, ]
    sub$burst_idx <- i
    sub
  })
  out <- out[order(out$.order), ]
  out$.order <- NULL

  # Assign frames: continuous within a burst, `burst_gap` pause between bursts
  rle_b <- rle(out$burst_idx)
  starts <- c(
    1,
    cumsum(rle_b$lengths[-length(rle_b$lengths)]) +
      1 +
      seq_len(length(rle_b$lengths) - 1) * burst_gap
  )
  frame_vec <- integer(nrow(out))
  pos <- 1
  for (k in seq_along(rle_b$lengths)) {
    n_k <- rle_b$lengths[k]
    frame_vec[pos:(pos + n_k - 1)] <- starts[k]:(starts[k] + n_k - 1)
    pos <- pos + n_k
  }
  out$frame <- frame_vec

  out$path_id <- id
  out
}


#' Animate one or two deer paths (stp_test-shaped) with a meteorite-tail
#' moving head. Single panel, colored by path_id.
#'
#' @param path   stp_test-like tibble (x1_, y1_, burst_, ...)
#' @param path2  Optional second path of the same shape
#' @param labels Legend labels for path / path2
#' @param colors Colors for path / path2
#' @param file   Output file (e.g. "plots/deer_42.mp4"). ".mp4" -> av,
#'               ".gif" -> gifski. If NULL, returns animation without saving.
#' @param fps    Frames per second (default 20)
#' @param duration Optional total duration (s); overrides fps-per-step
#' @param wake_length Tail length as fraction of total frames (default 0.05)
#' @param point_size Size of moving head
#' @param width,height,res Output dimensions (pixels) and resolution
#' @param burst_gap Empty frames inserted between bursts for tail fade-out
#' @param step_duration Seconds of animation per data step. Determines clip
#'                      length as n_steps * step_duration. Ignored if
#'                      `duration` is set.
animate_deer_path <- function(
  path,
  path2 = NULL,
  path3 = NULL,
  labels = c("observed", "simulated", "simulated2"),
  colors = c("orange", "blue", "green"),
  file = NULL,
  fps = 30,
  duration = NULL,
  step_duration = 0.15,
  wake_length = 0.05,
  point_size = 3,
  width = 800,
  height = 800,
  res = 150,
  burst_gap = 1
) {
  if (!requireNamespace("gganimate", quietly = TRUE)) {
    stop("Please install gganimate.")
  }

  # Workaround: amt::simulate_path returns n+1 rows (includes start point),
  # so simulated bursts are one step longer than observed bursts. Drop the
  # last row of each burst in sim paths if that pattern is detected.
  trim_sim <- function(sim_path, obs_path) {
    obs_counts <- table(obs_path$burst_)
    sim_counts <- table(sim_path$burst_)
    shared <- intersect(names(obs_counts), names(sim_counts))
    if (
      length(shared) > 0 &&
        all(sim_counts[shared] == obs_counts[shared] + 1)
    ) {
      sim_path <- sim_path %>%
        dplyr::group_by(burst_) %>%
        dplyr::slice(-dplyr::n()) %>%
        dplyr::ungroup()
    }
    sim_path
  }

  # Normalize inputs ---------------------------------------------------------
  df1 <- .normalize_stp(path, labels[1], burst_gap = burst_gap)
  dfs <- list(df1)
  path_levels <- labels[1]

  if (!is.null(path2)) {
    path2 <- trim_sim(path2, path)
    dfs[[length(dfs) + 1]] <- .normalize_stp(
      path2,
      labels[2],
      burst_gap = burst_gap
    )
    path_levels <- c(path_levels, labels[2])
  }
  if (!is.null(path3)) {
    path3 <- trim_sim(path3, path)
    dfs[[length(dfs) + 1]] <- .normalize_stp(
      path3,
      labels[3],
      burst_gap = burst_gap
    )
    path_levels <- c(path_levels, labels[3])
  }

  # Truncate all to the shared number of frames so they end together
  shared_n <- min(vapply(dfs, function(d) max(d$frame), numeric(1)))
  dfs <- lapply(dfs, function(d) d[d$frame <= shared_n, ])
  df <- dplyr::bind_rows(dfs)

  df$path_id <- factor(df$path_id, levels = path_levels)
  # Unique group per (path, burst) so lines/reveal break between bursts
  df$group_id <- interaction(df$path_id, df$burst_, drop = TRUE)

  n_steps <- max(df$frame)

  # Build animation ----------------------------------------------------------
  # Attach color directly (avoid scale_color_manual which can conflict with
  # shadow_wake's color interpolation on some gganimate versions).
  color_map <- stats::setNames(colors[seq_along(path_levels)], path_levels)
  df$color <- color_map[as.character(df$path_id)]

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = x, y = y, group = group_id, colour = color)
  ) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_colour_identity() +
    ggplot2::coord_equal() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none") +
    gganimate::transition_time(frame) +
    gganimate::shadow_wake(wake_length = wake_length)

  # Render -------------------------------------------------------------------
  anim_args <- list(
    plot = p,
    fps = fps,
    width = width,
    height = height,
    res = res,
    device = "png"
  )
  if (is.null(duration)) {
    duration <- n_steps * step_duration
  }
  anim_args$nframes <- max(round(duration * fps), 30)

  if (!is.null(file)) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "mp4") {
      if (!requireNamespace("av", quietly = TRUE)) {
        stop("Install 'av' for mp4 output.")
      }
      anim_args$renderer <- gganimate::av_renderer()
    } else if (ext == "gif") {
      if (!requireNamespace("gifski", quietly = TRUE)) {
        stop("Install 'gifski' for gif output.")
      }
      anim_args$renderer <- gganimate::gifski_renderer()
    } else {
      stop("file must end in .mp4 or .gif")
    }
  }

  anim <- do.call(gganimate::animate, anim_args)

  if (!is.null(file)) {
    dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
    gganimate::anim_save(file, animation = anim)
  }

  invisible(anim)
}

