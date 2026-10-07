# Simulation helpers for the xgboost path ------------------------------------
#
# The simulation keeps a season's real steps, endpoints and columns, and
# replaces only which endpoint of each step was used, drawn from a true
# model written down in scripts/xgb_sim/sim_truth_xgb.R. These are the two
# pieces the simulation scripts share.
#
# Part of the helper library; scripts source scripts/helper_functions.R,
# which sources every file in this folder.

#' Standardised columns the true habitat functions are written in
#'
#' Each numeric column is centred and scaled over every row of the season
#' table; the FAMD axis over the rows that have it, and 0 where it is NA
#' (outside forest, where it has no effect).
#'
#' @param d The real season table
#' @return List of standardised columns, landcover as its integer code, and
#'   `famd_on`, TRUE where the FAMD axis exists
xgb_sim_columns <- function(d) {
  std <- function(x) {
    (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
  }
  famd_on <- !is.na(d$famd3_end)
  list(
    elevation = std(d$elevation_end),
    forest_edge = std(d$forest_edge_end),
    ndvi = std(d$ndvi_end),
    landcover = d$landcover,
    northness = std(d$northness_end),
    eastness = std(d$eastness_end),
    famd3 = ifelse(famd_on, std(d$famd3_end), 0),
    famd_on = famd_on
  )
}

#' The true habitat functions, before their strengths are applied
#'
#' Elevation is a smooth hump, or in scenario D a sharp band; forest edge a
#' saturating decline; NDVI, northness and the FAMD axis straight lines;
#' landcover a value per class (forest preferred, grassland mildly,
#' developed avoided); and the two interactions the products of standardised
#' elevation with forest edge and with NDVI. Eastness has no function: its
#' true effect is zero.
#'
#' @param Z Columns from xgb_sim_columns()
#' @param scenario "A", "B", "C" or "D"
#' @return Named list of term values, one per row
xgb_sim_terms <- function(Z, scenario) {
  lc_effect <- stats::setNames(rep(0, length(LANDCOVER_LEVELS)),
                               LANDCOVER_LEVELS)
  lc_effect[c("forest", "grassland", "developed")] <- c(1, 0.5, -1)
  list(
    elevation = if (scenario == "D") {
      as.numeric(Z$elevation > -1 & Z$elevation < 0.4)
    } else {
      exp(-(Z$elevation + 0.3)^2 / (2 * 0.8^2))
    },
    forest_edge = -tanh(Z$forest_edge),
    ndvi = Z$ndvi,
    landcover = unname(lc_effect[Z$landcover + 1]),
    northness = Z$northness,
    famd3 = Z$famd3,
    elev_x_edge = Z$elevation * Z$forest_edge,
    elev_x_ndvi = Z$elevation * Z$ndvi
  )
}

#' Draw the used endpoint of every step from a true model
#'
#' Within each step the true scores become probabilities (the same
#' within-step softmax the model fits), and one endpoint is drawn from them.
#'
#' @param eta True score per row, rows sorted by stratum
#' @param sizes Stratum sizes from xgb_strata_sizes()
#' @param seed Seed for the draw
#' @return Integer vector, one per step: the row (in the table's order) of
#'   the endpoint drawn as used
xgb_sim_choose <- function(eta, sizes, seed) {
  p <- xgb_softmax(eta, sizes)
  sor <- rep(seq_along(sizes), sizes)
  set.seed(seed)
  u <- stats::runif(length(sizes))
  # Cumulative probability within each step; the drawn endpoint is the
  # first whose cumulative probability reaches that step's uniform draw.
  cum <- stats::ave(p, sor, FUN = cumsum)
  # A step's last endpoint closes its total at exactly 1, so rounding can
  # never leave a step with nothing drawn
  cum[cumsum(sizes)] <- 1
  hit <- cum >= u[sor]
  first_hit <- !duplicated(sor[hit])
  out <- which(hit)[first_hit]
  stopifnot(length(out) == length(sizes))
  out
}

#' A season table with simulated choices
#'
#' The real table with case_ replaced by the simulated choice, and each
#' step's rows reordered so the used endpoint comes first, as every table
#' in the xgb path is laid out.
#'
#' @param pool_path The real season table (data/xgb/pooled_start_<s>.rds)
#' @param choice_path An xgb_sim_choose() result saved with saveRDS()
#' @return The table, ready for fit_model_xgb.R
xgb_sim_pool <- function(pool_path, choice_path) {
  d <- readRDS(pool_path)
  used <- readRDS(choice_path)
  d$case_ <- 0L
  d$case_[used] <- 1L
  d <- d[order(d$stratum, -d$case_), ]
  sizes <- xgb_strata_sizes(d)
  first <- cumsum(c(1, utils::head(sizes, -1)))
  stopifnot(all(d$case_[first] == 1), sum(d$case_) == length(sizes))
  d
}
