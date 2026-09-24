# xgboost-specific helpers -----------------------------------------------------
#
# Exploratory variable importance for step-selection data, fit with
# gradient-boosted trees on the same likelihood the amt and GAM paths use.
# The output is a ranking of how much each candidate variable matters
# relative to the others, to decide what is worth carrying into the modelling
# paths. It is not a filter: the gates in filter_models_gam.R /
# filter_models_amt.R do that job.
#
# Four pieces make the trees fit a step-selection model rather than a plain
# classifier:
#
#   * The loss is the conditional logit: within each stratum (one observed
#     step plus its random points) the scores go through a softmax, and the
#     likelihood is the probability given to the observed endpoint. Written
#     out, this is the same likelihood as amt's clogit and mgcv's cox.ph
#     form, so a linear booster on this loss reproduces amt's coefficients.
#     Strata never enter as a feature; they only set the normalisation.
#
#   * Variables are split across four boosters grown in alternation, each
#     holding the others' current prediction fixed. Their sum is one
#     additive model. Movement and distance to the home-range centre are the
#     nuisance block, mirroring the GAM null, so an importance is what a
#     variable adds on top of that null. They get their own boosters and a
#     full budget because they are controlled for rather than compared; the
#     environmental candidates and the FAMD axes get equal budgets because
#     they are compared with each other.
#
#   * A booster can be "switched": its output is multiplied by an indicator
#     built at fit time from a column's NA pattern and anchored so it is 0
#     at a score of 0. That is how the FAMD scores enter — they exist only on
#     forest cells. Without it, xgboost's own missing-value branch turns the
#     NA pattern into a forest indicator and takes the forest effect away
#     from the landcover term.
#
#   * Importance is the drop in held-out conditional log-likelihood over
#     folds of whole steps. Trees may combine variables (max_depth 2, so any
#     root-to-leaf path is pairwise), which costs the decomposition: a tree
#     splitting on landcover then NDVI is damaged by permuting either, so
#     shared structure is charged to both and the values do not add up to
#     anything. Read them as a ranking.
#
# Part of the helper library split out of scripts/helper_functions.R, which
# now sources every file in this folder. Scripts keep sourcing that one
# aggregator, so nothing here needs to be sourced directly.

#' Stratum sizes of a step-selection table
#'
#' Rows must be sorted so each stratum is contiguous, which is what
#' prep_pool_xgb.R writes.
#'
#' @param d Step data with a `stratum` column
#' @return Integer vector of stratum sizes, in row order
xgb_strata_sizes <- function(d) {
  rle(d$stratum)$lengths
}

#' Softmax within each stratum
#'
#' Scores are shifted by each stratum's first row before exponentiating,
#' which keeps exp() in range without changing the result.
#'
#' @param eta Score (linear predictor) per row
#' @param sizes Stratum sizes from xgb_strata_sizes()
#' @return Probability per row; each stratum sums to 1
xgb_softmax <- function(eta, sizes) {
  first <- cumsum(c(1, utils::head(sizes, -1)))
  id <- rep(seq_along(sizes), sizes)
  ex <- exp(eta - rep(eta[first], sizes))
  ex / rep(rowsum(ex, id, reorder = FALSE)[, 1], sizes)
}

#' Conditional log-likelihood of the observed steps
#'
#' @param eta Score per row
#' @param case_ Used/available indicator (1 = observed step)
#' @param sizes Stratum sizes from xgb_strata_sizes()
#' @return Sum over strata of log P(observed endpoint)
xgb_cond_loglik <- function(eta, case_, sizes) {
  p <- xgb_softmax(eta, sizes)
  sum(log(p[case_ == 1]))
}

#' DMatrix with the categorical columns declared
#'
#' @param X Numeric matrix; categorical columns hold integer codes 0..k-1
#' @param categorical Names of the categorical columns
#' @param label Used/available indicator, for training matrices
#' @param group Stratum sizes, for training matrices
#' @return xgb.DMatrix
xgb_matrix <- function(X, categorical = character(0), label = NULL,
                       group = NULL) {
  xgboost::xgb.DMatrix(
    X,
    label = label,
    group = group,
    feature_types = ifelse(colnames(X) %in% categorical, "c", "q")
  )
}

#' Rescale an importance value for reporting
#'
#' The raw drop is a total over every observed step in the pool, which runs
#' to hundreds or thousands and means little on its own. Two rescalings:
#'
#'   * "deer" — the average drop per deer-year, the headline number. It
#'     reads against the pipeline's gate 3, where a model must beat the null
#'     by 3 log units for a deer, so a variable worth 2.3 per deer is worth
#'     about two thirds of that bar.
#'   * "steps100" — the drop per 100 observed steps. Deer differ in track
#'     length between seasons (winter deer contribute about four times as
#'     many steps each as autumn deer), so this is the one to reach for when
#'     comparing seasons.
#'
#' Neither rescaling makes seasons strictly comparable. Permutation
#' importance measures how hard the model leans on a variable, and a model
#' fitted to fewer steps over-fits more and leans harder on everything; pf
#' sits about three points of held-out score further from its in-sample fit
#' than fa or nb does, and its values are inflated to match. Ranks within a
#' season are safe; magnitudes across seasons are not.
#'
#' @param x Raw importance (a total over the pool)
#' @param n_deer,n_steps Deer-years and observed steps behind that total
#' @param per "deer", "steps100" or "total"
#' @return The rescaled value
xgb_scale_value <- function(x, n_deer, n_steps, per = "deer") {
  switch(
    per,
    deer = x / n_deer,
    steps100 = 100 * x / n_steps,
    total = x,
    stop("per must be 'deer', 'steps100' or 'total'")
  )
}

#' Axis label for a rescaled importance
#' @param per "deer", "steps100" or "total"
#' @return Label text
xgb_scale_label <- function(per = "deer") {
  switch(
    per,
    deer = "Drop in total log score per deer",
    steps100 = "Drop in total log score per 100 steps",
    total = "Drop in total log score",
    stop("per must be 'deer', 'steps100' or 'total'")
  )
}

#' Build the pooled step table for one season and year
#'
#' Generates `n_pts` uniform-disc random points per observed step for every
#' deer in `keys`, extracts the covariates, and stacks the deer into one
#' table sorted by stratum with the observed step first — the row order the
#' objective and the permutations rely on.
#'
#' Landcover is written as 0-based integer codes in LANDCOVER_LEVELS order,
#' for xgboost's categorical splits. FAMD columns keep their NA off forest;
#' a row missing any other listed variable is dropped, and with it any step
#' left without its observed point.
#'
#' `deer` is an index within this call, not an animal id; the animal is the
#' prefix of `key`. xgb_season_pool() renumbers when it stacks years.
#'
#' @param keys Deer keys (<id>_<season>_<year>) to include
#' @param tracks The filtered deer table, with a `key` column
#' @param rasters List of the loaded rasters: landcover, water, ndvi,
#'   landfire, topo
#' @param n_pts Random points per observed step
#' @param hab_vars Habitat columns to keep, including "landcover"
#' @param famd_vars FAMD columns to keep (NA off forest)
#' @param move_vars Movement columns to keep
#' @return The pooled step table
build_xgb_pool <- function(keys, tracks, rasters, n_pts, hab_vars,
                           famd_vars,
                           move_vars = c("sl_", "tod_day", "cos_ta")) {
  pooled <- purrr::map_dfr(seq_along(keys), function(i) {
    row <- tracks[tracks$key == keys[i], ]
    stopifnot(nrow(row) == 1)
    one <- make_random_pt_extraction(
      data = row, n_pts = n_pts, water = rasters$water, model = "nonp",
      stp_col = "stp", output_col = "rand"
    )
    one <- extract_step_variables(
      data = one, env = rasters$landcover, ndvi = rasters$ndvi,
      landfire = rasters$landfire, topo = rasters$topo,
      random_col = "rand", output_col = "var"
    )
    v <- as.data.frame(one$var[[1]])
    out <- tibble::tibble(
      key = keys[i],
      deer = i,
      step_id_ = v$step_id_,
      case_ = as.integer(v$case_),
      sl_ = v$sl_,
      tod_day = as.numeric(v$tod_start_ == "day"),
      cos_ta = cos(v$ta_),
      landcover = as.integer(factor(v$wiscland_end,
                                    levels = LANDCOVER_LEVELS)) - 1L
    )
    for (nm in setdiff(c(hab_vars, famd_vars), "landcover")) {
      out[[nm]] <- v[[nm]]
    }
    out
  })

  pooled |>
    dplyr::filter(dplyr::if_all(dplyr::all_of(c(move_vars, hab_vars)),
                                ~ !is.na(.))) |>
    dplyr::group_by(deer, step_id_) |>
    dplyr::filter(sum(case_) == 1, dplyr::n() >= 2) |>
    dplyr::ungroup() |>
    dplyr::arrange(deer, step_id_, dplyr::desc(case_)) |>
    dplyr::mutate(
      stratum = as.integer(factor(paste(deer, step_id_),
                                  unique(paste(deer, step_id_))))
    )
}

#' Stack one season's years into a single table
#'
#' Strata are renumbered over the whole season from `key`, which carries the
#' animal, the season and the year. `deer` is only an index within each
#' year's file, so stacking on it alone would merge two years of the same
#' animal into one stratum; it is renumbered here to run across the season.
#'
#' Nothing else needs adding. Year is constant within a stratum, so it
#' cancels in the softmax: there is no year term to fit and no random effect
#' to worry about. The annual landcover and NDVI rasters were already read
#' per deer-year when the per-year files were built, so stacking rows stacks
#' values each taken against its own year's raster.
#'
#' Strata need not be the same size. A few lose an available point that fell
#' outside a raster, and the softmax normalises over whatever each stratum
#' holds; what is required is that every stratum keeps its observed step and
#' at least one alternative.
#'
#' @param season Season code ("fa", "nb", "pf")
#' @param dir Folder of per-year files from prep_pool_xgb.R
#' @return One table for the season, with `year`, `animal` and a renumbered
#'   `deer` (one value per deer-year)
xgb_season_pool <- function(season, dir = "data/xgb") {
  files <- list.files(dir, sprintf("^pooled_%s_[0-9]{4}[.]rds$", season),
                      full.names = TRUE)
  if (!length(files)) {
    stop(sprintf("no pooled files for season '%s' in %s", season, dir))
  }
  parts <- lapply(sort(files), function(f) {
    p <- readRDS(f)$pooled
    p$year <- as.integer(sub(".*_([0-9]{4})[.]rds$", "\\1", f))
    p
  })
  d <- dplyr::bind_rows(parts)
  d$animal <- sub("_[a-z]+_[0-9]{4}$", "", d$key)
  d$deer <- as.integer(factor(d$key, levels = unique(d$key)))
  d$stratum <- as.integer(factor(paste(d$key, d$step_id_),
                                 levels = unique(paste(d$key,
                                                       d$step_id_))))

  # The fit needs each stratum contiguous with its observed step first,
  # which is how the per-year files are written and how binding preserves
  # them. Checked rather than assumed.
  sizes <- rle(d$stratum)$lengths
  stopifnot(
    length(sizes) == dplyr::n_distinct(d$stratum),
    min(sizes) >= 2,
    all(d$case_[cumsum(c(1, utils::head(sizes, -1)))] == 1),
    sum(d$case_) == length(sizes)
  )
  d
}

#' Booster specifications for the alternating fit
#'
#' Four boosters, in the order they take their turn each round:
#'   move  movement, always available, never ranked
#'   hr    distance to the home-range centre, on its own, never ranked
#'   hab   the environmental candidates
#'   famd  the ordination axes, switched off outside forest
#'
#' move and hr are the nuisance block and mirror the GAM null
#' (movement + s(HR_center_end)). Putting HR on its own is what lets it be
#' fitted properly: sharing the habitat booster it took about a seventh of
#' the trees and stayed under-fitted, and no amount of extra rounds fixed
#' that without feeding the same rounds to variables that only fit noise.
#'
#' The ranked boosters get every column at every split and no interaction
#' constraint, so a tree may combine variables; `max_depth` 2 keeps any
#' root-to-leaf path pairwise. Column sampling was dropped once HR left the
#' candidate pool: it existed to stop HR taking every tree, the ranking is
#' unchanged from one sampled column per split up to all of them, and the
#' held-out score moves by about a tenth of a percent across that range.
#'
#' @param move_vars Movement columns
#' @param hr_var Home-range-centre column
#' @param hab_vars Environmental candidate columns
#' @param famd_vars Columns for the switched booster (NA outside their
#'   domain); empty for none
#' @param switch_col Column whose NA pattern defines the switch; defaults to
#'   the first entry of `famd_vars`
#' @param move_groups Interaction groups within the movement booster, by
#'   column name; the default mirrors the amt movement block, where step
#'   length interacts with time of day and the turning angle does not
#' @param categorical Names of categorical columns
#' @param learning_rate,max_depth,min_child_weight,reg_lambda Boosting
#'   settings shared by all boosters
#' @param nthread Threads per booster
#' @return List of booster specs for fit_xgb_boosters()
make_xgb_specs <- function(
  move_vars,
  hr_var,
  hab_vars,
  famd_vars = character(0),
  switch_col = NULL,
  move_groups = list(c("sl_", "tod_day"), "cos_ta"),
  categorical = character(0),
  learning_rate = 0.05,
  max_depth = 2L,
  min_child_weight = 1,
  reg_lambda = 1,
  nthread = 1L
) {
  params_for <- function(constraints = NULL) {
    xgboost::xgb.params(
      booster = "gbtree",
      learning_rate = learning_rate,
      max_depth = max_depth,
      min_child_weight = min_child_weight,
      reg_lambda = reg_lambda,
      nthread = nthread,
      interaction_constraints = constraints
    )
  }
  # interaction_constraints take 0-based column positions.
  index_groups <- function(feats, groups) {
    lapply(groups, function(g) match(intersect(g, feats), feats) - 1)
  }

  specs <- list(
    list(
      name = "move",
      feats = move_vars,
      switched = FALSE,
      switch_col = NA_character_,
      categorical = categorical,
      params = params_for(index_groups(move_vars, move_groups))
    ),
    list(
      name = "hr",
      feats = hr_var,
      switched = FALSE,
      switch_col = NA_character_,
      categorical = categorical,
      params = params_for()
    ),
    list(
      name = "hab",
      feats = hab_vars,
      switched = FALSE,
      switch_col = NA_character_,
      categorical = categorical,
      params = params_for()
    )
  )
  if (length(famd_vars)) {
    specs[[length(specs) + 1]] <- list(
      name = "famd",
      feats = famd_vars,
      switched = TRUE,
      switch_col = if (is.null(switch_col)) famd_vars[1] else switch_col,
      categorical = categorical,
      params = params_for()
    )
  }
  specs
}

#' Raw margin of one booster over a range of its trees
#'
#' @param bst Fitted booster
#' @param dm DMatrix to predict on
#' @param from,to Tree range, 1-based and inclusive
#' @return Margin per row
xgb_raw <- function(bst, dm, from, to) {
  stats::predict(bst, dm, outputmargin = TRUE,
                 iterationrange = c(from, to))
}

#' Contribution of a switched booster
#'
#' The margin is anchored by subtracting its value at a score of 0 and then
#' zeroed wherever the switch is off, so the term is exactly 0 outside the
#' variables' domain and at the centre of the score space. The level that
#' separates inside from outside therefore stays with whatever term carries
#' it (landcover, in the pooled models).
#'
#' @param bst Fitted booster
#' @param dm DMatrix to predict on
#' @param dm0 One-row DMatrix of zeros, from xgb_zero_matrix()
#' @param on Logical switch per row
#' @param from,to Tree range, 1-based and inclusive
#' @return Contribution per row
xgb_switched_part <- function(bst, dm, dm0, on, from, to) {
  ifelse(on, xgb_raw(bst, dm, from, to) - xgb_raw(bst, dm0, from, to)[1], 0)
}

#' One-row matrix of zeros, the anchor point of a switched booster
#' @param spec A booster spec from make_xgb_specs()
#' @return xgb.DMatrix with one row
xgb_zero_matrix <- function(spec) {
  X <- matrix(0, 1, length(spec$feats),
              dimnames = list(NULL, spec$feats))
  xgb_matrix(X, spec$categorical)
}

#' Switch vector for a spec (TRUE where the booster's variables exist)
#' @param spec A booster spec
#' @param d Step data
#' @return Logical vector, all TRUE for unswitched boosters
xgb_switch <- function(spec, d) {
  if (!spec$switched) {
    return(rep(TRUE, nrow(d)))
  }
  !is.na(d[[spec$switch_col]])
}

#' Contribution of every booster on `d`
#'
#' @param boosters List of fitted boosters
#' @param d Step data
#' @param specs Booster specs
#' @param n_trees Trees to use per booster; defaults to all of them
#' @return List of contributions, one vector per booster
xgb_contributions <- function(boosters, d, specs, n_trees = NULL) {
  if (is.null(n_trees)) {
    n_trees <- vapply(boosters, xgboost::xgb.get.num.boosted.rounds,
                      numeric(1))
  }
  lapply(seq_along(specs), function(b) {
    X <- as.matrix(d[, specs[[b]]$feats])
    dm <- xgb_matrix(X, specs[[b]]$categorical)
    if (specs[[b]]$switched) {
      xgb_switched_part(boosters[[b]], dm, xgb_zero_matrix(specs[[b]]),
                        xgb_switch(specs[[b]], d), 1, n_trees[b])
    } else {
      xgb_raw(boosters[[b]], dm, 1, n_trees[b])
    }
  })
}

#' Grow one tree on an existing booster
#'
#' Updating in place keeps the cost per round flat. xgb.iter.update() is
#' xgboost-internal, so this falls back to continuing with xgb.train() if a
#' future version drops it; the fit is the same, only slower.
#'
#' @param bst Booster to extend, or NULL to start one
#' @param dm Training DMatrix
#' @param params Booster parameters
#' @param n_done Trees already grown (the iteration number xgboost expects)
#' @param obj Objective returning list(grad, hess)
#' @return The booster, extended by one tree
xgb_grow_one <- function(bst, dm, params, n_done, obj) {
  if (is.null(bst)) {
    return(xgboost::xgb.train(params, dm, nrounds = 1, objective = obj,
                              verbose = 0))
  }
  update_one <- get0("xgb.iter.update", envir = asNamespace("xgboost"))
  if (is.function(update_one)) {
    update_one(bst, dm, n_done, obj)
    return(bst)
  }
  xgboost::xgb.train(params, dm, nrounds = 1, objective = obj,
                     xgb_model = bst, verbose = 0)
}

#' Fit the alternating boosters
#'
#' Each round grows one tree per booster, in turn, with the other boosters'
#' current prediction held fixed. Boosters are updated in place and only the
#' newest tree is predicted each round, so run time grows linearly with the
#' number of rounds rather than quadratically.
#'
#' Every objective ignores xgboost's own `preds` argument and uses the score
#' this function tracks. That is what lets a switched booster be trained on
#' the anchored, switched-off score rather than on its raw output.
#'
#' Each round draws a bag of whole steps and learns only from those. This is
#' Friedman's stochastic gradient boosting, here purely as a regulariser —
#' the honest score comes from xgb_cv_importance()'s folds, not from the
#' out-of-bag steps. Steps are bagged whole: sampling points inside a step
#' would break the within-stratum gradients, which sum to zero over a
#' stratum and must keep doing so. Set bag_frac = 1 to disable it.
#'
#' @param d Step data, sorted by stratum with the observed step first
#' @param specs Booster specs from make_xgb_specs()
#' @param n_rounds Trees per booster
#' @param bag_frac Share of steps each round learns from
#' @param verbose_every Print progress every this many rounds (0 = silent)
#' @return List with `boosters` (raw serialised), `n_trees` and the final
#'   in-sample conditional log-likelihood
fit_xgb_boosters <- function(d, specs, n_rounds, bag_frac = 0.632,
                             verbose_every = 0) {
  y <- d$case_
  sizes <- xgb_strata_sizes(d)
  strat_of_row <- rep(seq_along(sizes), sizes)
  bagging <- bag_frac < 1
  on <- lapply(specs, xgb_switch, d = d)
  X <- lapply(specs, function(s) as.matrix(d[, s$feats]))
  dm <- lapply(seq_along(specs), function(b) {
    xgb_matrix(X[[b]], specs[[b]]$categorical, label = y, group = sizes)
  })
  dm_pred <- lapply(seq_along(specs), function(b) {
    xgb_matrix(X[[b]], specs[[b]]$categorical)
  })
  dm0 <- lapply(specs, xgb_zero_matrix)
  contrib <- lapply(specs, function(s) rep(0, nrow(d)))
  bst <- vector("list", length(specs))
  n_trees <- integer(length(specs))
  t0 <- Sys.time()

  for (i in seq_len(n_rounds)) {
    in_bag <- rep(TRUE, length(sizes))
    if (bagging) {
      in_bag <- rep(FALSE, length(sizes))
      in_bag[sample.int(length(sizes), round(bag_frac * length(sizes)))] <-
        TRUE
    }
    row_in_bag <- in_bag[strat_of_row]
    for (b in seq_along(specs)) {
      eta_now <- Reduce(`+`, contrib)
      w <- row_in_bag * if (specs[[b]]$switched) on[[b]] else 1
      # Gradient and Hessian of the conditional logit, zeroed where a
      # switched booster's variables do not exist, so those rows shape no
      # split while still counting as available points in their stratum.
      obj <- function(preds, dtrain) {
        p <- xgb_softmax(eta_now, sizes)
        list(grad = w * (p - y), hess = w * p * (1 - p))
      }
      bst[[b]] <- xgb_grow_one(bst[[b]], dm[[b]], specs[[b]]$params,
                               n_trees[b], obj)
      n_trees[b] <- n_trees[b] + 1L
      nt <- n_trees[b]
      contrib[[b]] <- contrib[[b]] + if (specs[[b]]$switched) {
        xgb_switched_part(bst[[b]], dm_pred[[b]], dm0[[b]], on[[b]], nt, nt)
      } else {
        xgb_raw(bst[[b]], dm_pred[[b]], nt, nt)
      }
    }
    if (verbose_every > 0 && i %% verbose_every == 0) {
      cat(sprintf(
        "  round %d / %d  %.0f s  logLik %.1f\n", i, n_rounds,
        as.numeric(difftime(Sys.time(), t0, units = "secs")),
        xgb_cond_loglik(Reduce(`+`, contrib), y, sizes)
      ))
    }
  }

  # The running per-tree totals must match a full prediction: up to a
  # constant for an ordinary booster (each tree's margin carries base_score,
  # which cancels within a stratum) and exactly for a switched one.
  for (b in seq_along(specs)) {
    full <- if (specs[[b]]$switched) {
      xgb_switched_part(bst[[b]], dm_pred[[b]], dm0[[b]], on[[b]], 1,
                        n_trees[b])
    } else {
      xgb_raw(bst[[b]], dm_pred[[b]], 1, n_trees[b])
    }
    gap <- contrib[[b]] - full
    if (diff(range(gap)) > 1e-4) {
      stop(sprintf(
        paste("Booster '%s': the running per-tree total drifted from a",
              "full prediction."),
        specs[[b]]$name
      ))
    }
    if (specs[[b]]$switched && any(contrib[[b]][!on[[b]]] != 0)) {
      stop(sprintf("Booster '%s' is not 0 where its switch is off",
                   specs[[b]]$name))
    }
  }

  list(
    boosters = lapply(bst, xgboost::xgb.save.raw),
    n_trees = stats::setNames(n_trees, vapply(specs, `[[`, character(1),
                                              "name")),
    loglik = xgb_cond_loglik(Reduce(`+`, contrib), y, sizes)
  )
}

#' Permutation importance from random folds of whole steps
#'
#' Steps are dealt at random into `n_folds` folds. Each fold is scored by a
#' model fitted on the other folds, so a held-out step never reached any
#' round of the model scoring it, and that model carries its full complement
#' of trees. Every step is held out exactly once, so summing the folds'
#' drops gives a total over all of them.
#'
#' A variable is permuted among the rows of the same stratum, which breaks
#' the link between its value and which point was used while leaving each
#' step's set of available values untouched. A switched booster's variables
#' are permuted only where the switch is on, so the score's effect is
#' measured separately from the on/off contrast.
#'
#' Fold at the stratum level, never the row: a step split across folds would
#' leave the softmax denominator inconsistent with the availability set that
#' defines it.
#'
#' Two things to keep in mind when reading the output. A held-out step's
#' neighbours in time are in the training set, and the conditional
#' likelihood removes whatever is constant within a step, so this leaks less
#' than it would for an unconditional model but not nothing. And with trees
#' free to combine variables the drops no longer decompose: shared structure
#' is charged to every variable involved.
#'
#' Give the data a column of noise and it lands at zero here, which is where
#' to read "carries nothing" from. Zero itself is not that mark for every
#' scheme — scoring by the out-of-bag steps instead puts a null column
#' several units *below* zero, the further below the more trees it was
#' given, which is why that route was dropped.
#'
#' @param d Step data the model is fit on
#' @param specs Booster specs
#' @param n_rounds Trees per booster
#' @param bag_frac Passed to fit_xgb_boosters()
#' @param n_folds Folds of whole steps
#' @param n_perm Permutations averaged per variable per fold
#' @param seed Seed for the fold draw
#' @param verbose Print a line per fold
#' @return data.frame(variable, booster, cv), most important first, with the
#'   held-out log-likelihood in the "ll_cv" attribute
xgb_cv_importance <- function(d, specs, n_rounds, bag_frac = 0.632,
                              n_folds = 5, n_perm = 3, seed = 1,
                              verbose = FALSE) {
  strata <- unique(d$stratum)
  set.seed(seed)
  fold_of_stratum <- sample(rep_len(seq_len(n_folds), length(strata)))
  fold_of_row <- fold_of_stratum[match(d$stratum, strata)]

  ranked <- which(vapply(specs, function(s) s$name %in% c("hab", "famd"),
                         logical(1)))
  vars <- unlist(lapply(specs[ranked], `[[`, "feats"))
  owner <- vapply(vars, function(v) {
    which(vapply(specs, function(s) v %in% s$feats, logical(1)))[1]
  }, numeric(1))

  total <- stats::setNames(numeric(length(vars)), vars)
  ll_total <- 0

  for (f in seq_len(n_folds)) {
    tr <- d[fold_of_row != f, ]
    te <- d[fold_of_row == f, ]
    fit <- fit_xgb_boosters(tr, specs, n_rounds = n_rounds,
                            bag_frac = bag_frac, verbose_every = 0)
    bst <- lapply(fit$boosters, xgboost::xgb.load.raw)
    n_trees <- vapply(bst, xgboost::xgb.get.num.boosted.rounds, numeric(1))

    sizes_te <- xgb_strata_sizes(te)
    on_te <- lapply(specs, xgb_switch, d = te)
    X_te <- lapply(specs, function(s) as.matrix(te[, s$feats]))
    contrib <- xgb_contributions(bst, te, specs, n_trees)
    base <- Reduce(`+`, contrib)
    ll <- xgb_cond_loglik(base, te$case_, sizes_te)
    ll_total <- ll_total + ll
    if (verbose) {
      cat(sprintf("  fold %d / %d  held-out logLik %.1f\n", f, n_folds, ll))
    }

    # Permuting a variable only touches the booster that owns it, so only
    # that booster's contribution is recomputed — and on the held-out fold
    # every tree is in play, so that is one prediction, not one per tree.
    for (j in seq_along(vars)) {
      b <- owner[j]
      rows <- which(on_te[[b]])
      total[j] <- total[j] + mean(replicate(n_perm, {
        Xp <- X_te[[b]]
        sh <- rows[order(te$stratum[rows], stats::runif(length(rows)))]
        Xp[rows, vars[j]] <- Xp[sh, vars[j]]
        dm <- xgb_matrix(Xp, specs[[b]]$categorical)
        newc <- if (specs[[b]]$switched) {
          xgb_switched_part(bst[[b]], dm, xgb_zero_matrix(specs[[b]]),
                            on_te[[b]], 1, n_trees[b])
        } else {
          xgb_raw(bst[[b]], dm, 1, n_trees[b])
        }
        ll - xgb_cond_loglik(base - contrib[[b]] + newc, te$case_,
                             sizes_te)
      }))
    }
  }

  out <- data.frame(
    variable = vars,
    booster = vapply(specs[owner], `[[`, character(1), "name"),
    cv = as.numeric(total),
    row.names = NULL
  )
  attr(out, "ll_cv") <- ll_total
  out[order(-out$cv), ]
}

#' Every split, and every parent-child split pair, of one booster
#'
#' A tree using two variables is not yet an interaction. The unit is a
#' root-to-leaf path: the root splits on A, a child splits on B, so the tree
#' is saying how B matters depends on which side of A's threshold you are.
#' At max_depth 2 every interaction is one such parent-child pair.
#'
#' This reads the fitted trees and nothing else. That also sets the limit:
#' a greedy learner always finds a second split, there is no "no
#' interaction" option once a tree has depth, and the gains in the dump are
#' in-sample.
#'
#' @param bst Fitted booster
#' @return List of two data frames: `splits` (one row per internal node) and
#'   `pairs` (one row per parent-child pair of internal nodes)
xgb_tree_structure <- function(bst) {
  # with_stats carries the gain and cover of each split; without it the
  # dump has structure only.
  js <- jsonlite::fromJSON(
    paste(xgboost::xgb.dump(bst, dump_format = "json", with_stats = TRUE),
          collapse = ""),
    simplifyVector = FALSE
  )
  splits <- list()
  pairs <- list()

  # A categorical split carries a set of levels rather than a cut point, so
  # its threshold is recorded as NA.
  threshold_of <- function(node) {
    sc <- node$split_condition
    if (is.null(sc) || length(sc) != 1 || !is.numeric(sc)) {
      return(NA_real_)
    }
    as.numeric(sc)
  }
  is_internal <- function(node) !is.null(node$split)

  walk <- function(node, tree_id) {
    if (!is_internal(node)) {
      return(invisible(NULL))
    }
    splits[[length(splits) + 1]] <<- data.frame(
      tree = tree_id,
      depth = if (is.null(node$depth)) NA_integer_ else node$depth,
      variable = node$split,
      threshold = threshold_of(node),
      gain = if (is.null(node$gain)) NA_real_ else node$gain
    )
    for (kid in node$children) {
      if (is_internal(kid)) {
        pairs[[length(pairs) + 1]] <<- data.frame(
          tree = tree_id,
          parent = node$split,
          side = if (identical(kid$nodeid, node$yes)) "yes" else "no",
          child = kid$split,
          child_gain = if (is.null(kid$gain)) NA_real_ else kid$gain
        )
      }
      walk(kid, tree_id)
    }
    invisible(NULL)
  }

  for (i in seq_along(js)) {
    walk(js[[i]], i)
  }
  list(
    splits = if (length(splits)) do.call(rbind, splits) else NULL,
    pairs = if (length(pairs)) do.call(rbind, pairs) else NULL
  )
}

#' Pair frequency against what independence would give
#'
#' Expected counts come from the parent and child marginals of the pair
#' table itself, so the ratio asks whether the model chose *this* child
#' under *this* parent more often than its overall taste for each would
#' explain. A ratio near 1 is what an additive model looks like once it has
#' been forced to grow depth-2 trees.
#'
#' Read it against a noise column rather than against 1. On the pooled fits
#' the largest ratio anywhere is noise pairing with itself, and a variable
#' splitting on itself — refinement of one curve, not an interaction — is
#' the most enriched pattern throughout. What separates the real pairs is
#' gain share, not ratio.
#'
#' @param pairs The `pairs` frame from xgb_tree_structure()
#' @return Ordered pairs with observed, expected, ratio and gain share
xgb_pair_lift <- function(pairs) {
  n <- nrow(pairs)
  p_parent <- table(pairs$parent) / n
  p_child <- table(pairs$child) / n
  tab <- pairs |>
    dplyr::group_by(parent, child) |>
    dplyr::summarise(observed = dplyr::n(),
                     gain = sum(child_gain, na.rm = TRUE),
                     .groups = "drop")
  tab$expected <- n * as.numeric(p_parent[tab$parent]) *
    as.numeric(p_child[tab$child])
  tab$ratio <- tab$observed / tab$expected
  tab$gain_share <- tab$gain / sum(tab$gain)
  tab[order(-tab$ratio), ]
}
