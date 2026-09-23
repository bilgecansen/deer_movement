# xgboost-specific helpers -----------------------------------------------------
#
# Exploratory variable importance for step-selection data, fit with
# gradient-boosted trees on the same likelihood the amt and GAM paths use.
# The output is a ranking of how much each candidate variable matters
# relative to the others, to decide what is worth carrying into the modelling
# paths. It is not a filter: the gates in filter_models_gam.R /
# filter_models_amt.R do that job.
#
# Three pieces make the trees fit a step-selection model rather than a plain
# classifier:
#
#   * The loss is the conditional logit: within each stratum (one observed
#     step plus its random points) the scores go through a softmax, and the
#     likelihood is the probability given to the observed endpoint. Written
#     out, this is the same likelihood as amt's clogit and mgcv's cox.ph
#     form, so a linear booster on this loss reproduces amt's coefficients.
#     Strata never enter as a feature; they only set the normalisation.
#
#   * Variables are split across several boosters that are grown in
#     alternation, each holding the others' current prediction fixed. Their
#     sum is one additive model. Movement keeps its own booster so it is
#     always available, while habitat variables compete for one randomly
#     chosen column per tree.
#
#   * A booster can be "switched": its output is multiplied by an indicator
#     built at fit time from a column's NA pattern and anchored so it is 0
#     at a score of 0. That is how the FAMD scores enter — they exist only on
#     forest cells. Without it, xgboost's own missing-value branch turns the
#     NA pattern into a forest indicator and takes the forest effect away
#     from the landcover term.
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
#'   * "deer" — the average drop per deer, the headline number. It reads
#'     against the pipeline's gate 3, where a model must beat the null by 3
#'     log units for a deer, so a variable worth 2.3 per deer is worth about
#'     two thirds of that bar.
#'   * "steps100" — the drop per 100 observed steps. Deer differ in track
#'     length between seasons (winter deer contribute about four times as
#'     many steps each as autumn deer), so this is the one to use when
#'     comparing season-years.
#'
#' @param x Raw importance (a total over the pool)
#' @param n_deer,n_steps Deer and observed steps behind that total
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
#' Call it once for a fixed table (prep_pool_xgb.R), or once per replicate
#' with a different seed to redraw the random points
#' (replicate_importance_xgb.R).
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

#' Booster specifications for the alternating fit
#'
#' Movement keeps every column at every split. Habitat and FAMD trees are
#' single-variable (one interaction group per column), and with
#' `one_col_per_tree` each tree is offered one randomly chosen column, so a
#' strong variable cannot crowd the others out of every tree.
#'
#' @param move_vars Movement columns, always available
#' @param hab_vars Habitat columns
#' @param famd_vars Columns for the switched booster (NA outside their
#'   domain); empty for none
#' @param switch_col Column whose NA pattern defines the switch; defaults to
#'   the first entry of `famd_vars`
#' @param move_groups Interaction groups within the movement booster, by
#'   column name; the default mirrors the amt movement block, where step
#'   length interacts with time of day and the turning angle does not
#' @param categorical Names of categorical columns
#' @param learning_rate,max_depth,min_child_weight,reg_lambda Boosting
#'   settings shared by all three boosters
#' @param one_col_per_tree Sample a single column per tree in the habitat and
#'   FAMD boosters
#' @param nthread Threads per booster
#' @return List of booster specs for fit_xgb_boosters()
make_xgb_specs <- function(
  move_vars,
  hab_vars,
  famd_vars = character(0),
  switch_col = NULL,
  move_groups = list(c("sl_", "tod_day"), "cos_ta"),
  categorical = character(0),
  learning_rate = 0.05,
  max_depth = 2L,
  min_child_weight = 1,
  reg_lambda = 1,
  one_col_per_tree = TRUE,
  nthread = 1L
) {
  params_for <- function(constraints, colsample) {
    xgboost::xgb.params(
      booster = "gbtree",
      learning_rate = learning_rate,
      max_depth = max_depth,
      min_child_weight = min_child_weight,
      reg_lambda = reg_lambda,
      nthread = nthread,
      interaction_constraints = constraints,
      colsample_bytree = colsample
    )
  }
  # interaction_constraints take 0-based column positions.
  index_groups <- function(feats, groups) {
    lapply(groups, function(g) match(intersect(g, feats), feats) - 1)
  }
  singles <- function(feats) as.list(seq_along(feats) - 1)
  one_col <- function(feats) {
    if (one_col_per_tree) 1 / length(feats) else 1
  }

  specs <- list(list(
    name = "move",
    feats = move_vars,
    switched = FALSE,
    switch_col = NA_character_,
    categorical = categorical,
    params = params_for(index_groups(move_vars, move_groups), 1)
  ))
  specs[[length(specs) + 1]] <- list(
    name = "hab",
    feats = hab_vars,
    switched = FALSE,
    switch_col = NA_character_,
    categorical = categorical,
    params = params_for(singles(hab_vars), one_col(hab_vars))
  )
  if (length(famd_vars)) {
    specs[[length(specs) + 1]] <- list(
      name = "famd",
      feats = famd_vars,
      switched = TRUE,
      switch_col = if (is.null(switch_col)) famd_vars[1] else switch_col,
      categorical = categorical,
      params = params_for(singles(famd_vars), one_col(famd_vars))
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
#' Each round draws a bag of whole steps and learns only from those, so every
#' step is left out of roughly (1 - bag_frac) of the rounds. The trees from
#' the rounds that left a step out give that step an out-of-bag score, which
#' xgb_perm_importance() uses. Steps are bagged whole: sampling points inside
#' a step would break the within-stratum gradients that keep a step-constant
#' column unusable on its own. Set bag_frac = 1 to disable it.
#'
#' @param d Step data, sorted by stratum with the observed step first
#' @param specs Booster specs from make_xgb_specs()
#' @param n_rounds Trees per booster
#' @param bag_frac Share of steps each round learns from
#' @param verbose_every Print progress every this many rounds (0 = silent)
#' @return List with `boosters` (raw serialised), `n_trees`, the final
#'   in-sample conditional log-likelihood and the `bag` matrix (rounds x
#'   steps, NULL when bag_frac = 1)
fit_xgb_boosters <- function(d, specs, n_rounds, bag_frac = 0.632,
                             verbose_every = 0) {
  y <- d$case_
  sizes <- xgb_strata_sizes(d)
  strat_of_row <- rep(seq_along(sizes), sizes)
  bagging <- bag_frac < 1
  bag <- if (bagging) {
    matrix(FALSE, n_rounds, length(sizes))
  } else {
    NULL
  }
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
      bag[i, ] <- in_bag
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
    loglik = xgb_cond_loglik(Reduce(`+`, contrib), y, sizes),
    bag = bag
  )
}

#' Columns each tree of a booster splits on
#'
#' @param bst Fitted booster
#' @return List with one character vector per tree; empty for a leaf-only
#'   tree
xgb_tree_columns <- function(bst) {
  js <- jsonlite::fromJSON(
    paste(xgboost::xgb.dump(bst, dump_format = "json"), collapse = ""),
    simplifyVector = FALSE
  )
  lapply(js, function(tree) {
    cols <- character(0)
    walk <- function(node) {
      if (!is.null(node$split)) {
        cols <- c(cols, node$split)
        for (kid in node$children) {
          cols <- c(cols, walk(kid))
        }
      }
      cols
    }
    unique(walk(tree))
  })
}

#' Contribution of a single tree
#'
#' @param bst Fitted booster
#' @param spec Its spec
#' @param X Feature matrix to score
#' @param on Switch vector for the booster
#' @param i Tree index, 1-based
#' @return Contribution per row
xgb_tree_part <- function(bst, spec, X, on, i) {
  dm <- xgb_matrix(X, spec$categorical)
  if (spec$switched) {
    xgb_switched_part(bst, dm, xgb_zero_matrix(spec), on, i, i)
  } else {
    xgb_raw(bst, dm, i, i)
  }
}

#' Permutation importance, one value per variable
#'
#' A variable's values are permuted among the rows of the same stratum, which
#' breaks the link between the value and which point was used while leaving
#' each step's set of available values untouched. A switched booster's
#' variables are permuted only among the rows where the switch is on, so the
#' score's effect is measured separately from the on/off contrast. The score
#' is the resulting drop in conditional log-likelihood — the same units as
#' delta_logp — averaged over `n_perm` permutations.
#'
#' Two versions are returned. `in_sample` scores every step with the whole
#' model, so it includes whatever noise the trees fitted. `out_of_bag` scores
#' each step with only the trees from the rounds that left it out, which is
#' the version to rank by: a variable the trees merely fitted noise on scores
#' about zero or below. Out-of-bag needs the `bag` matrix from
#' fit_xgb_boosters(); without it only `in_sample` is returned.
#'
#' Because only a third or so of the trees score any given step, the
#' out-of-bag model is under-fitted and its numbers are smaller throughout.
#' Read them as a ranking within one fit, not as a measure of prediction.
#'
#' The scored boosters must grow single-variable trees (make_xgb_specs() does
#' this for the habitat and FAMD boosters), so that a variable's trees can be
#' rescored on their own; this is checked.
#'
#' @param boosters Raw serialised boosters from fit_xgb_boosters()
#' @param d Step data the model was fit on
#' @param specs Booster specs
#' @param bag Bag matrix from fit_xgb_boosters(), or NULL
#' @param vars Variables to score; defaults to every habitat and switched
#'   column (movement is not permuted)
#' @param n_perm Permutations to average over
#' @return data.frame(variable, booster, in_sample[, out_of_bag]), most
#'   important first
xgb_perm_importance <- function(boosters, d, specs, bag = NULL, vars = NULL,
                                n_perm = 3) {
  bst <- lapply(boosters, xgboost::xgb.load.raw)
  n_trees <- vapply(bst, xgboost::xgb.get.num.boosted.rounds, numeric(1))
  sizes <- xgb_strata_sizes(d)
  strat_of_row <- rep(seq_along(sizes), sizes)
  on <- lapply(specs, xgb_switch, d = d)
  X <- lapply(specs, function(s) as.matrix(d[, s$feats]))

  if (is.null(vars)) {
    vars <- unlist(lapply(specs[-1], `[[`, "feats"))
  }
  owner <- vapply(vars, function(v) {
    which(vapply(specs, function(s) v %in% s$feats, logical(1)))[1]
  }, numeric(1))

  tree_cols <- vector("list", length(specs))
  for (b in unique(owner)) {
    tree_cols[[b]] <- xgb_tree_columns(bst[[b]])
    if (any(lengths(tree_cols[[b]]) > 1)) {
      stop(sprintf(
        paste("Booster '%s' has trees using more than one column;",
              "importance needs single-variable trees."),
        specs[[b]]$name
      ))
    }
  }

  base_in <- Reduce(`+`, xgb_contributions(bst, d, specs, n_trees))
  ll_in <- xgb_cond_loglik(base_in, d$case_, sizes)
  # Out-of-bag score: for every row, the trees of the rounds that left its
  # step out of the bag.
  base_oob <- NULL
  if (!is.null(bag)) {
    base_oob <- rep(0, nrow(d))
    for (i in seq_len(nrow(bag))) {
      keep <- !bag[i, strat_of_row]
      for (b in seq_along(specs)) {
        if (i <= n_trees[b]) {
          part <- xgb_tree_part(bst[[b]], specs[[b]], X[[b]], on[[b]], i)
          base_oob[keep] <- base_oob[keep] + part[keep]
        }
      }
    }
    ll_oob <- xgb_cond_loglik(base_oob, d$case_, sizes)
  }

  drop_for <- function(v, b) {
    trees <- which(vapply(tree_cols[[b]],
                          function(cc) identical(cc, v), logical(1)))
    if (!length(trees)) {
      return(c(in_sample = 0, out_of_bag = 0))
    }
    rows <- which(on[[b]])
    rowMeans(replicate(n_perm, {
      Xp <- X[[b]]
      shuffled <- rows[order(d$stratum[rows], stats::runif(length(rows)))]
      Xp[rows, v] <- Xp[shuffled, v]
      delta_in <- rep(0, nrow(d))
      delta_oob <- rep(0, nrow(d))
      for (i in trees) {
        change <- xgb_tree_part(bst[[b]], specs[[b]], Xp, on[[b]], i) -
          xgb_tree_part(bst[[b]], specs[[b]], X[[b]], on[[b]], i)
        delta_in <- delta_in + change
        if (!is.null(bag)) {
          keep <- !bag[i, strat_of_row]
          delta_oob[keep] <- delta_oob[keep] + change[keep]
        }
      }
      c(
        in_sample = ll_in - xgb_cond_loglik(base_in + delta_in, d$case_,
                                            sizes),
        out_of_bag = if (is.null(bag)) {
          NA_real_
        } else {
          ll_oob - xgb_cond_loglik(base_oob + delta_oob, d$case_, sizes)
        }
      )
    }))
  }

  scores <- vapply(seq_along(vars),
                   function(j) drop_for(vars[j], owner[j]), numeric(2))
  out <- data.frame(
    variable = vars,
    booster = vapply(specs[owner], `[[`, character(1), "name"),
    in_sample = scores["in_sample", ],
    row.names = NULL
  )
  if (!is.null(bag)) {
    out$out_of_bag <- scores["out_of_bag", ]
  }
  rank_by <- if (is.null(bag)) out$in_sample else out$out_of_bag
  out[order(-rank_by), ]
}
