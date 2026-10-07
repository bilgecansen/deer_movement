#' @description
#' MANUAL WALKTHROUGH of the xgboost path: one season, one round of boosting,
#' one tree opened up, one held-out fold, one shuffled variable, one split
#' scored on held-out steps, and the decision rules applied to the stored
#' results. Straight-line code meant to be stepped through line by line, and
#' a quick check that the code works as intended. Companion to
#' walkthrough_gam.R and walkthrough_amt.R.
#'
#' The helpers in scripts/helpers/helpers_xgb.R are INLINED here. Library
#' calls (xgboost, survival) stay as calls.
#'
#' There are no loops in sections 1-10. Every loop in the xgb path is shown
#' for one iteration: one round of the round x booster loop in
#' fit_xgb_boosters(), one fold of xgb_cv_importance(), one variable and one
#' permutation of its importance loop, one tree of the replay in
#' xgb_heldout_split_gain(). Once section 5 has built one round by hand,
#' section 8 lets fit_xgb_boosters() run a few dozen rounds for the fold
#' model. Section 11 audits every stored result file, and that is the one
#' place this file iterates.
#'
#' WHY THE ASSERTIONS MATTER. Two kinds of stopifnot() run here:
#'   * Against a production helper on the same inputs. An inlined copy is a
#'     second implementation and can drift; if it does, you would be
#'     reviewing code that is not what runs. These fail loudly the day a
#'     helper changes.
#'   * Against a design decision: yearlings are gone, the home-range centre
#'     reaches the habitat trees only in rsf_hr, the threshold is the same in
#'     every script, the stored results came from the current pools and the
#'     current settings. These catch a misunderstanding rather than a typo.
#' Do not delete either kind.
#'
#' Sections:
#'   0  configuration, read from the production scripts
#'   1  the season table                          (prep_start_xgb.R)
#'   2  the loss: conditional logit               (xgb_softmax, ..._loglik)
#'   3  which way each score should move: p - y
#'   4  columns into boosters                     (make_xgb_specs)
#'   5  one round of alternating boosting         (fit_xgb_boosters, 1 round)
#'   6  inside one tree: leaf values and gain by hand
#'   7  why a start-of-step column is never a root split
#'   8  one held-out fold                         (xgb_cv_importance, 1 fold)
#'   9  permutation importance, one variable      (xgb_cv_importance, 1 job)
#'  10  held-out gain of one split                (xgb_heldout_split_gain)
#'  11  from numbers to decisions, and an audit of results/xgb/
#'
#' Usage: open it and step through. Running it top to bottom also works, in
#' under a minute, and stops at the first assertion that fails.

suppressPackageStartupMessages({
  library(xgboost)
  library(survival)
  library(tidyverse)
})

# The helpers are loaded ONLY so the assertions can compare against them,
# and for the few library-sized pieces named where they are used.
source("scripts/helper_functions.R")


# 0 ---- Configuration ---------------------------------------------------------
# The settings that define the analysis are READ from the production scripts
# rather than retyped, so this file cannot quietly walk through a different
# model from the one that was fit.
config_value <- function(script, name) {
  pattern <- sprintf("^%s <- ", name)
  line <- grep(pattern, readLines(script), value = TRUE)
  stopifnot(length(line) == 1)
  eval(parse(text = sub(pattern, "", line)))
}
FIT_SCRIPT <- "scripts/xgb/fit_model_xgb.R"
N_ROUNDS <- config_value(FIT_SCRIPT, "N_ROUNDS")
LEARNING_RATE <- config_value(FIT_SCRIPT, "LEARNING_RATE")
MAX_DEPTH <- config_value(FIT_SCRIPT, "MAX_DEPTH")
BAG_FRAC <- config_value(FIT_SCRIPT, "BAG_FRAC")
N_FOLDS <- config_value(FIT_SCRIPT, "N_FOLDS")
SEED <- config_value(FIT_SCRIPT, "SEED")
# The L2 penalty on leaf values, make_xgb_specs()'s default
LAMBDA <- formals(make_xgb_specs)$reg_lambda

# Choices for the walkthrough itself
SEASON <- "pf" # the smallest season, so every section runs in seconds
CONFIG <- "full" # the one model type with all five boosters
FOLD <- 1L # the held-out fold walked through
WALK_ROUNDS <- 50L # rounds for the fold model in section 8
N_THREAD <- 4L
# Used to show why a stratum-constant column cannot be a root split
START_COL <- "forest_edge_start"

cat(sprintf(
  "production settings: %d rounds per booster, learning rate %g, depth %d,",
  N_ROUNDS,
  LEARNING_RATE,
  MAX_DEPTH
))
cat(sprintf(" lambda %g, bag %g, %d folds\n", LAMBDA, BAG_FRAC, N_FOLDS))


# 1 ---- The season table (prep_start_xgb.R) ----------------------------------
# One row per candidate endpoint. A stratum is one observed step: the point
# the deer reached (case_ = 1) plus the random points it could have reached
# instead (case_ = 0). Every year of the season is stacked into one table.
d <- readRDS(sprintf("data/xgb/pooled_start_%s.rds", SEASON))

# The noise columns, drawn exactly as fit_model_xgb.R draws them: one value
# per row, so they vary within a stratum like any end-point column.
set.seed(SEED)
d$shadow_gauss <- stats::rnorm(nrow(d))
set.seed(SEED + 1L)
d$shadow_famd <- stats::rnorm(nrow(d))

# Strata are contiguous and the observed step comes first. Everything below
# relies on that layout.
sizes <- rle(d$stratum)$lengths
n_strata <- length(sizes)
strat_of_row <- rep(seq_len(n_strata), sizes)
first_row <- cumsum(c(1, utils::head(sizes, -1)))
y <- d$case_
stopifnot(
  n_strata == dplyr::n_distinct(d$stratum),
  all(y[first_row] == 1),
  sum(y) == n_strata
)
cat(sprintf(
  "%s: %d deer-years, %d animals, %d strata (observed steps), %d rows\n",
  SEASON,
  dplyr::n_distinct(d$key),
  dplyr::n_distinct(d$animal),
  n_strata,
  nrow(d)
))
cat(sprintf("candidates per stratum: %d to %d\n", min(sizes), max(sizes)))

# DESIGN: yearlings are gone. Deer collared at 20 months are dropped from
# every track that starts before they turn two.
ages <- deer_track_age(readRDS("library/SW_filtered_deer.RData"))
young <- ages$key[!is.na(ages$age_months) & ages$age_months < 24]
stopifnot(!any(unique(d$key) %in% young))

# Two kinds of column, told apart by how they vary across one step's
# candidates:
#   end-point columns  (elevation_end, landcover, ...) are measured where each
#                      candidate ends, so they differ from candidate to
#                      candidate.
#   start columns      (forest_edge_start, ...) are measured where the step
#                      began. All the candidates of a step begin at the same
#                      point, so a start column has one value per step,
#                      repeated on every candidate.
# Each noise column copies the kind it stands in for: shadow_gauss (habitat)
# gets a new value on every row, like an end-point column; shadow_start
# (start block) gets one value per step, like a start column.
constant_within <- function(x) all(x == x[first_row][strat_of_row])
stopifnot(
  constant_within(d[[START_COL]]),
  constant_within(d$shadow_start),
  !constant_within(d$shadow_gauss)
)

# The FAMD axes describe forest stands, so they have values only on cells
# the stand data cover (almost all of them forest) and are NA everywhere
# else. The FAMD booster is switched on for rows that have values and off
# for rows that are NA (section 5). To tell the two apart it looks at a
# single column, famd1_end. That shortcut is right only if all five axes are
# NA on exactly the same rows, so check it: a row with famd1_end present but
# another axis missing would be switched on with a hole in its data.
famd_cols <- grep("^famd[0-9]+_end$", names(d), value = TRUE)
stopifnot(all(vapply(
  famd_cols,
  function(v) {
    identical(is.na(d[[v]]), is.na(d[[famd_cols[1]]]))
  },
  logical(1)
)))
cat(sprintf(
  "FAMD present on %.1f%% of rows\n",
  100 * mean(!is.na(d[[famd_cols[1]]]))
))


# 2 ---- The loss: conditional logit (xgb_softmax, xgb_cond_loglik) -----------
# Every candidate gets a score, eta. Within a stratum the scores go through a
# softmax: p = exp(eta) / sum(exp(eta)) over that stratum's candidates, so the
# candidates of one step share probability 1 between them. The log
# likelihood is the log probability of the candidate the deer actually
# chose, summed over steps. Points from different steps never compete.
#
# Subtracting each stratum's first score before exponentiating keeps exp()
# in range and changes nothing: a constant cancels in the ratio.
cond_softmax <- function(eta, sizes) {
  sor <- rep(seq_along(sizes), sizes)
  first <- cumsum(c(1, utils::head(sizes, -1)))
  ex <- exp(eta - eta[first][sor])
  ex / rowsum(ex, sor, reorder = FALSE)[sor, 1]
}
cond_loglik <- function(eta, y, sizes) {
  sum(log(cond_softmax(eta, sizes)[y == 1]))
}

set.seed(SEED)
eta_test <- stats::rnorm(nrow(d))
stopifnot(
  all.equal(cond_softmax(eta_test, sizes), xgb_softmax(eta_test, sizes)),
  all.equal(
    cond_loglik(eta_test, y, sizes),
    xgb_cond_loglik(eta_test, y, sizes)
  )
)

# The null: every score 0, every candidate equally likely, so each step
# contributes -log(number of candidates).
null_ll <- cond_loglik(rep(0, nrow(d)), y, sizes)
stopifnot(all.equal(null_ll, -sum(log(sizes))))
cat(sprintf(
  "null log-likelihood %.1f (%.3f per step)\n",
  null_ll,
  null_ll / n_strata
))

# The same likelihood amt maximises. survival::clogit (what amt::fit_clogit
# calls) evaluated at a fixed coefficient on step length, without
# iterating, gives the same number as our loss with eta = b * sl_.
d_sub <- d[d$stratum %in% unique(d$stratum)[1:300], ]
sizes_sub <- rle(d_sub$stratum)$lengths
b <- -0.002
cl <- suppressWarnings(survival::clogit(
  case_ ~ sl_ + strata(stratum),
  data = d_sub,
  init = b,
  control = survival::coxph.control(iter.max = 0)
))
stopifnot(all.equal(
  cl$loglik[1],
  cond_loglik(b * d_sub$sl_, d_sub$case_, sizes_sub)
))
cat(sprintf("clogit and our loss agree at b = %g: %.6f\n", b, cl$loglik[1]))


# 3 ---- Which way each score should move: p - y ------------------------------
# Boosting never works on the log-likelihood directly. Before each tree it
# asks one question of every candidate: should this score go up or down, and
# how hard? The answer is p - y, where p is the probability the model
# currently gives the candidate and y is 1 for the point the deer chose, 0
# otherwise.
#
# At the very start, before any tree, every score is 0 and each of a step's
# 101 candidates has p = 1/101:
#   the chosen point    p - y = 1/101 - 1 = -0.990   push its score up, hard
#   each random point   p - y = 1/101 - 0 = +0.010   push it down, gently
p0 <- 1 / sizes[strat_of_row]
g0 <- p0 - y
# The first step: its chosen point, then three of its random points
print(round(g0[1:4], 4))

# THE LOAD-BEARING FACT. The one push up exactly balances the hundred pushes
# down: -0.990 + 100 x 0.010 = 0. That holds for every step at every stage of
# fitting, because a step's probabilities always add to 1 and exactly one of
# its points was chosen. Section 7 shows what follows from it.
stopifnot(max(abs(rowsum(g0, strat_of_row))) < 1e-12)

# Alongside the push, xgboost uses p * (1 - p): how much room that
# probability still has to move. Section 6 shows how the two combine into
# the value a tree adds.
h0 <- p0 * (1 - p0)

# 4 ---- Columns into boosters (fit_model_xgb.R, make_xgb_specs) --------------
# The model is a sum of separate boosters, each with its own columns:
#   move      step length, time of day and turning angle      nuisance
#   modifier  movement again, with start-of-step columns       ranked
#   hr        distance to the home-range centre                nuisance
#   hab       the environmental candidates at the endpoint     ranked
#   famd      the forest ordination axes, switched off outside forest
# move + hr mirror the GAM null, so whatever the ranked blocks add is what a
# variable adds on top of movement and the home-range centre.
MOVE_VARS <- c("sl_", "tod_day", "cos_ta")
HR_VAR <- "HR_center_end"
FAMD_VARS <- c(famd_cols, "shadow_famd")
START_VARS <- c(
  intersect(c("wiscland_start", "ndvi_start", "forest_edge_start"), names(d)),
  "shadow_start"
)
ID_COLS <- c("key", "deer", "animal", "year", "step_id_", "case_", "stratum")
HAB_VARS <- setdiff(
  names(d),
  c(MOVE_VARS, HR_VAR, FAMD_VARS, START_VARS, ID_COLS)
)
CATEGORICAL <- c("landcover", "wiscland_start")

# DESIGN: the habitat candidates are end-point columns plus their noise column.
stopifnot(setequal(
  HAB_VARS,
  c(
    "landcover",
    "forest_edge_end",
    "elevation_end",
    "northness_end",
    "eastness_end",
    "shadow_gauss",
    intersect("ndvi_end", names(d))
  )
))
cat(sprintf("habitat: %s\n", paste(HAB_VARS, collapse = ", ")))

# The specs for CONFIG, written out (make_xgb_specs inlined). Every booster
# grows depth-2 trees with the same learning rate and penalty.
params_for <- function(constraints = NULL) {
  xgboost::xgb.params(
    booster = "gbtree",
    learning_rate = LEARNING_RATE,
    max_depth = MAX_DEPTH,
    min_child_weight = 1,
    reg_lambda = LAMBDA,
    nthread = N_THREAD,
    interaction_constraints = constraints
  )
}
block <- function(
  name,
  feats,
  params,
  switched = FALSE,
  switch_col = NA_character_
) {
  list(
    name = name,
    feats = feats,
    switched = switched,
    switch_col = switch_col,
    categorical = CATEGORICAL,
    params = params
  )
}
specs <- list(
  # Constraints are 0-based column positions: step length may combine with
  # time of day (0 with 1); the turning angle (2) stands alone, as in amt.
  block("move", MOVE_VARS, params_for(list(c(0, 1), 2))),
  block("modifier", c("sl_", "cos_ta", START_VARS), params_for()),
  block("hr", HR_VAR, params_for()),
  block("hab", HAB_VARS, params_for()),
  block(
    "famd",
    FAMD_VARS,
    params_for(),
    switched = TRUE,
    switch_col = FAMD_VARS[1]
  )
)
specs_of <- function(cfg) {
  make_xgb_specs(
    MOVE_VARS,
    HR_VAR,
    HAB_VARS,
    FAMD_VARS,
    config = cfg,
    start_vars = START_VARS,
    categorical = CATEGORICAL,
    learning_rate = LEARNING_RATE,
    max_depth = MAX_DEPTH,
    nthread = N_THREAD
  )
}
stopifnot(CONFIG == "full", all.equal(specs, specs_of(CONFIG)))
names(specs) <- vapply(specs, `[[`, "", "name")

# DESIGN: the six model types differ only in what the table below says.
# blocks_of() describes a list of booster specs in the terms the design was
# agreed in. Section 11 uses it again on the stored results. A type with no
# habitat block (null) has none of the habitat properties.
blocks_of <- function(s) {
  nm <- vapply(s, `[[`, "", "name")
  hab <- if ("hab" %in% nm) s[[which(nm == "hab")]] else NULL
  tibble::tibble(
    habitat = !is.null(hab),
    modifier = "modifier" %in% nm,
    hr_in_hab = !is.null(hab) && HR_VAR %in% hab$feats,
    one_variable_trees = !is.null(hab) &&
      !is.null(hab$params$interaction_constraints)
  )
}
describe <- function(cfg) {
  dplyr::bind_cols(tibble::tibble(type = cfg), blocks_of(specs_of(cfg)))
}
design <- dplyr::bind_rows(
  describe("null"),
  describe("main"),
  describe("rsf"),
  describe("rsf_hr"),
  describe("start"),
  describe("full")
)
print(as.data.frame(design), row.names = FALSE)
stopifnot(identical(
  design,
  tibble::tibble(
    type = c("null", "main", "rsf", "rsf_hr", "start", "full"),
    habitat = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE),
    modifier = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
    hr_in_hab = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
    one_variable_trees = c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE)
  )
))

# 5 ---- One round of alternating boosting (fit_xgb_boosters, one round) ------
# Each round grows ONE tree per booster, in turn. Before each tree the
# current score of every row (the sum of all boosters so far) goes through
# the softmax, giving fresh gradients, and the new tree is fit to those. So
# the boosters share one loss and take turns reducing it.
#
# Each round learns from a bag of whole steps (63.2%), a regulariser only.
# Rows outside the bag get weight 0. Steps are bagged whole: splitting a
# stratum would break the zero-sum gradients.
set.seed(SEED)
in_bag <- rep(FALSE, n_strata)
in_bag[sample.int(n_strata, round(BAG_FRAC * n_strata))] <- TRUE
row_in_bag <- in_bag[strat_of_row]

dm_of <- function(data, s, label = NULL, group = NULL) {
  xgboost::xgb.DMatrix(
    as.matrix(data[, s$feats]),
    label = label,
    group = group,
    feature_types = ifelse(s$feats %in% s$categorical, "c", "q")
  )
}
# Grow one tree on one booster. Left to itself, xgboost would work out the
# pushes (section 3) from its own trees only. But each booster here sees
# only some of the columns, and a candidate's probability depends on the
# score from ALL of them added together. So the pushes are worked out here,
# from that combined score (eta_now), and handed to xgboost through `obj`;
# xgboost's own version, `preds`, is never used. That is what makes five
# separate boosters one model.
#
# w is each row's weight: 0 for steps outside this round's bag, and for the
# FAMD booster also 0 outside forest. A row with weight 0 pushes nothing, so
# it has no say in where the tree splits.
grow_one_tree <- function(s, eta_now, w) {
  p <- cond_softmax(eta_now, sizes)
  # xgboost always calls obj(preds, dtrain); both are deliberately unused.
  obj <- function(preds, dtrain) {
    list(grad = w * (p - y), hess = w * p * (1 - p))
  }
  bst <- xgboost::xgb.train(
    s$params,
    dm_of(d, s, label = y, group = sizes),
    nrounds = 1,
    objective = obj,
    verbose = 0
  )
  list(bst = bst, p = p, w = w)
}
# A tree's contribution is its leaf value for each row (plus a constant that
# xgboost adds to every row, which cancels within a stratum).
#
# outputmargin = TRUE asks for that raw score, the "margin", before xgboost
# applies any transformation of its own. With one of its built-in
# objectives it would otherwise return, say, a logistic-transformed
# probability. Our transformation is the probability among the candidates
# of one step (section 2), which xgboost cannot do because it does not know
# what a step is, so we always take the raw score and transform it
# ourselves. Every predict() below asks for it for the same reason.
margin_of <- function(bst, data, s) {
  stats::predict(
    bst,
    dm_of(data, s),
    outputmargin = TRUE,
    iterationrange = c(1, 1)
  )
}

# Why the scores are ADDED, booster after booster. Each tree learns a
# correction to the score the trees before it have already set, not a score
# of its own. The move tree gives every candidate a score; the pushes that
# follow say what is still wrong after move; the modifier tree is fit to
# those pushes, so its leaf values are amounts to add on top of move's. The
# model's score is therefore move's part plus the modifier's correction.
#
# One step, three candidates, the deer chose A:
#                                          A     B     C     P(A)
#   start, every score 0                  0     0     0     0.33
#   after the move tree                   0.3   0     0     0.40
#   push for A: 0.40 - 1 = -0.60, so the modifier tree adds +0.1 to A
#   summed, c_move + c_mod                0.4   0     0     0.43
#   modifier alone, c_mod                 0.1   0     0     0.36
# The +0.1 means "a further 0.1 on top of move's 0.3", not "A deserves 0.1".
# Keeping only c_mod would leave A worse off than after move alone, applying
# a correction to scores it never saw. Nothing is counted twice: move's
# score enters the total once, and the modifier tree only used it as the
# starting point it was correcting.
#
# This is how a single booster already works: tree 2 is fit to what tree 1
# left, and the prediction is tree 1 + tree 2. Here the running total spans
# five boosters instead of one, which is also why the objective ignores
# xgboost's own `preds`: those hold only that booster's trees.
eta <- rep(0, nrow(d))
t_move <- grow_one_tree(specs$move, eta, row_in_bag)
c_move <- margin_of(t_move$bst, d, specs$move)

eta <- c_move
t_mod <- grow_one_tree(specs$modifier, eta, row_in_bag)
c_mod <- margin_of(t_mod$bst, d, specs$modifier)

eta <- c_move + c_mod
t_hr <- grow_one_tree(specs$hr, eta, row_in_bag)
c_hr <- margin_of(t_hr$bst, d, specs$hr)

eta <- c_move + c_mod + c_hr
t_hab <- grow_one_tree(specs$hab, eta, row_in_bag)
c_hab <- margin_of(t_hab$bst, d, specs$hab)

# The FAMD booster is SWITCHED. Outside forest its rows get weight 0, so they
# shape no split, and its contribution there is forced to 0. Inside forest
# the contribution is anchored by subtracting the prediction at all-zero
# axes, so the booster carries what the axes add WITHIN forest and the
# forest/non-forest contrast stays with landcover.
on_famd <- !is.na(d[[specs$famd$switch_col]])
eta <- c_move + c_mod + c_hr + c_hab
t_famd <- grow_one_tree(specs$famd, eta, row_in_bag * on_famd)
zero_row <- matrix(
  0,
  1,
  length(specs$famd$feats),
  dimnames = list(NULL, specs$famd$feats)
)
zero_dm <- xgboost::xgb.DMatrix(
  zero_row,
  feature_types = rep("q", length(specs$famd$feats))
)
c_famd <- ifelse(
  on_famd,
  margin_of(t_famd$bst, d, specs$famd) -
    stats::predict(t_famd$bst, zero_dm, outputmargin = TRUE)[1],
  0
)
eta_round1 <- c_move + c_mod + c_hr + c_hab + c_famd
ll_round1 <- cond_loglik(eta_round1, y, sizes)

# Two checks that the pushes behave as section 3 said, now at a real stage
# of fitting rather than the start:
#   * when the habitat tree was grown, after three trees had already moved
#     the scores, each step's pushes still summed to zero;
#   * the round made the fit better than the null. Had p - y pointed the
#     wrong way, every tree would have pushed scores the wrong way and the
#     log-likelihood would have fallen.
stopifnot(
  max(abs(rowsum(t_hab$p - y, strat_of_row))) < 1e-10,
  ll_round1 > null_ll
)

# Against production: the same round from fit_xgb_boosters().
set.seed(SEED)
fit1 <- fit_xgb_boosters(d, specs, n_rounds = 1, bag_frac = BAG_FRAC)
prod_c <- xgb_contributions(
  lapply(fit1$boosters, xgboost::xgb.load.raw),
  d,
  specs
)
stopifnot(
  all.equal(list(c_move, c_mod, c_hr, c_hab, c_famd), prod_c),
  all.equal(ll_round1, fit1$loglik)
)
cat(sprintf(
  "after one round of %d trees: log-likelihood %.1f (null %.1f)\n",
  length(specs),
  ll_round1,
  null_ll
))


# 6 ---- Inside one tree: leaf values and gain by hand ------------------------
# The habitat tree just grown, read from the model's own JSON. Nodes are
# numbered from 1 here; a leaf has no children.
js <- jsonlite::fromJSON(
  rawToChar(xgboost::xgb.save.raw(t_hab$bst, raw_format = "json")),
  simplifyVector = FALSE
)
tree_js <- js$learner$gradient_booster$model$trees[[1]]
feature_names <- unlist(js$learner$feature_names)
tr <- list(
  left = unlist(tree_js$left_children) + 1L,
  right = unlist(tree_js$right_children) + 1L,
  variable = feature_names[unlist(tree_js$split_indices) + 1L],
  weight = unlist(tree_js$base_weights),
  hess = unlist(tree_js$sum_hessian),
  gain = unlist(tree_js$loss_changes)
)
tr$variable[tr$left == 0L] <- NA_character_
stopifnot(all.equal(tr, xgb_tree_arrays(t_hab$bst)[[1]]))

kids <- c(tr$left[1], tr$right[1])
cat(sprintf(
  "habitat tree 1: root splits on %s; children split on %s\n",
  tr$variable[1],
  paste(
    ifelse(is.na(tr$variable[kids]), "(leaf)", tr$variable[kids]),
    collapse = " and "
  )
))

# Which leaf each row lands in
leaf <- as.integer(stats::predict(
  t_hab$bst,
  dm_of(d, specs$hab),
  predleaf = TRUE
)) +
  1L

# The gradients this tree was grown on: the softmax of the score before it,
# times each row's weight (0 outside the bag).
g_hab <- t_hab$w * (t_hab$p - y)
h_hab <- t_hab$w * t_hab$p * (1 - t_hab$p)

# A leaf's value is how far to move the score of every row in it: the
# "adjustment" those rows get. It comes from two sums over the leaf:
#   G  the total push, p - y. Negative means the scores should go up.
#   H  how fast that push fades as the scores move: raise a candidate's
#      score by one unit and its push changes by about p * (1 - p).
# Moving every score in the leaf by -G / H would make the leaf's pushes add
# up to about zero, using up the push. Two brakes are applied to that:
#   lambda         added to H. It shrinks the move for a leaf with little H,
#                  where only a few rows would decide a large move.
#   learning rate  only this share (5%) of the move is taken, so later trees
#                  can still correct it.
# value = -learning_rate * G / (H + lambda)
# (Textbooks call -G / H a Newton step, after Newton's method for finding
# where something reaches zero.)
#
# `rows` counts every row that lands in the leaf. Only the bagged ones carry
# a push, so only they shape G and H.
G_leaf <- rowsum(g_hab, leaf)
H_leaf <- rowsum(h_hab, leaf)
leaf_ids <- as.integer(rownames(G_leaf))
by_hand <- -LEARNING_RATE * G_leaf[, 1] / (H_leaf[, 1] + LAMBDA)
stopifnot(
  all.equal(unname(by_hand), tr$weight[leaf_ids], tolerance = 1e-4),
  all.equal(unname(H_leaf[, 1]), tr$hess[leaf_ids], tolerance = 1e-4)
)
print(
  data.frame(
    leaf = leaf_ids,
    rows = as.vector(table(leaf)),
    G = G_leaf[, 1],
    H = H_leaf[, 1],
    value = by_hand
  ),
  row.names = FALSE,
  digits = 4
)

# Each row's contribution from this tree IS its leaf's value, up to the
# constant xgboost adds to every row.
stopifnot(diff(range(c_hab - tr$weight[leaf])) < 1e-5)

# How the tree chose its splits. A split is scored by how much better the
# fit gets if its two sides get their own leaf value instead of one shared
# value:
#   gain = GL^2/(HL + lambda) + GR^2/(HR + lambda) - G^2/(H + lambda)
# where GL, HL are the sums of the pushes and of p * (1 - p) on the left,
# GR, HR the same on the right, and G, H the two sides together. xgboost
# tries every column and every cut point, and keeps the largest gain.
split_gain <- function(left, g, h) {
  gl <- sum(g[left])
  hl <- sum(h[left])
  gr <- sum(g[!left])
  hr <- sum(h[!left])
  gl^2 / (hl + LAMBDA) + gr^2 / (hr + LAMBDA) - (gl + gr)^2 / (hl + hr + LAMBDA)
}

# The root split's gain, rebuilt from the G and H of each side
below <- function(k) if (tr$left[k] == 0L) k else c(tr$left[k], tr$right[k])
goes_left <- leaf %in% below(kids[1])
root_gain <- split_gain(goes_left, g_hab, h_hab)
stopifnot(abs(root_gain - tr$gain[1]) < 1e-4 * root_gain)
cat(sprintf(
  "root gain %.3f by hand, %.3f reported by xgboost\n",
  root_gain,
  tr$gain[1]
))
# xgboost stores a leaf's value already multiplied by the learning rate and
# an internal node's weight without it. The difference matters later:
#   * Section 10 asks what a node would have added had the tree stopped there
#     instead of splitting it in two. That value is the learning rate times
#     the node's stored weight. The check below confirms it: it equals the
#     leaf value worked out by hand from every row under the node, as if
#     the node were a leaf.
#   * xgb_check_tree() in the helpers rebuilds every split's gain from the
#     stored weights. To get G back it has to undo the learning rate on
#     leaves but not on internal nodes.
k_inner <- kids[tr$left[kids] != 0L][1] # a child that splits again
under_inner <- leaf %in% below(k_inner)
stopifnot(
  !is.na(k_inner),
  all.equal(
    -LEARNING_RATE *
      sum(g_hab[under_inner]) /
      (sum(h_hab[under_inner]) + LAMBDA),
    LEARNING_RATE * tr$weight[k_inner],
    tolerance = 1e-4
  )
)


# 7 ---- Why a start-of-step column is never a root split ---------------------
# Section 6 scored a split by summing the pushes on each side. Apply the same
# score to the start of fitting (the pushes of section 3, every score 0) and
# it shows why the modifier booster's start-of-step columns can only ever
# appear BELOW a split on step length or turning angle.

# A start column is constant within a stratum, so any split on it sends
# WHOLE strata to each side. Each side is then a sum of whole steps, whose
# pushes sum to zero (section 3): GL = GR = 0 and the gain is exactly 0.
on_start <- d[[START_COL]] < stats::median(d[[START_COL]])
gain_start_root <- split_gain(on_start, g0, h0)
stopifnot(abs(sum(g0[on_start])) < 1e-8, abs(gain_start_root) < 1e-10)

# A column that varies within a stratum can separate a stratum's chosen
# point from its alternatives, so it earns gain.
short <- d$sl_ < stats::median(d$sl_)
gain_sl_root <- split_gain(short, g0, h0)
gain_noise_root <- split_gain(d$shadow_gauss < 0, g0, h0)

cat(sprintf(
  "root gain: %s %.2e | sl_ %.2f | per-row noise %.2f\n",
  START_COL,
  gain_start_root,
  gain_sl_root,
  gain_noise_root
))
# Noise earns gain too: gain is measured on the data the split was chosen
# on. Section 10 measures it on held-out steps instead.

# Below a split on step length each stratum is cut in two, and the halves no
# longer cancel. A stratum's short half sums to the probability the model
# gave its short candidates, minus 1 if the deer took a short step.
stopifnot(abs(sum(g0[short & on_start])) > 1)

# So a start column CAN now find gain, if steps that start in one place take
# short steps more often than steps that start elsewhere. That is the
# interaction between movement and where the step began. Whether it does
# depends on the column; the per-step noise column shows what none looks
# like. A split with no signal scores slightly BELOW zero, because lambda
# penalises the extra leaf, and xgboost never makes such a split.
numeric_start <- setdiff(START_VARS, "wiscland_start")
under_sl <- vapply(
  numeric_start,
  function(v) {
    split_gain(d[[v]][short] < stats::median(d[[v]]), g0[short], h0[short])
  },
  numeric(1)
)
cat("gain of a median split under sl_ < median:\n")
print(round(under_sl, 2))

# And in the model actually fit: the modifier tree grown in section 5 splits
# first on a movement column, whatever start columns it was offered.
mod_root <- xgb_tree_arrays(t_mod$bst)[[1]]$variable[1]
stopifnot(mod_root %in% c("sl_", "cos_ta"))
cat(sprintf("modifier tree 1 splits first on %s\n", mod_root))


# 8 ---- One held-out fold (xgb_cv_importance, one fold) ----------------------
# Whole steps are dealt at random into folds. Each fold is scored by a model
# fit on the other folds, so a held-out step never reached any tree that
# scores it. Summed over the folds, every step is held out exactly once.
set.seed(SEED)
strata <- unique(d$stratum)
fold_of_stratum <- sample(rep_len(seq_len(N_FOLDS), length(strata)))
fold_of_row <- fold_of_stratum[match(d$stratum, strata)]

# Folds are cut by stratum, never by row: a step split across folds would
# leave its softmax without the candidates that define it.
held <- rowsum(as.integer(fold_of_row == FOLD), strat_of_row)[, 1]
stopifnot(all(held == 0 | held == sizes))

tr_d <- d[fold_of_row != FOLD, ]
te <- d[fold_of_row == FOLD, ]
sizes_te <- rle(te$stratum)$lengths

# Section 5 built one round by hand; here fit_xgb_boosters() runs that
# round WALK_ROUNDS times (production runs N_ROUNDS).
fit_fold <- fit_xgb_boosters(
  tr_d,
  specs,
  n_rounds = WALK_ROUNDS,
  bag_frac = BAG_FRAC
)
bst <- lapply(fit_fold$boosters, xgboost::xgb.load.raw)
names(bst) <- names(specs)

# Each booster's contribution on the held-out steps, all trees at once
all_trees <- function(b, s) {
  stats::predict(b, dm_of(te, s), outputmargin = TRUE)
}
on_te <- !is.na(te[[specs$famd$switch_col]])
c_te <- list(
  move = all_trees(bst$move, specs$move),
  modifier = all_trees(bst$modifier, specs$modifier),
  hr = all_trees(bst$hr, specs$hr),
  hab = all_trees(bst$hab, specs$hab),
  famd = ifelse(
    on_te,
    all_trees(bst$famd, specs$famd) -
      stats::predict(bst$famd, zero_dm, outputmargin = TRUE)[1],
    0
  )
)
stopifnot(all.equal(unname(c_te), xgb_contributions(bst, te, specs)))

eta_te <- Reduce(`+`, c_te)
ll_te <- cond_loglik(eta_te, te$case_, sizes_te)
null_te <- -sum(log(sizes_te))
null_tr <- -sum(log(rle(tr_d$stratum)$lengths))
cat(sprintf(
  "fold %d, %d rounds: gain over null per step %.3f held out, %.3f in sample\n",
  FOLD,
  WALK_ROUNDS,
  (ll_te - null_te) / length(sizes_te),
  (fit_fold$loglik - null_tr) / sum(tr_d$case_)
))
# The in-sample gain is larger: part of what the trees learned is specific
# to the steps they were fit on. Only the held-out number counts.

# 9 ---- Permutation importance, one variable (xgb_cv_importance, one job) ----
# Shuffle one column among the candidates of each stratum and see how much
# the held-out log-likelihood drops. Each step keeps its own set of
# available values; only the link between a value and the choice is broken.
# Only the booster that owns the column is recomputed.
shuffle_within <- function(X, v, rows) {
  sh <- rows[order(te$stratum[rows], stats::runif(length(rows)))]
  X[rows, v] <- X[sh, v]
  X
}
predict_X <- function(b, s, X) {
  stats::predict(
    b,
    xgboost::xgb.DMatrix(
      X,
      feature_types = ifelse(s$feats %in% s$categorical, "c", "q")
    ),
    outputmargin = TRUE
  )
}
X_hab <- as.matrix(te[, specs$hab$feats])
all_rows <- seq_len(nrow(te))

set.seed(SEED)
Xp <- shuffle_within(X_hab, "elevation_end", all_rows)
# Each stratum still holds exactly the values it had, in a new order
key_of <- function(X) {
  X[order(te$stratum, X[, "elevation_end"]), "elevation_end"]
}
stopifnot(identical(key_of(X_hab), key_of(Xp)))
drop_elev <- ll_te -
  cond_loglik(
    eta_te - c_te$hab + predict_X(bst$hab, specs$hab, Xp),
    te$case_,
    sizes_te
  )
Xn <- shuffle_within(X_hab, "shadow_gauss", all_rows)
drop_noise <- ll_te -
  cond_loglik(
    eta_te - c_te$hab + predict_X(bst$hab, specs$hab, Xn),
    te$case_,
    sizes_te
  )
cat(sprintf(
  "drop when shuffled: elevation %.2f, noise %.2f (this fold)\n",
  drop_elev,
  drop_noise
))

# A start column is constant within a stratum, so shuffling it there changes
# nothing and scores exactly zero. Start columns are instead dealt between
# strata, a whole step at a time, and stay constant within each stratum.
X_mod <- as.matrix(te[, specs$modifier$feats])
stopifnot(identical(shuffle_within(X_mod, START_COL, all_rows), X_mod))

te_first <- cumsum(c(1, utils::head(sizes_te, -1)))
te_sor <- rep(seq_along(sizes_te), sizes_te)
Xa <- X_mod
Xa[, START_COL] <- X_mod[te_first, START_COL][
  sample.int(length(sizes_te))
][te_sor]
stopifnot(all(Xa[, START_COL] == Xa[te_first, START_COL][te_sor]))
drop_start <- ll_te -
  cond_loglik(
    eta_te - c_te$modifier + predict_X(bst$modifier, specs$modifier, Xa),
    te$case_,
    sizes_te
  )
cat(sprintf(
  "%s dealt between strata: drop %.2f (this fold)\n",
  START_COL,
  drop_start
))
# Production repeats this for every ranked column, N_PERM times in each of
# the N_FOLDS folds, and sums: that total over deer-years is the importance
# plot.

# 10 ---- Held-out gain of one split (xgb_heldout_split_gain, one tree) -------
# The gain of sections 6 and 7 is measured on the steps the split was chosen
# on, so even noise earns some. Here the same question is asked of held-out
# steps.
# Replay the fold model's trees IN THE ORDER THEY WERE GROWN, so each tree
# meets the held-out score as it stood just before it. For the first
# habitat tree that means round 1's move, modifier and hr trees.
#
# The gain here is the change in the log-likelihood itself, not the formula
# from section 6. That formula is a shortcut: it treats each candidate's
# probability as if it moved on its own, when in fact a step's
# probabilities move together, because they must add to 1. Good enough for
# choosing splits; this is the exact number.
tree1 <- function(b, s) {
  stats::predict(b, dm_of(te, s), outputmargin = TRUE, iterationrange = c(1, 1))
}
before <- tree1(bst$move, specs$move) +
  tree1(bst$modifier, specs$modifier) +
  tree1(bst$hr, specs$hr)

th <- xgb_tree_arrays(bst$hab)[[1]] # section 6 inlined this
one_tree <- xgboost::xgb.slice.Booster(bst$hab, 1, 1)
leaf_te <- as.integer(stats::predict(
  one_tree,
  dm_of(te, specs$hab),
  predleaf = TRUE
)) +
  1L

# Collapse node k: every row below it gets the node's single value, as if
# the tree had stopped there.
collapse <- function(vals, k) {
  if (th$left[k] != 0L) {
    vals[c(th$left[k], th$right[k])] <- LEARNING_RATE * th$weight[k]
  }
  vals
}
ll_with <- function(vals) {
  cond_loglik(before + vals[leaf_te], te$case_, sizes_te)
}
th_kids <- c(th$left[1], th$right[1])
k <- th_kids[th$left[th_kids] != 0L][1] # the first child that splits

# A child split's held-out gain: the tree as grown, minus the tree with that
# split collapsed. The root's: the tree cut back to its first split, minus
# no tree at all.
heldout_child <- ll_with(th$weight) - ll_with(collapse(th$weight, k))
heldout_root <- ll_with(collapse(collapse(th$weight, th_kids[1]), th_kids[2])) -
  cond_loglik(before, te$case_, sizes_te)

# Against production, which replays every tree of every booster
n_trees <- vapply(bst, xgboost::xgb.get.num.boosted.rounds, numeric(1))
replay <- xgb_heldout_split_gain(bst, specs, te, n_trees)
prod_rows <- replay$splits[
  replay$splits$booster == "hab" &
    replay$splits$tree == 1,
]
stopifnot(
  all.equal(heldout_root, prod_rows$heldout[prod_rows$depth == 0]),
  all.equal(heldout_child, prod_rows$heldout[prod_rows$depth == 1][1]),
  all.equal(replay$ll, ll_te, tolerance = 1e-6)
)
cat(sprintf(
  "habitat tree 1, held out: root %s %.3f | %s under it %.3f\n",
  th$variable[1],
  heldout_root,
  th$variable[k],
  heldout_child
))

# What the interaction plot is made of. Every child split's held-out gain is
# credited to its (root, child) pair; B under A and A under B are the same
# pair; production sums them over all trees and all folds and divides by
# deer-years. This fold alone, after WALK_ROUNDS rounds:
fold_pairs <- replay$splits |>
  dplyr::filter(depth == 1, booster == "hab", parent != variable) |>
  dplyr::mutate(a = pmin(parent, variable), b = pmax(parent, variable)) |>
  dplyr::group_by(a, b) |>
  dplyr::summarise(heldout = sum(heldout), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(heldout))
print(as.data.frame(utils::head(fold_pairs, 3)), row.names = FALSE, digits = 3)


# 11 ---- From numbers to decisions, and an audit of results/xgb/ -------------
# Three decisions read off the same scale, held-out log score per deer-year,
# against the same threshold:
#   model type   the simplest type within THRESHOLD of the best
#   variable     its importance must clear THRESHOLD
#   pair         its interaction log score must clear THRESHOLD
# DESIGN: the threshold and the complexity order are the same in every
# script that applies them, and the runner fits every type they rank.
plot_scripts <- c(
  "scripts/xgb/plot_models_xgb.R",
  "scripts/xgb/plot_importance_xgb.R",
  "scripts/xgb/plot_interactions_xgb.R"
)
thresholds <- vapply(plot_scripts, config_value, numeric(1), name = "THRESHOLD")
orders <- lapply(plot_scripts, config_value, name = "COMPLEXITY")
THRESHOLD <- thresholds[[1]]
COMPLEXITY <- orders[[1]]
SEASONS <- config_value("scripts/xgb/run_models_xgb.R", "SEASONS")
CONFIGS <- config_value("scripts/xgb/run_models_xgb.R", "CONFIGS")
stopifnot(
  all(thresholds == THRESHOLD),
  all(vapply(orders, identical, logical(1), COMPLEXITY)),
  setequal(CONFIGS, COMPLEXITY)
)
cat(sprintf(
  "threshold %g per deer; simplest first: %s\n",
  THRESHOLD,
  paste(COMPLEXITY, collapse = " < ")
))

files <- list.files(
  "results/xgb",
  pattern = "^compare_.*[.]rds$",
  full.names = TRUE
)
if (!length(files)) {
  cat("\nNo results in results/xgb/ yet; the audit needs run_models_xgb.R.\n")
} else {
  # The one loop in this file: a line per stored result.
  res <- purrr::map(files, readRDS)
  audit <- purrr::map_dfr(res, function(r) {
    dplyr::bind_cols(
      tibble::tibble(
        season = r$season,
        config = r$config,
        held_out = r$ll_cv,
        n_deer = r$n_deer_years,
        n_steps = r$n_steps,
        rounds = r$n_trees / length(r$specs),
        # Every block shares these settings, and every type has a movement
        # block, so read them from that.
        learning_rate = r$specs[[1]]$params$learning_rate,
        max_depth = r$specs[[1]]$params$max_depth,
        has_split_gain = !is.null(r$split_gain)
      ),
      blocks_of(r$specs)
    )
  })
  expected <- dplyr::left_join(
    tidyr::expand_grid(season = SEASONS, config = CONFIGS),
    dplyr::rename(design, config = type),
    by = "config"
  )
  joined <- dplyr::left_join(
    expected,
    audit,
    by = c("season", "config"),
    suffix = c("", "_fit")
  )
  stopifnot(
    # every season x type was fit, once
    nrow(audit) == nrow(expected),
    !anyNA(joined$held_out),
    # with the agreed blocks
    all(joined$habitat == joined$habitat_fit),
    all(joined$modifier == joined$modifier_fit),
    all(joined$hr_in_hab == joined$hr_in_hab_fit),
    all(joined$one_variable_trees == joined$one_variable_trees_fit),
    # with the production settings
    all(audit$rounds == N_ROUNDS),
    all(audit$learning_rate == LEARNING_RATE),
    all(audit$max_depth == MAX_DEPTH),
    # by the code that scores splits on held-out steps
    all(audit$has_split_gain)
  )
  # Every type of a season was fit to the same table, and for the walked
  # season that table is the current one, yearlings dropped.
  same_pool <- audit |>
    dplyr::group_by(season) |>
    dplyr::summarise(k = dplyr::n_distinct(paste(n_deer, n_steps)))
  walked <- audit[audit$season == SEASON, ]
  stopifnot(
    all(same_pool$k == 1),
    all(walked$n_deer == dplyr::n_distinct(d$key)),
    all(walked$n_steps == n_strata)
  )

  # The selection, coded the two ways the plotting scripts code it
  sel_a <- audit |>
    dplyr::group_by(season) |>
    dplyr::mutate(per_deer = (held_out - max(held_out)) / n_deer) |>
    dplyr::filter(per_deer > -THRESHOLD) |>
    dplyr::slice_min(match(config, COMPLEXITY), n = 1) |>
    dplyr::ungroup()
  sel_b <- audit |>
    dplyr::group_by(season) |>
    dplyr::mutate(
      within = (held_out - max(held_out)) / n_deer > -THRESHOLD,
      selected = config == COMPLEXITY[min(match(config[within], COMPLEXITY))]
    ) |>
    dplyr::filter(selected) |>
    dplyr::ungroup()
  stopifnot(identical(sel_a$config, sel_b$config))

  # What clears the threshold in each season's selected model
  cleared <- purrr::map_dfr(seq_len(nrow(sel_a)), function(i) {
    r <- res[[which(
      audit$season == sel_a$season[i] &
        audit$config == sel_a$config[i]
    )]]
    vars <- r$importance |>
      dplyr::filter(!(booster == "modifier" & scheme == "within")) |>
      dplyr::mutate(per_deer = cv / r$n_deer_years) |>
      dplyr::filter(per_deer >= THRESHOLD, !grepl("^shadow", variable))
    pairs <- r$split_gain$pairs |>
      dplyr::filter(parent != child) |>
      dplyr::mutate(a = pmin(parent, child), b = pmax(parent, child)) |>
      dplyr::group_by(a, b) |>
      dplyr::summarise(
        per_deer = sum(heldout) / r$n_deer_years,
        .groups = "drop"
      )
    tibble::tibble(
      season = r$season,
      model = r$config,
      variables = paste(vars$variable, collapse = ", "),
      best_pair = if (nrow(pairs)) {
        sprintf("%.2f", max(pairs$per_deer))
      } else {
        "(no pairs)"
      },
      pairs_cleared = sum(pairs$per_deer >= THRESHOLD)
    )
  })
  cat(
    sprintf(
      "\n%s (%s)\n  variables over %g: %s\n  best pair %s, pairs over %g: %d\n",
      cleared$season,
      cleared$model,
      THRESHOLD,
      cleared$variables,
      cleared$best_pair,
      THRESHOLD,
      cleared$pairs_cleared
    ),
    sep = ""
  )
}

cat(
  "\nAll assertions passed:",
  "this walkthrough computes what the helpers compute,",
  "and the stored results match the agreed design.\n"
)
