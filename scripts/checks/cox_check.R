#' @description
#' Does the xgb path's objective reproduce xgboost's built-in Cox objective
#' (survival:cox)?
#'
#' The conditional logit is the Cox partial likelihood with one event per
#' stratum: mgcv fits it as cox.ph with cbind(times, stratum), and
#' survival::clogit() fits it through coxph(). xgboost's survival:cox has no
#' strata, so the xgb path writes its own objective: push p - y and
#' H p (1 - p), with p computed within each step. survival:cox computes the
#' same push and H over its risk sets, so on data where a risk set is one
#' step the two should be the same thing. Three comparisons, 10 trees each:
#'
#'   1  One step, used endpoint first. survival:cox's risk set is exactly
#'      that step, so both must grow identical trees.
#'   2  The same step, used endpoint last. Every row has the same time, and
#'      survival:cox counts an event only for rows it reaches after that
#'      event, so the available rows get no push and its trees cannot split.
#'      Ours does not depend on row order. The xgb tables always put the used
#'      endpoint first, the order in which the two agree.
#'   3  Two steps. survival:cox compares each used endpoint with both steps'
#'      endpoints, ours with its own step's only, so the trees must differ.
#'
#' Input: data/xgb/pooled_start_<season>.rds
#'
#' Configuration: edit the block below before running.

# Configuration ---------------------------------------------------------------
SEASON <- "pf"
FEATS <- c("sl_", "elevation_end", "forest_edge_end", "ndvi_end")
PARAMS <- list(max_depth = 2, eta = 0.3, lambda = 1, min_child_weight = 0,
               tree_method = "exact", nthread = 1)
N_ROUNDS <- 10L
# Thresholds and leaf values are stored as single-precision floats
TOL <- 1e-6

# Load packages ---------------------------------------------------------------
library(xgboost)

# helper functions
source("scripts/helper_functions.R")

d <- readRDS(sprintf("data/xgb/pooled_start_%s.rds", SEASON))

# The two fits ----------------------------------------------------------------
# Built-in Cox: label +1 is an event at time 1, -1 is censored at time 1.
# base_score 1 (a start score of log 1 = 0) stops xgboost estimating one
# from the labels; base_margin already starts every row at 0.
fit_cox <- function(s) {
  dm <- xgb.DMatrix(as.matrix(s[, FEATS]),
                    label = ifelse(s$case_ == 1, 1, -1),
                    base_margin = rep(0, nrow(s)))
  xgb.train(c(PARAMS, objective = "survival:cox", base_score = 1), dm,
            nrounds = N_ROUNDS, verbose = 0)
}

# Ours: the objective of fit_xgb_boosters(), with the score tracked outside
# xgboost the same way
fit_ours <- function(s) {
  sizes <- xgb_strata_sizes(s)
  dm <- xgb.DMatrix(as.matrix(s[, FEATS]), label = s$case_,
                    base_margin = rep(0, nrow(s)))
  eta_now <- rep(0, nrow(s))
  bst <- NULL
  for (i in seq_len(N_ROUNDS)) {
    obj <- function(preds, dtrain) {
      p <- xgb_softmax(eta_now, sizes)
      list(grad = p - s$case_, hess = p * (1 - p))
    }
    bst <- xgb_grow_one(bst, dm, PARAMS, i - 1L, obj)
    eta_now <- predict(bst, dm, outputmargin = TRUE)
  }
  bst
}

# Comparing trees -------------------------------------------------------------
# A model's text dump split into its wording (variables and structure) and
# its numbers (node ids, thresholds, leaf values)
NUM_RE <- "-?[0-9]+[.]?[0-9]*(e[-+]?[0-9]+)?"
dump_parts <- function(bst) {
  txt <- paste(xgb.dump(bst), collapse = "\n")
  list(words = gsub(NUM_RE, "#", txt),
       numbers = as.numeric(regmatches(txt, gregexpr(NUM_RE, txt))[[1]]))
}
same_trees <- function(a, b) {
  pa <- dump_parts(a)
  pb <- dump_parts(b)
  identical(pa$words, pb$words) &&
    max(abs(pa$numbers - pb$numbers)) < TOL
}
largest_gap <- function(a, b) {
  max(abs(dump_parts(a)$numbers - dump_parts(b)$numbers))
}
n_splits <- function(bst) {
  sum(grepl("<", xgb.dump(bst)))
}
show_first_tree <- function(a, b) {
  first <- function(bst) {
    x <- xgb.dump(bst)
    x[2:(which(x == "booster[1]") - 1)]
  }
  fa <- first(a)
  fb <- first(b)
  n <- max(length(fa), length(fb))
  length(fa) <- n
  length(fb) <- n
  cat(sprintf("    %-44s %s\n", c("ours", fa), c("survival:cox", fb)),
      sep = "")
}

# 1. One step, used endpoint first ---------------------------------------------
one <- d[d$stratum == d$stratum[1], ]
stopifnot(nrow(one) == 101, one$case_[1] == 1, sum(one$case_) == 1)
ours_1 <- fit_ours(one)
cox_1 <- fit_cox(one)
cat("\n1. one step, used endpoint first\n")
cat(sprintf("   largest gap in any threshold or leaf value: %.1g\n",
            largest_gap(ours_1, cox_1)))
cat("   first tree:\n")
show_first_tree(ours_1, cox_1)
stopifnot(same_trees(ours_1, cox_1), n_splits(ours_1) > 0)

# 2. The same step, used endpoint last -----------------------------------------
one_last <- one[c(2:nrow(one), 1), ]
ours_2 <- fit_ours(one_last)
cox_2 <- fit_cox(one_last)
cat("\n2. the same step, used endpoint last\n")
cat(sprintf("   splits in all %d trees: ours %d, survival:cox %d\n", N_ROUNDS,
            n_splits(ours_2), n_splits(cox_2)))
stopifnot(same_trees(ours_2, ours_1), n_splits(cox_2) == 0)

# 3. Two steps -----------------------------------------------------------------
two <- d[d$stratum %in% unique(d$stratum)[1:2], ]
stopifnot(sum(two$case_) == 2)
ours_3 <- fit_ours(two)
cox_3 <- fit_cox(two)
cat("\n3. two steps\n")
cat("   first tree:\n")
show_first_tree(ours_3, cox_3)
stopifnot(!same_trees(ours_3, cox_3))

cat("\nAll Cox checks passed: on one step our objective is survival:cox.\n")
