#' @description
#' Explore the fitted deer-movement GAMs (output of fit_GAM.R).
#'
#' Section 1 (here): the k / shrinkage audit. fit_GAM.R writes one row per smooth
#' per (deer, model) fit to results/k_audit_gam.rds, each labelled with a status
#' from gam_smooth_diag():
#'   * "ok"          smooth is supported and below its basis ceiling
#'   * "k-bound"     edf pressed against k' -> consider a higher k
#'   * "near-linear" collapsed to the penalty null space (a straight line)
#'   * "removed"     shrunk to ~zero (no support; expected for cc / re / fs terms
#'                   under ordinary REML, see gam_smooth_diag docs)
#' We draw one bar plot per smooth showing how often it lands in each status
#' across all deer x model fits, and save it to plots/.
#'
#' Section 2 (forthcoming): relative selection strength (RSS) plots built from
#' the saved results_gam_<key>.rds models. Stubbed at the bottom for now.
#'
#' Usage:  Rscript scripts/gam/explore_gam_results.R

# Load packages ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(mgcv) # predict.gam on the saved cox.ph fits (Section 2)
  library(patchwork) # assemble the RSS panels (Section 2)
})

# prepare_gam_data(), LANDCOVER_LEVELS (used in Section 2)
source("scripts/helper_functions.R")

# Config ----------------------------------------------------------------------
audit_path <- "results/k_audit_gam.rds"
out_dir <- "plots"
dir.create(out_dir, showWarnings = FALSE)

if (!file.exists(audit_path)) {
  stop(sprintf("Audit file not found: %s (run fit_GAM.R first).", audit_path))
}

# Status order + palette. Ordered roughly worst-to-best so the bars read
# consistently; colours flag the "needs attention" states (k-bound = raise k;
# removed = dropped).
status_levels <- c("ok", "k-bound", "near-linear", "removed")
status_cols <- c(
  "ok" = "#4daf4a",
  "k-bound" = "#ff7f00",
  "near-linear" = "#377eb8",
  "removed" = "#999999"
)

# Human-readable facet titles, in the order we want the panels laid out
# (movement, then home-range terms, then resource terms). Any smooth not listed
# here falls back to its raw term label, so the script survives formula changes.
smooth_labels <- c(
  "s(tod_):sl_" = "Movement: s(tod) x sl (cc)",
  "s(HR_edge_end)" = "HR edge",
  "s(HR_center_end)" = "HR center",
  "s(ndvi_end)" = "NDVI: global response",
  "s(ndvi_end,wiscland_end)" = "NDVI x landcover (fs)",
  "s(wiscland_end)" = "Landcover intercept (re, winter)",
  "s(HR_center_end,wiscland_end)" = "HR center x landcover (fs)"
)

# Load audit ------------------------------------------------------------------
audit <- readRDS(audit_path)

n_deer <- dplyr::n_distinct(audit$key)
n_fits <- nrow(dplyr::distinct(audit, key, model))
cat(sprintf(
  "Loaded %d smooth rows from %d deer (%d deer x model fits).\n",
  nrow(audit),
  n_deer,
  n_fits
))

# Order smooths (known first, then any extras) and the status factor.
smooth_order <- c(
  names(smooth_labels),
  setdiff(sort(unique(audit$smooth)), names(smooth_labels))
)
facet_labels <- smooth_labels
extra <- setdiff(smooth_order, names(facet_labels))
facet_labels[extra] <- extra # unknown smooths label themselves

audit <- audit |>
  mutate(
    smooth = factor(smooth, levels = smooth_order),
    status = factor(status, levels = status_levels)
  )

# Console summary: smooth x status counts ------------------------------------
status_table <- audit |>
  count(smooth, status, .drop = FALSE) |>
  pivot_wider(names_from = status, values_from = n, values_fill = 0)
cat("\nStatus frequency per smooth:\n")
print(as.data.frame(status_table), row.names = FALSE)

# Bar plot: one panel per smooth ---------------------------------------------
counts <- audit |>
  count(smooth, status, .drop = FALSE)

p <- ggplot(counts, aes(status, n, fill = status)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = ifelse(n > 0, n, "")),
    vjust = -0.3,
    size = 2.8
  ) +
  facet_wrap(
    ~smooth,
    scales = "free_y",
    labeller = as_labeller(facet_labels)
  ) +
  scale_fill_manual(values = status_cols, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "GAM smooth k / shrinkage audit",
    subtitle = sprintf(
      "Status frequency per smooth across %d deer x model fits",
      n_fits
    ),
    x = NULL,
    y = "Number of fits",
    fill = "Status"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

out_file <- file.path(out_dir, "gam_k_audit_status.png")
ggsave(out_file, p, width = 10, height = 7, dpi = 150)
cat(sprintf("\nSaved status bar plot -> %s\n", out_file))

# Relative selection strength (RSS) plots ------------------------------------
# For one well-behaved non-winter deer (all smooths "ok"), show the RSS each
# fitted smoother implies. RSS(x) = exp(f(x) - f(x_ref)) -- selection at
# covariate value x relative to a reference, other covariates held fixed. Each
# term's contribution is pulled with predict(type = "terms"):
#   * tod x sl   (model 1): s(tod_, by = sl_) movement modulation, shown at 3
#                           step lengths (it scales with sl_), centred over the day
#   * HR edge    (model 2): s(HR_edge_end), relative to its median
#   * HR center  (model 3): s(HR_center_end), relative to its median
#   * NDVI x LC  (model 5): global s(ndvi_end) + per-class fs together (the GS
#                           block), i.e. landcover + NDVI selection vs forest

# Pick the deer: first non-winter deer whose smooths are all "ok". Override by
# setting rss_key to a specific "<id>_<season>_<year>".
rss_key <- NA_character_
if (is.na(rss_key)) {
  rss_key <- audit |>
    mutate(season = map_chr(str_split(key, "_"), 2)) |>
    filter(season != "nb") |>
    group_by(key) |>
    summarise(all_ok = all(status == "ok"), .groups = "drop") |>
    filter(all_ok) |>
    pull(key) |>
    first()
}
if (is.na(rss_key)) {
  stop("No non-winter deer with all-ok smooths found.")
}
cat(sprintf("\nRSS plots for deer: %s\n", rss_key))

res <- readRDS(sprintf("results/results_gam_%s.rds", rss_key))
dat <- prepare_gam_data(
  readRDS(sprintf("data/tracks/data_%s.rds", rss_key))[["stp.var"]][[1]]
)

# Fill any covariate a model needs but a panel doesn't vary, at representative
# (median / reference) values, so predict() always has every column.
fill_const <- function(df) {
  cst <- list(
    sl_ = median(dat$sl_),
    ta_ = 0,
    tod_ = 12,
    HR_edge_end = median(dat$HR_edge_end),
    HR_center_end = median(dat$HR_center_end),
    ndvi_end = median(dat$ndvi_end, na.rm = TRUE),
    wiscland_end = factor("forest", levels = LANDCOVER_LEVELS)
  )
  for (nm in names(cst)) if (is.null(df[[nm]])) df[[nm]] <- cst[[nm]]
  df
}

# Summed contribution of the term(s) whose label matches `pattern`.
term_contrib <- function(gam, nd, pattern) {
  tt <- predict(gam, fill_const(nd), type = "terms")
  cols <- grep(pattern, colnames(tt), value = TRUE)
  rowSums(tt[, cols, drop = FALSE])
}

theme_rss <- theme_bw(base_size = 10) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# Panel 1: tod x sl (model 1). The term s(tod_)*sl_ scales with step length, so
# show it at the 25/50/75th sl percentiles; each curve is centred on its own
# daily geometric mean (RSS = 1 = that step length's average over the day).
sl_q <- quantile(dat$sl_, c(.25, .5, .75))
g1 <- tidyr::expand_grid(
  tod_ = seq(0, 24, length.out = 200),
  sl_lab = factor(names(sl_q), levels = names(sl_q))
) |>
  mutate(sl_ = sl_q[as.character(sl_lab)])
g1$contrib <- term_contrib(res[["1"]]$gam, g1, "tod_")
g1 <- g1 |>
  group_by(sl_lab) |>
  mutate(rss = exp(contrib - mean(contrib))) |>
  ungroup()
p1 <- ggplot(g1, aes(tod_, rss, colour = sl_lab)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_line(linewidth = 0.8) +
  scale_colour_brewer(palette = "Dark2", name = "step length\n(percentile)") +
  scale_x_continuous(breaks = seq(0, 24, 6)) +
  labs(title = "tod x sl  (model 1)", x = "time of day (h)", y = "RSS") +
  theme_rss

# Panel 2: HR edge (model 2), relative to its median.
g2 <- tibble(HR_edge_end = seq(
  min(dat$HR_edge_end), max(dat$HR_edge_end),
  length.out = 200
))
ref2 <- term_contrib(
  res[["2"]]$gam, tibble(HR_edge_end = median(dat$HR_edge_end)), "HR_edge"
)
g2$rss <- exp(term_contrib(res[["2"]]$gam, g2, "HR_edge") - ref2)
p2 <- ggplot(g2, aes(HR_edge_end, rss)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_line(linewidth = 0.8, colour = "#1b9e77") +
  labs(title = "HR edge  (model 2)", x = "HR_edge_end", y = "RSS (vs median)") +
  theme_rss

# Panel 3: HR center (model 3), relative to its median.
g3 <- tibble(HR_center_end = seq(
  min(dat$HR_center_end), max(dat$HR_center_end),
  length.out = 200
))
ref3 <- term_contrib(
  res[["3"]]$gam, tibble(HR_center_end = median(dat$HR_center_end)), "HR_center"
)
g3$rss <- exp(term_contrib(res[["3"]]$gam, g3, "HR_center") - ref3)
p3 <- ggplot(g3, aes(HR_center_end, rss)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_line(linewidth = 0.8, colour = "#d95f02") +
  labs(
    title = "HR center  (model 3)",
    x = "HR_center_end",
    y = "RSS (vs median)"
  ) +
  theme_rss

# Panel 4: the NDVI x landcover GS block (model 5) = global s(ndvi_end) +
# per-class fs deviations, evaluated for every landcover class over the observed
# NDVI range, relative to forest at median NDVI. In the clean GS every class gets
# its own NDVI curve, and the fs also carries the class intercept -- so vertical
# offset is landcover selection and the shape is the class's NDVI response.
ndvi_rng <- range(dat$ndvi_end, na.rm = TRUE)
g4 <- tidyr::expand_grid(
  ndvi_end = seq(ndvi_rng[1], ndvi_rng[2], length.out = 200),
  wiscland_end = factor(LANDCOVER_LEVELS, levels = LANDCOVER_LEVELS)
)
ref4 <- term_contrib(
  res[["5"]]$gam,
  tibble(
    wiscland_end = factor("forest", levels = LANDCOVER_LEVELS),
    ndvi_end = median(dat$ndvi_end, na.rm = TRUE)
  ),
  "ndvi_end"
)
g4$rss <- exp(term_contrib(res[["5"]]$gam, g4, "ndvi_end") - ref4)
p4 <- ggplot(g4, aes(ndvi_end, rss, colour = wiscland_end)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_line(linewidth = 0.8) +
  scale_colour_viridis_d(option = "turbo", name = "landcover") +
  scale_y_log10() +
  labs(
    title = "NDVI x landcover  (model 5)",
    x = "NDVI",
    y = "RSS (vs forest, log)"
  ) +
  theme_rss

# Assemble + save -------------------------------------------------------------
fig <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = sprintf("GAM relative selection strength - deer %s", rss_key),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )
out_rss <- file.path(out_dir, sprintf("gam_rss_%s.png", rss_key))
ggsave(out_rss, fig, width = 12, height = 9, dpi = 150)
cat(sprintf("Saved RSS plots -> %s\n", out_rss))
