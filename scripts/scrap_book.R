env_local <- terra::unwrap(deer_input$crop_env)
ndvi_local <- terra::unwrap(deer_input$crop_ndvi)

stp_data = deer_input$stp_test
x_median = deer_input$x_median
y_median = deer_input$y_median
env_test = env_local
ndvi_test = ndvi_local
issf_train = model_sims[[2]]

median_pt <- terra::vect(
  cbind(x_median, y_median),
  crs = terra::crs(env_test[[1]])
)
env_test$HR <- NA
env_test$HR <- terra::distance(env_test$HR, median_pt) / 1000

bursts <- unique(stp_data$burst_)

burst_data <- stp_data |> dplyr::filter(burst_ == 1)
current_month <- NULL
step_results <- vector("list", nrow(burst_data))


# Update NDVI only when month changes
mo <- lubridate::month(burst_data$t1_[1])

env_test$ndvi <- terra::resample(
  ndvi_test[[mo]],
  env_test,
  method = "near"
)


# Build start point from observed location
start_pt <- burst_data[1, c("x1_", "y1_", "t1_")] |>
  amt::make_track(
    .x = x1_,
    .y = y1_,
    .t = t1_,
    crs = terra::crs(env_test)
  ) |>
  amt::make_start() |>
  amt::mutate(dt = lubridate::hours(4))

# Build redistribution kernel
kernel <- tryCatch(
  amt::redistribution_kernel(
    x = issf_train,
    map = env_test,
    fun = function(xy, map) {
      xy |>
        amt::extract_covariates(map, where = "both") |>
        amt::time_of_day(
          include.crepuscule = FALSE,
          where = "both"
        ) |>
        dplyr::mutate(
          tod_start_day = as.integer(tod_start_ == "day"),
          tod_start_night = as.integer(tod_start_ == "night"),
          days = lubridate::yday(t2_) - min(lubridate::yday(t2_)) + 1
        )
    },
    start = start_pt,
    landscape = "discrete",
    as.rast = FALSE
  ),
  error = function(e) NULL
)

z <- amt::redistribution_kernel(
  x = issf_train,
  map = env_test,
  fun = function(xy, map) {
    xy |>
      amt::extract_covariates(map, where = "both") |>
      amt::time_of_day(
        include.crepuscule = FALSE,
        where = "both"
      ) |>
      dplyr::mutate(
        tod_start_day = as.integer(tod_start_ == "day"),
        tod_start_night = as.integer(tod_start_ == "night"),
        days = lubridate::yday(t2_) - min(lubridate::yday(t2_)) + 1
      )
  },
  start = start_pt,
  landscape = "discrete",
  as.rast = T
)

obs_pt <- cbind(burst_data$x2_[1], burst_data$y2_[1])
p <- terra::extract(kernel_rast$redistribution.kernel, obs_pt)[1, 2]

lp <- if (!is.na(p) && p > 0) log(p) else NA_real_

library(tidyverse)

# Read all loglik results
loglik_files <- list.files(
  "results",
  pattern = "results_loglik_\\d+\\.rds",
  full.names = TRUE
)

results_all <- map_dfr(loglik_files, function(f) {
  row_no <- as.integer(str_extract(basename(f), "\\d+"))
  df <- readRDS(f)
  df$deer_row <- row_no
  df
})

# Keep only models that beat null by >= 4
candidates <- results_all %>%
  filter(!is.na(delta_logp), delta_logp >= 4)

# For deer with no candidates, keep the null
null_only_deer <- setdiff(
  unique(results_all$deer_row),
  unique(candidates$deer_row)
)

results_filtered <- bind_rows(
  candidates,
  results_all %>% filter(deer_row %in% null_only_deer, model == 2)
) %>%
  select(deer_row, model, total_logp, n_steps, delta_logp)

deer_mvt <- readRDS("data_deer_1_119.rds")
stp <- deer_mvt$stp.var.train[[117]]

# What time range does this deer cover?
range(stp$t1_)

# What NDVI dates are available?
ndvi <- rast("NDVI_2018.tif")
terra::time(ndvi)
