library(tidyverse)
library(amt)

load("sw_deer_tracks.rda")

z <- sw_deer_tracks[1, ]$data[[1]]

z2 <- z %>%
  track_resample() %>%
  filter_min_n_burst(n = 3) %>%
  steps_by_burst() %>%
  time_of_day(include.crepuscule = F)

z3 <- z2 %>%
  random_steps(n = 9) %>%
  time_of_day(include.crepuscule = F) %>%
  mutate(log_sl_ = log(sl_))
