# runtime on 1 core with 187 GB of memory: ~ 35 min w 100 tracks, ~ 8 h w 1024
library('ctmm')    # for continuous-time movement modeling
library('dplyr')   # for data wrangling (e.g., %>%)
library('purrr')   # for functional programming (e.g., map(), map_dbl())
library('tidyr')   # for data wrangling (e.g., nested tibbles)
source('analysis/simulations/movement-model.R') # for consistency btwn scripts

tracks <- readRDS('simulations/tracks.rds') # simulated tracks
types <- c('constant', 'linear', 'cyclical', 'drifting', 'erratic')# trend types

days <- readRDS('simulations/days.rds')

days_summarized <-
  days %>% # times at which animals reached satiety on each day in each panel
  # add labelled telemetry tracks for each day
  left_join(tracks, by = 'day') %>%
  # filter tracks by when animal became full
  mutate(tel = map2(tel, t_expl,
                    \(.t, .max) data.frame(.t) %>% filter(t <= .max))) %>%
  # group by panel, animal, and day
  group_by(mean, variance, animal, day) %>%
  # group by panel and animal (not by day)
  group_by(mean, variance, animal) %>%
  mutate(t_start = lag(t_expl * 2), # start at end of previous explor. + return
         t_start = if_else(is.na(t_start), 0, t_start), # start 1st day at t = 0
         t_start = cumsum(t_start)) %>%
  # drop groups
  ungroup()

#' splitting into two sections bc `unnest()` can't work with very long vectors
N <- nrow(days_summarized)

days_summarized <-
  # expand nested datasets into a single big dataset
  bind_rows(unnest(days_summarized[1:floor(N / 2), ], tel),
            unnest(days_summarized[(floor(N / 2) + 1):N, ], tel)) %>%
  mutate(
    # change t to be the date and time of consecutive days, not "time of day"
    t = t + t_start,
    # make timestamp reflect the new values of t
    timestamp = as.POSIXct(t, tz = 'UTC',
                           origin = as.POSIXct('1970-01-01 00:00 UTC')),
    # id for each animal (needed to create a telemetry object)
    individual.local.identifier = paste(mean, ' mean,',
                                        variance, 'variance,',
                                        'animal', animal)) %>%
  # drop unnecessary columns
  select(-c(cell_id, new_cell, full, t_expl, t_start)) %>%
  # move the time columns to the beginning
  relocate(t, .after = day) %>%
  # create separated (nested) datasets for each animal in each panel
  nest(tel = -c(mean, variance, animal, mu, sigma2))

saveRDS(days_summarized, 'simulations/days-summarized.rds')

if(FALSE) {
  # diagnostics to ensure time was added up correctly
  library('ggplot2')
  theme_set(theme_bw())
  
  days_summarized <- days_summarized %>%
    mutate(t_max = map_dbl(tel, \(x) x %>%
                             group_by(day) %>%
                             summarize(t_max = max(t)) %>%
                             pull(t_max) %>%
                             mean()),
           t_min = map_dbl(tel, \(x) min(x$t)))
  
  unique(days_summarized$t_min)
  
  ggplot(days_summarized) +
    facet_grid(variance ~ mean) +
    geom_line(aes(animal, t_max), color = 'darkorange', lwd = 1)
  
  days_summarized %>%
    filter(mean == 'constant', variance == 'constant') %>%
    slice(1) %>%
    pull(tel) %>%
    first() %>%
    ggplot() +
    geom_path(aes(x, y, color = t, group = day))
}
