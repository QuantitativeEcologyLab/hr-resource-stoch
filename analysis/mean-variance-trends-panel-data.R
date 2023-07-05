#' created a separate script for `d55` since it's used in multiple scripts

library('dplyr') # for data wrangling
library('tidyr') # for data wrangling (expand_grid())
library('purrr') # for functional programming (e.g., map_*())
source('functions/mean-variance-trends.R')  # functions to generate means and variances

N <- 200 # number of animals
types <- c('constant', 'linear', 'cyclical', 'drifting', 'erratic') # functions types

d55 <-
  expand_grid(mean = types,
              variance = types) %>%
  # assign trends in means and variances depending on the type
  mutate(t = map2(mean, variance,
                  \(.x, .y) tibble(animal = 1:N,
                                   mu = get(paste0('mean_', .x))(animal),
                                   sigma2 = get(paste0('variance_', .y))(animal)))) %>%
  unnest(t) %>%
  mutate(mean = factor(mean, levels = types),
         variance = factor(variance, levels = types))
