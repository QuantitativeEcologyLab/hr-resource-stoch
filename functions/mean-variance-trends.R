library('dplyr') # for data wrangling (e.g., %>%) and tibble()

# mean functions
# mean should be far from zero because a small mean with large variance implies the scale
# parameter is very large, which results in a heavy right tail and a large range in the
# simulated values.
# This issue would be clearer when parameterizing in terms of shape and scale, but the
# parameterization is often harder to plot and visualize.
min_mu <- 100
max_mu <- 200

#' animals reach their `REQUIRED` amount of energy after eating many times: 
#' `REQUIRED` should be high enough that all animals are forced to eat multiple times
#' `REQUIRED` should *not* be so high that HR estimates become too variable
#' the variance in HR depends on the variance in total food, i.e. sum of variances
REQUIRED <- max_mu * 10

mean_constant <- function(generations) {
  rep(mean(c(min_mu, max_mu)), length(generations))
}

mean_linear <- function(generations) {
  b <- (max_mu - min_mu) / (200 - 1)
  
  generations * b + min_mu - b
}

mean_cyclical <- function(generations) {
  sin(generations / 10) * (max_mu - min_mu) / 2 + mean(c(min_mu, max_mu))
}

mean_drifting <- function(generations, seed = 8) {
  set.seed(seed = seed)
  mu <- cumsum(rnorm(length(generations)))
  mu <- mu - min(mu) # set minimum to 0
  mu <- mu / max(mu) # shrink the range to [0, 1]
  mu <- mu * (max_mu - min_mu) # scale to the range of mu
  mu <- mu + min_mu # shift to the right intercept
  return(mu)
}

mean_erratic <- function(generations) {
  case_when(generations < quantile(generations, 0.25)~mean(c(min_mu, max_mu)),
            generations < quantile(generations, 2/3) ~ min_mu,
            TRUE ~ max_mu)
}

# variance functions
min_s2 <- 50^2
max_s2 <- 500^2

variance_constant <- function(generations) {
  rep(mean(c(min_s2, max_s2)), length(generations))
}

variance_linear <- function(generations) {
  b <- (max_s2 - min_s2) / (200 - 1)
  
  generations * b + min_s2 - b
}

variance_cyclical <- function(generations) {
  ((sin(generations / 15) + 1)/2)^3 * (max_s2 - min_s2) + min_s2
}

variance_drifting <- function(generations, seed = 1) {
  set.seed(seed)
  sigma2 <- length(generations) %>%
    rnorm() %>%
    cumsum() %>%
    abs()
  sigma2 <- sigma2 - min(sigma2) # ranges [0, Inf)
  sigma2 <- sigma2 / max(sigma2) # to [0, 1]
  sigma2 <- sigma2 * (max_s2 - min_s2) + min_s2 # to [min_s2, max_s2]
  return(sigma2)
}

variance_erratic <- function(generations) {
  case_when(generations < quantile(generations, 0.2) ~ (max_s2 + min_s2) / 1.75,
            generations < quantile(generations, 0.6) ~ max_s2,
            TRUE ~ min_s2)
}

# check ranges of functions
if(FALSE) {
  library('purrr')   # for functional programming
  library('ggplot2') # for fancy plots
  source('functions/rgamma2.R') # rgamma() parameterized by mean and variance
  theme_set(theme_bw())
  expand.grid(animal = 1:200,
              mean = c('constant', 'linear', 'cyclical', 'drifting', 'erratic'),
              variance = c('constant', 'linear', 'cyclical', 'drifting', 'erratic')) %>%
    mutate(mu = purrr::map2_dbl(mean, animal, \(.mean, .animal) {
      foo <- get(paste0('mean_', .mean))
      foo(1:200)[.animal]
    }),
    sigma2 = purrr::map2_dbl(variance, animal, \(.variance, .animal) {
      foo <- get(paste0('variance_', .variance))
      foo(1:200)[.animal]
      }),
    cv = sqrt(sigma2) / mu) %>%
    ggplot(aes(animal, cv)) +
    facet_grid(variance ~ mean) +
    geom_line() +
    scale_x_continuous(sec.axis = sec_axis(trans = 'identity',
                                           name = 'Mean', breaks = NULL)) +
    scale_y_continuous(limits = c(0, NA),
                       sec.axis = sec_axis('identity', 'Variance', NULL))
  
  tidyr::expand_grid(mu = c(min_mu, max_mu),
                     s2 = c(min_s2, max_s2)) %>%
    mutate(cv = sqrt(s2) / mu,
           k = s2 / mu^2, # shape parameter
           theta = s2 / mu, # scale parameter
           sample = map2(mu, s2, \(..x, ..y) rgamma2(mu = ..x, sigma2 = ..y, N = 1e6)),
           min = map_dbl(sample, min),
           max = map_dbl(sample, max),
           range = map_dbl(sample, \(x) range(x) %>% diff()))
}
