# runtime on 1 core with 187 GB of memory: < 3 minutes for 100 tracks, < 20 min for 1e3
library('ctmm')  # for continuous-time movement modeling
library('dplyr') # for data wrangling (e.g., %>%)
library('purrr') # for functional programming (e.g., map(), map_dbl())
library('tidyr') # for data wrangling (e.g., nested tibbles)
source('analysis/figures/mean-variance-trends-panel-data.R') # trends in E(U) and V(U)
source('functions/rgamma2.R') # rgamma() parameterized by mean and variance

tracks <- readRDS('simulations/labelled-tracks.rds') %>% # movement tracks
  filter(day <= 100) # 100 are sufficient for stable estimates
types <- c('constant', 'linear', 'cyclical', 'drifting', 'erratic') # trend types

# using for loops to reduce computation costs
for(i in 1:(length(types))) {
  mean_fcn <- get(paste0('mean_', types[i]))
  for(j in 1:(length(types))) {
    var_fcn <- get(paste0('variance_', types[j]))
    # calculate the number of visits necessary for each animal in each day
    days <-
      # tibble of animals and movement dataset
      d55 %>%
      filter(mean == types[i], variance == types[j]) %>%
      mutate(d = list(tracks)) %>%
      unnest(d) %>% # unnest the datasets so we have a single, large tibble
      # generate the food for each row from a gamma distribution
      mutate(food = rgamma2(mu = mu, sigma2 = sigma2, N = n()),
             # if the animal visits a new cell, it finds food, otherwise it doesn't
             food = if_else(new_cell, food, 0)) %>%
      # end the movement once the animal has reached satiety
      group_by(day, animal) %>%
      # calculate the total number of visits, total calories, and if animal is full
      mutate(satiety = cumsum(food), # use for diagnostics if animals don't get full
             full = satiety >= REQUIRED) %>% # did the animal reach the required needs?
      filter(full, ! duplicated(full)) %>% # take the 1st row where the animal is full
      rename(t_expl = t) %>% # to avoid duplicated colnames with tracks
      # remove unneeded columns (also avoids duplicated colnames with tracks)
      dplyr::select(-c(x, y, vx, vy, timestamp, longitude, latitude, food))
    
    # save the simulation
    saveRDS(days,
            paste0('simulations/days-5-by-5/', types[i], '-mean-',
                   types[j], '-variance-days-exploration.rds'))
  } # close `for` loop of variances
} # close `for` loop of means

# create a single tibble from all the simulations
days <-
  expand_grid(mean = types, variance = types) %>%
  mutate(mean = factor(mean, levels = types),
         variance = factor(variance, levels = types),
         d = map2(mean, variance, # add the simulated tracks
                  \(x, y) paste0('simulations/days-5-by-5/',
                                 x, '-mean-',
                                 y, '-variance-days-exploration.rds') %>%
                    readRDS() %>%
                    select(-c(mean, variance)))) %>%
  unnest(d)

if(any(! days$full)) {
  warning(paste('CAUTION:', sum(days$full), 'animals did *NOT* reach satiety!'))
}

saveRDS(days, 'simulations/days.rds') # save all simulations together

if(FALSE) {
  # as a check, plot mean exploration time for each scenario +/- 2 SE
  library('ggplot2') # for fancy plots
  theme_set(theme_bw())
  
  days %>%
    group_by(mean, variance, animal) %>%
    summarize(est = mean(t_expl),
              se = sd(t_expl) / sqrt(n()),
              lwr = max(0, est - 2 * se), # keep > 0
              upr = est + 2 * se,
              rows = n(),
              .groups = 'keep') %>%
    ggplot() +
    facet_grid(variance ~ mean) +
    geom_ribbon(aes(animal, ymin = lwr, ymax = upr), alpha = 0.4) +
    geom_line(aes(animal, est)) +
    labs(x = NULL, y = 'Estimated time to satiety')
}
