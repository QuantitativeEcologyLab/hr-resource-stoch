# runtime on 1 core with 187 GB of memory: < 3 minutes for 100 tracks, < 20 min for 2^10
library('ctmm')  # for continuous-time movement modeling
library('dplyr') # for data wrangling (e.g., %>%)
library('purrr') # for functional programming (e.g., map(), map_dbl())
library('tidyr') # for data wrangling (e.g., nested tibbles)
library('ggplot2') # for fancy plots
source('analysis/figures/default-figure-styling.R')
source('analysis/simulations/mean-variance-trends-panel-data.R')
source('functions/rgamma2.R') # rgamma() parameterized by mean and variance

tracks <- readRDS('simulations/labelled-tracks.rds') %>% # movement tracks
  filter(day <= 200) # ~100 are sufficient for stable estimates
types <- c('constant', 'linear', 'cyclical', 'drifting', 'erratic') # trend types

# using for loops to reduce computation costs
for(i in 1:(length(types))) {
  mean_fcn <- get(paste0('mean_', types[i]))
  cat('Running', types[i], 'mean\n')
  for(j in 1:(length(types))) {
    cat('  Running', types[j], 'variance\n')
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

# as a check, plot mean exploration time for each scenario +/- 2 SE
for(filename in list.files('simulations', pattern = '*-replicates.rds')) {
  readRDS(filename) %>%
    mutate(mean = paste(mean, 'mean'),
           variance = paste(variance, 'var')) %>%
    group_by(mean, variance, animal) %>%
    summarize(est = mean(t_expl),
              # add 2 SEs around the mean
              se = sd(t_expl) / sqrt(n()),
              lwr = max(0, est - 2 * se),
              upr = est + 2 * se,
              .groups = 'drop') %>%
    ggplot() +
    coord_cartesian(ylim = c(2e3, 12e3)) +
    facet_grid(variance ~ mean) +
    geom_ribbon(aes(animal, ymin = lwr, ymax = upr), fill = 'black') +
    geom_line(aes(animal, est), color = '#CCBB44') +
    labs(x = 'Time', y = 'Estimated time to satiety',
         title = paste(n_distinct(days$day), 'replicates'))
  
  ggsave(paste0('figures/5-by-5-sensitivity-analysis/time-to-satiety-',
                n_distinct(days$day), '-replicates.png'),
         width = 10, height = 5.625, units = 'in', dpi = 'print', bg = 'white')
}

all <-
  tibble(d = map_dfr(list.files('simulations', pattern = '*-replicates.rds',
                                full.names = TRUE),
                     function(fn) {
                       readRDS(fn) %>%
                         mutate(mean = paste(mean, 'mean'),
                                variance = paste(variance, 'var'),
                                replicates = max(day)) %>%
                         group_by(mean, variance, animal, replicates) %>%
                         summarize(est = mean(t_expl),
                                   .groups = 'drop')
                     })) %>%
  unnest(d) %>%
  mutate(est = est / (1 %#% 'hours')) %>%
  pivot_wider(names_from = replicates, values_from = est) %>%
  pivot_longer(cols = as.character(2^(4:9)), values_to = 'est',
               names_to = 'replicates') %>%
  mutate(diff_est = est - `1024`,
         rel_est = est / `1024`,
         replicates = factor(replicates, levels = 2^(4:9)))

ggplot(all) +
  facet_wrap(~ replicates) +
  geom_line(aes(animal, diff_est, group = paste(mean, variance)), alpha = 0.1) +
  geom_hline(yintercept = 0, color = 'red') +
  labs(x = 'Time', y = 'Difference in average daily time required (hours)')

ggsave('figures/5-by-5-sensitivity-analysis/time-to-satiety-difference.png',
       width = 10, height = 5.625, units = 'in', dpi = 'print', bg = 'white')

ggplot(all) +
  facet_wrap(~ replicates) +
  geom_line(aes(animal, rel_est, group = paste(mean, variance)), alpha = 0.1) +
  geom_hline(yintercept = 1, color = 'red') +
  scale_y_continuous(trans = 'log2', breaks = 2^c(-0.5, 0, 0.5),
                     labels = parse(text = c('2^-0.5', '0', '2^0.5')))+
  labs(x = 'Time',
       y = expression(bold(Relative~difference~'in'~time~to~satiety~(log[2]~scale))))

ggsave('figures/5-by-5-sensitivity-analysis/time-to-satiety-rel-difference.png',
       width = 10, height = 5.625, units = 'in', dpi = 'print', bg = 'white')
