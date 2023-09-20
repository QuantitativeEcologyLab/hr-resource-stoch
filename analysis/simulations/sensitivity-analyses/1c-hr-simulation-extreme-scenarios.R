library('ctmm')    # for continuous-time movement modeling
library('dplyr')   # for data wrangling (e.g., %>%)
library('tidyr')   # for data wrangling (e.g., nested tibbles)
library('ggplot2') # for fancy plots
source('analysis/figures/default-figure-styling.R')
source('functions/rgamma2.R') # rgamma() parameterized by mean and variance
source('analysis/mean-variance-trends-panel-data.R') # means & variances
source('analysis/simulations/movement-model.R') # for consistency across scripts
source('functions/get_hr.R') # for extracting gaussian home range
set.seed(1) # for consistent results

tels <- readRDS('simulations/tracks.rds') # list of telemetry tracks
tracks <- readRDS('simulations/labelled-tracks.rds') # tibble of tracks
MAX_T <- max(tracks$t) # maximum amount of exploration time

WORST <- filter(d55, mu == min(mu)) %>% # lowest mean resources
  filter(sigma2 == max(sigma2)) %>% # with highest variance
  slice(1) # take the first row only
BEST <- filter(d55, mu == max(mu)) %>% # highest mean resources
  filter(sigma2 == min(sigma2)) %>% # with lowest variance
  slice(1) # take the first row only

days <-
  # modify WORST and BEST to follow the syntax used in 'offspring-simulations.R'
  transmute(bind_rows(WORST, BEST),
            animal,
            mu,
            sigma2,
            d = list(tracks),
            scenario = c('Worst case', 'Best case')) %>%
  unnest(d) %>% # unnest the datasets so we have a single, large tibble
  select(-timestamp) %>%
  # generate the food for each row from a gamma distribution
  mutate(food = rgamma2(mu = mu, sigma2 = sigma2, N = n()),
         # if the animal visits a new cell, it finds food, otherwise it doesn't
         food = if_else(new_cell, food, 0)) %>%
  # end the movement once the animal has reached satiety
  group_by(day, animal, scenario) %>%
  # calculate the total number of visits, total calories, and if animal is full
  mutate(satiety = cumsum(food), # use for diagnostics if animals don't get full
         full = satiety >= REQUIRED) %>% # did the animal reach the its needs?
  filter(cumsum(full) <= 1) %>% # full only once
  ungroup()

if(FALSE) {
  # check if the ends of each day are correct and make sense
  days_end <-
    days %>%
    group_by(scenario, day) %>%
    filter(full, ! duplicated(full)) %>% # take the 1st row where animal is full
    rename(t_expl = t) %>% # to avoid duplicated colnames with tracks
    # remove unneded columns (also avoids duplicated colnames with tracks)
    dplyr::select(-c(x, y, vx, vy, longitude, latitude, food))
  
  # check max fraction of time used (should be < 1)
  max(days_end$t_expl) / MAX_T
  
  # are animals are full only once/day? (should be == 1)
  sum(days_end$full) / (max(days$day) * 2)
  
  # plot of satiety over time by animal
  ggplot(days, aes(t, satiety, group = day)) +
    facet_wrap(~ scenario) +
    geom_line(alpha = 0.05) +
    geom_point(aes(t_expl), days_end, alpha = 0.1) +
    geom_hline(yintercept = REQUIRED, color = 'red') +
    geom_vline(xintercept = MAX_T, color = 'blue')
  
  # check distribution of animals
  ggplot(days_end, aes(scenario, t_expl)) +
    geom_hline(yintercept = MAX_T, color = 'red') +
    geom_violin(fill = 'forestgreen', alpha = 0.3) +
    geom_boxplot(fill = NA) +
    labs(x = '', y = 'Exploration time')
  
  # check home ranges of animals
  ggplot(days) +
    facet_grid(. ~ scenario) +
    coord_equal() +
    geom_hex(aes(longitude, latitude)) +
    scale_fill_distiller('Count', type = 'seq', na.value = 'transparent') +
    theme(legend.position = 'top')
}

# single estimates that eventually converge to the asymptote ----
days_summarized <-
  days %>%
  # find how long it took to reach satiety
  group_by(scenario, day) %>%
  nest(tel_day = -c(scenario, day)) %>%
  mutate(t_expl = map_dbl(tel_day, \(d) max(d$t))) %>%
  # add days sequentially
  group_by(scenario) %>%
  mutate(t_start = lag(2 * t_expl), # add the return time before next "day"
         t_start = if_else(is.na(t_start), 0, t_start), # start at 0, not NA
         t_start = cumsum(t_start), # make start times comsecutive
         tel_day = map2(day, t_expl,
                        \(i, te) tels$tel[[i]] %>% # extract tel for the day
                          data.frame() %>% # for filtering
                          filter(t <= te))) %>% # end tracks at satiety
  unnest(tel_day) %>% # make one big dataset
  mutate(t = t + t_start, # make times consecutive
         individual.local.identifier = scenario, # ctmm identifier
         timestamp = as.POSIXct(t, origin = '2000-01-01')) %>% # use new times 
  ungroup() # remove grouping by scenario

if(FALSE) {
  # check times are adding up correctly
  # best case should require less time
  # blue and red lines should have same length (within the pair)
  # black lines should be horizontal
  days_summarized %>%
    filter(day <= 10) %>%
    ggplot(aes(day, timestamp)) +
    facet_wrap(~ scenario, scales = 'free_y') +
    geom_line() +
    geom_line(aes(group = day), color = 'blue', lwd = 30) +
    geom_line(aes(day, timestamp + t_expl, group = day), color = 'red', lwd =30)
  
  days_summarized %>%
    filter(day <= 10) %>%
    ggplot() +
    facet_wrap(~ scenario) +
    coord_equal() +
    geom_path(aes(x, y, group = day), alpha = 0.5) +
    geom_point(aes(0, 0)) +
    geom_point(aes(x, y), filter(days_summarized, day <= 10, t == 0),
               color = 'red')
}

# estimate saturation curve of home range size over number of days
saturation_days <-
  expand_grid(n_days = (2^seq(1, log2(1e3), by = 0.2)) %>% round() %>% unique(),
              case = unique(days_summarized$scenario)) %>%
  mutate(data = map2(n_days, case,
                     \(.n, .case) filter(days_summarized,
                                         day <= .n,
                                         scenario == .case)),
         tel = map(data, as.telemetry), # convert data to telemetry for modeling
         theta = map(tel, \(x) ctmm.guess(data = x, interactive = FALSE)),
         m = map(1:n(), \(i) {
           cat('Fitting model', i, '\n')
           ctmm.fit(tel[[i]], theta[[i]])
         }), # fit movement model
         sigma = map_dbl(m, \(.m) ctmm:::area.covm(.m$sigma)), # var(position)
         hr = get_hr(.sigma = sigma, quantile = 0.95)) # Gaussian home range

saveRDS(saturation_days, 'simulations/hr-saturation-days.rds')

saturation_days %>%
  select(case, n_days, sigma, hr) %>%
  readr::write_csv('simulations/hr-saturation-days.csv')

# to avoid having to run all code above
saturation_days <- readRDS('simulations/hr-saturation-days.rds')

p_hr_days <-
  ggplot(saturation_days, aes(n_days, hr)) +
  facet_wrap(~ case, nrow = 1) +
  geom_vline(xintercept = 100, color = 'red') +
  geom_smooth(method = 'gam', color = 'black',
              formula = y ~ s(x, bs = 'cs', k = 10),
              method.args = list(family = Gamma(link = 'log'))) +
  geom_point(alpha = 0.3) +
  scale_x_continuous(expression(bold(Number~of~days~sampled~(log[2]~scale))),
                     trans = 'log2', breaks = c(2, 16, 128, 1024),
                     limits = c(2, 1100)) +
  scale_y_continuous(trans = 'log2',
                     expression(bold(atop(Estimated~'space-use',
                                          requirements~(log[2]~scale)))))

ggsave('figures/hr-over-days.png', plot = p_hr_days,
       width = 6, height = 3, dpi = 'print')
