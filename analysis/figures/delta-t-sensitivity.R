library('terra')   # for working with rasters
library('ctmm')    # for generating movement models
library('dplyr')   # for data wrangling (%>%, mutate(), slice(), filter(), etc.)
library('tidyr')   # for data wrangling (unnest, etc.)
library('purrr')   # for functional programming (map_*(), etc.)
library('ggplot2') # for fancy plots
library('cowplot') # for multi-panel fancy plots
source('analysis/simulations/movement-model.R') # for consistency between scripts
theme_set(theme_bw() + # change default ggplot theme
            theme(legend.position = 'none', strip.background = element_blank(),
                  strip.text = element_blank()))

# simulate some movement
get_tracks <- function(.seed) {
  simulate(model, # ctmm movement model
           t = 0:(max(SAMPLES)), # 12 hours of data, sampled every second
           seed = .seed, # set a seed for consistent results
           complete = TRUE, # prevents missing projection warning
           crs = PROJECTION) # projection CRS string
}
tels <- tibble(seed = 1:3, tel = map(seed, get_tracks)) %>%
  mutate(seed = factor(seed))

# plot intervals between encounters
intervals <-
  transmute(tels,
            seed,
            ints = map(tel, \(.tel) {
              .tel %>%
                data.frame() %>%
                mutate(cell_id = cellFromXY(object = HABITAT,
                                            xy = SpatialPoints.telemetry(.tel)),
                       new_cell = c(1, diff(cell_id)) != 0) %>%
                filter(new_cell) %>%
                mutate(interval = c(NA, diff(t)))
            })) %>%
  unnest(ints)

# check the main quantiles
quantile(intervals$interval, na.rm = TRUE)

# histogram of number of events as by the time interval between events
p_hist <-
  ggplot(intervals, aes(interval)) +
  facet_wrap(~ seed, ncol = 1) +
  geom_histogram(aes(fill = seed), binwidth = 30, center = 15, lwd = 0.1,
                 color = 'black', na.rm = TRUE) +
  scale_fill_viridis_d('Seed', ) +
  labs(x = 'Interval between encounters (s)', y = 'Count') +
  theme(panel.grid = element_blank())

# find number of encounters at each sampling time
d <- expand_grid(delta_t = seq(1, 300),
                 seed = 1:3) %>%
  mutate(data = map2(delta_t, seed, \(.dt, .seed) {
    .tel <- tels$tel[[.seed]]
    filtered <- .tel[.tel$t %% .dt == 0, ]
    filtered %>%
      data.frame() %>%
      mutate(cell_id = cellFromXY(object = HABITAT,
                                  xy = SpatialPoints.telemetry(filtered)),
             new_cell = c(1, diff(cell_id)) != 0)
  }),
  seed = factor(seed, levels = 1:3),
  encounters = map_int(data, \(.d) sum(.d$new_cell)),
  frac_encounters = map_dbl(data, \(.d) mean(.d$new_cell)))

p_events <-
  ggplot(d, aes(delta_t, encounters)) +
  facet_wrap(~ seed, ncol = 1) +
  geom_point(alpha = 0.3) +
  geom_smooth(aes(color = seed), method = 'gam', formula = y ~ s(x)) +
  scale_color_viridis_d('Seed') +
  labs(x = expression(paste('Sampling interval (', Delta, 't, s)')),
       y = 'Number of encounters') +
  scale_x_continuous(breaks = c(2^seq(0, 8, by = 2), 30), trans = 'log2') +
  theme(panel.grid = element_blank())

# check how many events are detected for each delta_t
p_prop <-
  d %>%
  select(delta_t, encounters, seed) %>% # only keep necessary columns
  # calculate the fraction of events lost by coarsening the sampling
  pivot_wider(names_from = delta_t, values_from = encounters) %>%
  mutate(max = `1`) %>%
  pivot_longer(-c(seed, max), names_to = 'delta_t', values_to ='encounters') %>%
  mutate(delta_t = as.numeric(delta_t), # from character back to numeric
         frac = encounters / max) %>% # calculate proportion of events detected
  ggplot(aes(delta_t, frac)) +
  facet_grid(seed ~ .) +
  geom_vline(xintercept = 30) +
  geom_point(alpha = 0.5) +
  geom_smooth(aes(color = seed), method = 'gam', formula = y ~ s(x)) +
  scale_color_viridis_d('Seed') +
  scale_x_continuous(trans = 'log2') + # to show the effect of doubling delta t
  labs(x = expression(paste('Sampling interval (', Delta, 't, s)')),
       y = 'Prop. of encounters detected    ') +
  theme(panel.grid = element_blank())

# plot a part of the track with multiple encounters in a short time interval
consecutive_encounters <- intervals %>%
  filter(seed == 1) %>% # only use track 1
  filter(interval <= 30) %>% # only keep pairs of events within 30 s
  filter(c(0, diff(t)) < 25 | c(diff(t), 0) < 25) %>% # keep groups within 25 s
  mutate(gross_t = round(t, -2)) %>% # round t to the hundreds
  group_by(gross_t) %>% # count how many encounters occurred in 100-s intervals
  mutate(count = n()) %>%
  ungroup() %>%
  filter(count == max(count)) # only keep interval with the most events

tel <- tels$tel[[1]]
tel_short <- tel[tel$t <= max(consecutive_encounters$t) + 500, ]

track <- tel_short %>%
  data.frame() %>%
  mutate(cell_id = cellFromXY(HABITAT, SpatialPoints.telemetry(tel_short)),
         new_cell = c(1, diff(cell_id)) != 0)

thin <- function(dt) {
  tel_short %>%
    data.frame() %>%
    filter(t %% dt == 0) %>%
    mutate(cell_id = cellFromXY(HABITAT,
                                tel_short[tel_short$t %% dt == 0, ] %>%
                                  SpatialPoints.telemetry()) %>%
             suppressWarnings(),
           new_cell = c(1, diff(cell_id)) != 0)
}

tracks <- tibble(delta_t = c(1, 5, 10, 20, 30, 60, 90, 120),
                 d = map(delta_t, thin),
                 group = factor(delta_t)) %>%
  mutate(group = factor(group, levels = group[order(delta_t)])) %>%
  unnest(d)

p_thinning <-
  ggplot(tracks, aes(x, y)) +
  facet_wrap(~ group,
             labeller = label_bquote(paste(Delta, 't = ',
                                           .(levels(tracks$group)[group]), s)),
                                     nrow = 2) +
  coord_equal(xlim = c(-0.1, 0.4), ylim = c(-0.5, 0.1)) +
  geom_tile(data = as.data.frame(HABITAT, xy = TRUE) %>%
              filter(abs(x) < 5, abs(y) < 5), fill = 'transparent', color = 'black') +
  geom_path(color = '#440154', lwd = 1) +
  geom_point(aes(color = new_cell, size = new_cell, shape = new_cell)) +
  geom_point(aes(color = new_cell, size = new_cell, shape = new_cell),
             filter(tracks, new_cell)) +
  scale_x_continuous(name = NULL, breaks = NULL) +
  scale_y_continuous(name = NULL, breaks = NULL) +
  scale_color_brewer('New cell', labels = c('No', 'Yes'), type = 'qual', palette = 6) +
  scale_size_manual('New cell', labels = c('No', 'Yes'), values = c(1, 2)) +
  scale_shape_manual('New cell', labels = c('No', 'Yes'), values = c(19, 15)) +
  theme_bw() +
  theme(legend.position = 'none')

full_tracks <- tels %>%
  mutate(tel = map(tel, data.frame)) %>%
  unnest(tel)
full_tracks <- mutate(full_tracks,
                      cell_id = cellFromXY(object = HABITAT,
                                           xy = SpatialPoints.telemetry(tels$tel)),
                      new_cell = c(1, diff(cell_id)) != 0)
p_tracks <-
  ggplot(full_tracks, aes(x, y)) +
  coord_equal(xlim = range(full_tracks$x), ylim = range(full_tracks$y)) +
  geom_tile(data = as.data.frame(HABITAT, xy = TRUE) %>%
              filter(x >= min(full_tracks$x) - 20, x <= max(full_tracks$x) + 20,
                     y >= min(full_tracks$y) - 20, y <= max(full_tracks$y) + 20),
            fill = 'transparent', color = '#00000020') + # black with some alpha
  geom_path(aes(group = seed), color = 'grey', lwd = 1) +
  geom_path(aes(group = seed, color = seed)) +
  scale_color_viridis_d() +
  scale_x_continuous(NULL, breaks = NULL) +
  scale_y_continuous('', breaks = NULL) +
  theme_bw() +
  theme(legend.position = 'none', panel.grid = element_blank())

# check that feeding events occur when moving to a new cell
p_tracks +
  geom_point(aes(alpha = new_cell), size = 1) +
  scale_alpha_manual(values = c(0, 1))

# create the final figure
panels <-
  plot_grid(plot_grid(get_legend(p_hist + theme(legend.position = 'top')),
                      get_legend(p_thinning + theme(legend.position = 'top'))),
            plot_grid(p_hist, p_prop, p_thinning, p_tracks,
                      labels = c('a.', 'b.', 'c.', 'd.'), ncol = 2),
            ncol = 1, rel_heights = c(1, 10))

ggsave('figures/thinning-examples.png', plot = panels,
       width = 7, height = 5, scale = 1.2, bg = 'white')
