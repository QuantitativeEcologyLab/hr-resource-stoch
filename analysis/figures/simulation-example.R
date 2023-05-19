library('ctmm')    # for movement modeling
library('dplyr')   # for data wrangling (%>%, mutate())
library('purrr')   # for functional programming (e.g., map*())
library('tidyr')   # for data wrangling (e.g., unnest())
library('ggplot2') # for plotting
library('cowplot') # for plots in grids
source('functions/label_visits.R') # label food encounters based on a raster
source('functions/pull_parameter.R') # to extract parameters from ctmm and akde objects
source('analysis/figures/default-figure-styling.R') # defaults for figures (theme, size)
theme_set(theme_map() + theme(legend.position = 'none'))
select <- dplyr::select # don't use raster::select as a default

REQUIRED <- 8 # amount of food required for satiety
DELTA_T <- 0.1 # time between measurements
SAMPLES <- seq(0, 30, by = DELTA_T) # sampling times
DIM <- 15 # dimensions of raster (number of cells per side)
HABITAT <- matrix(data = 1, nrow = DIM, ncol = DIM) %>% # raster of patches
  raster(xmx = 20, xmn = -20, ymx = 20, ymn = -20)
HABITAT_tbl <- rasterToPoints(HABITAT) %>% as_tibble() # raster in tibble format

# resource abundance palette
LOW <- 'white'
HIGH <- '#AA3377'

m <- ctmm(tau = c(Inf, 1), sigma = 1, mu = c(0, 0)) # infinitely diffusive model
tracks <-
  tibble(day = c(-10:-8, -6, -4:-2, 0:5), # -7, 7, 8, 9 give very linear movement
         tel = map(day, \(i) {
           simulate(m, # ctmm movement model
                    t = SAMPLES, # sampling times in seconds
                    seed = i, # consistent but different track each day
                    complete = TRUE) # add lat, long, and timestamp to track
         })) %>%
  mutate(labelled = map(tel,\(x) label_visits(.tel = x, .habitat = HABITAT)),
         truncated = map(labelled, \(x) x %>%
                           # animal eats whenever it crosses to a new cell
                           mutate(food = rgamma(n = n(), shape = 2, rate = 1),
                                  food = if_else(new_cell, food, 0),
                                  satiety = cumsum(food),
                                  full = satiety > REQUIRED) %>%
                           # only keep "not full" and first "full"
                           filter((! full) | (! duplicated(full))) %>%
                           mutate(t_full = t[which(full & ! duplicated(full))])),
         t_max = map_dbl(truncated, \(x) max(x$t)),
         t_start = lag(t_max * 2), # allow time for return
         t_start = if_else(is.na(t_start), 0, t_start))

# ensure the animal got full each day
transmute(tracks, day,
          full = map_lgl(truncated, \(x) any(x$full, na.rm = TRUE)))

# dataset of full tracks not truncated
labelled <- tracks %>%
  select(day, labelled) %>%
  unnest(labelled)

tracks <-
  tracks %>%
  select(- tel, - labelled) %>%
  mutate(t_start = cumsum(t_start)) %>%
  unnest(truncated) %>% # ungroup the data for each day
  mutate(t = t + t_start) %>% # change `t` so days occur consecutively
  # change timestamp to reflect the new t, add a dummy identifier for the animal
  mutate(timestamp = as.POSIXct(t, tz='UTC', origin = as.POSIXct('1970-01-01 00:00 UTC')),
         individual.local.identifier = 'example_animal') %>%
  relocate(c(t_start, t), .after = day) # time columns together

# ensure time occurs sequentially
ggplot(tracks, aes(timestamp, t, color = factor(day))) +
  geom_line() +
  theme_bw()

# example day
d <- filter(labelled, day == -9)
track <- filter(tracks, day == -9)

# animal leaves to explore each day and eat at points (starts at square)
ggplot(mapping = aes(x, y)) +
  coord_equal() +
  geom_path(data = track) +
  geom_point(data = filter(track, new_cell), size = 2) +
  geom_point(aes(0, 0), pch = 18, size = 4) +
  geom_point(aes(0, 0), pch = 18, size = 2, color = 'yellow')

# plot all tracks
ggplot(tracks) +
  coord_equal() +
  geom_path(aes(x, y, color = factor(day))) +
  geom_point(aes(x, y, color = factor(day))) +
  theme(legend.position = 'right')

# estimate home range
tel <- as.telemetry(tracks) # change data to telemetry structure
theta <- ctmm.guess(data = tel, interactive = FALSE) # variogram
model <- ctmm.fit(tel, theta) # fit movement model (OUF)
AKDE <- akde(tel, model) # HR autocorrelated kernel density estimate

# 95% quantile autocorrelated kernel density estimates
akde_95 <-
  SpatialPolygonsDataFrame.UD(AKDE, level.UD = 0.95) %>%
  fortify() %>%
  as_tibble()

# create the figure ----
HABITAT_tbl <-
  rasterToPoints(HABITAT) %>% # convert raster to a data.frame
  as_tibble() %>% # convert to a tibble for ease of use
  filter(x >= min(akde_95$long) - 1, x <= max(akde_95$long),
         y >= min(akde_95$lat) - 1, y <= max(akde_95$lat) + 1)

# a) movement from the first "day"
p_a <-
  ggplot(d, aes(x, y)) +
  coord_equal() +
  geom_tile(data = HABITAT_tbl, fill = 'transparent', color = 'black')+
  geom_path(lwd = 1, linejoin = 'round', color = 'grey') +
  geom_point(aes(0, 0), pch = 18, size = 4) +
  geom_point(aes(0, 0), pch = 18, size = 2, color = 'yellow')

# b) encounters (green dots)
p_b <-
  p_a +
  geom_point(data = filter(d, new_cell), size = 1, color = 'darkgreen') +
  geom_point(aes(0, 0), pch = 18, size = 4) +
  geom_point(aes(0, 0), pch = 18, size = 2, color = 'yellow')

# c) truncated path with satiety
p_c <-
  ggplot(track, aes(x, y)) +
  coord_equal() +
  geom_tile(data = HABITAT_tbl, fill = 'transparent', color = 'black') +
  geom_path(aes(color = satiety, group = day), lwd = 1, lineend ='round') +
  geom_point(data = filter(track, new_cell), size = 1, color = 'darkgreen') +
  geom_point(aes(0, 0), pch = 18, size = 4) +
  geom_point(aes(0, 0), pch = 18, size = 2, color = 'yellow') +
  scale_color_gradient('Satiety', low = LOW, high = HIGH,
                       limits = c(0, REQUIRED), breaks = c(0, REQUIRED),
                       labels = c('Low', 'High'))

# d) all tracks with satiety
p_d <-
  ggplot(tracks, aes(x, y)) +
  coord_equal() +
  geom_tile(data = HABITAT_tbl, fill = 'transparent', color = 'black') +
  geom_path(aes(color = satiety, group = day), lwd = 1, lineend='round') +
  geom_point(aes(0, 0), pch = 18, size = 4) +
  geom_point(aes(0, 0), pch = 18, size = 2, color = 'yellow') +
  scale_color_gradient('Satiety', low = LOW, high = HIGH,
                       limits = c(0, REQUIRED), breaks = c(0, REQUIRED),
                       labels = c('Low', 'High'))

# e) HR estimate
p_e <-
  ggplot(tracks, aes(x, y)) +
  coord_equal() +
  geom_tile(data = HABITAT_tbl, fill = 'transparent', color = 'black') +
  geom_polygon(aes(long, lat), color = 'black', fill = pal[3], alpha = 0.3,
               filter(akde_95, id == 'example_animal 95% est'))+
  geom_path(aes(color = satiety, group = day), lwd = 1, lineend ='round') +
  geom_polygon(aes(long, lat), color = 'black', fill = 'transparent',
               filter(akde_95, id == 'example_animal 95% est'))+
  geom_point(aes(0, 0), pch = 18, size = 4) +
  geom_point(aes(0, 0), pch = 18, size = 2, color = 'yellow') +
  scale_color_gradient('Satiety', low = LOW, high = HIGH,
                       limits = c(0, REQUIRED), breaks = c(0, REQUIRED),
                       labels = c('Low', 'High'))

p <-
  plot_grid(p_a, p_b, p_c, p_d, p_e,
          get_legend(
            p_e +
              scale_color_gradient('Satiety', low = LOW, high = HIGH,
                                   limits = 0:1, breaks = 0:1,
                                   labels = c('Empty', 'Full')) +
              theme_bw() +
              theme(legend.key.height = unit(0.5, 'cm'),
                    text = element_text(face = 'bold', size = 18),
                    legend.key = element_rect(color = 'black'))),
          nrow = 1, labels = c('a.', 'b.', 'c.', 'd.', 'e.', ''),
          rel_widths = c(rep(1, 5), 0.4), label_size = 22)

ggsave('figures/simulation-example.png', plot = p, width = 20,
       height = 3, scale = 1, bg = 'white', dpi = 600)
