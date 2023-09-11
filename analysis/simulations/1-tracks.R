# does not require large amounts of RAM, takes < 1 minute
library('ctmm')   # for generating movement models
library('terra') # for working with raster data
library('dplyr')  # for data wrangling
library('purrr')  # for functional programming
source('functions/label_visits.R') # function for when animal moves to a new cell & eats
source('analysis/simulations/movement-model.R') # for consistency between scripts

N_DAYS <- 2^10 # number of "days" (i.e., simulations with different seeds)
SIMULATE_TRACKS <- TRUE # set to "TRUE" if you want to generate new tracks

# extracts simulated tracks from a ctmm movement model for given sample times
get_tracks <- function(day, times = SAMPLES) {
  simulate(model, # ctmm movement model
           t = times, # sampling times in seconds
           seed = day, # set a seed for consistent results (different track each day)
           complete = TRUE, # add lat, long, and timestamp to the telemetry object
           crs = PROJECTION) # CRS projection string
}

# generate simulated tracks (will be truncated at satiety later)
if(SIMULATE_TRACKS) {
  tracks <- tibble(day = 1:N_DAYS, # a simulation for each day
                   tel = map(.x = day, # set a seed for consistent results
                             .f = get_tracks)) # function to generate tracks
  saveRDS(tracks, file = 'simulations/tracks.rds')
} else tracks <- readRDS('simulations/tracks.rds')

# find patch visits and calories consumed from the tracks
tracks <- transmute(tracks, # drop tel column
                    day, # keep day column
                    track = map(.x = tel, # add a column of full tracks
                                .f = \(x) label_visits(.tel = x, .habitat = HABITAT)))

# make a single, large tibble
tracks <- tidyr::unnest(tracks, track)
saveRDS(tracks, file = 'simulations/labelled-tracks.rds')

# check if grid any animals moved out of the habitat raster
if(any(is.na(tracks$new_cell))) {
  warning(paste('CAUTION:', sum(tracks$new_cell), 'animals left the habitat raster!'))
}

# check scaling is consistent and sampling frequency is sufficiently high
if(FALSE) {
  library('ggplot2')
  theme_set(theme_bw())
  
  test <- readRDS('simulations/tracks.rds')
  
  track <- test$tel[[2]] %>%
    data.frame() %>%
    mutate(cell_id = cellFromXY(object = HABITAT,
                                xy = SpatialPoints.telemetry(test$tel[[2]]) %>%
                                  as.data.frame() %>%
                                  select(x, y)),
           new_cell = c(1, diff(cell_id)) != 0)
  mean(track$new_cell) # fraction of observations that are feeding events
  
  track_short <- filter(track, t < 1e4)
  
  # convert the raster to a data.frame to ensure tracks don't go too far
  test_habitat <- as.data.frame(HABITAT, xy = TRUE)
  test_habitat_small <- filter(test_habitat,
                               x >= min(track$x), x <= max(track$x),
                               y >= min(track$y), y <= max(track$y))
  
  # plot a portion of a track
  ggplot(track_short, aes(x, y)) +
    coord_equal(xlim = range(track_short$x), ylim = range(track_short$y)) +
    geom_tile(data = test_habitat_small, fill = 'transparent', color = 'black') +
    geom_path() +
    geom_point(aes(color = new_cell, size = new_cell)) +
    geom_point(aes(0, 0), color = 'darkorange', size = 4) +
    scale_color_brewer('New cell', type = 'qual', palette = 6) +
    scale_size_manual('New cell', values = c(1, 2))
  
  # plot an entire track
  ggplot(track, aes(x, y)) +
    coord_equal(xlim = range(track$x), ylim = range(track$y)) +
    geom_tile(data = test_habitat, fill = 'transparent', color = 'black') +
    geom_path() +
    geom_point(aes(color = new_cell, size = new_cell)) +
    geom_point(aes(0, 0), color = 'darkorange', size = 4) +
    scale_color_brewer('New cell', type = 'qual', palette = 6) +
    scale_size_manual('New cell', values = c(1, 2))
  
  # create a tibble of all the tracks
  test <- test %>%
    mutate(tel = map(tel, data.frame)) %>%
    tidyr::unnest(tel) %>%
    group_by(day) %>%
    mutate(new_cell = cellFromXY(object = HABITAT,
                                 xy = SpatialPoints.telemetry(tibble(x, y)) %>%
                                   as.data.frame() %>%
                                   select(x, y)) %>%
             suppressWarnings(),
           new_cell = c(0, diff(new_cell)) != 0)
  
  # plot all the tracks (can take some time with 1024 tracks...)
  test %>%
    ggplot(aes(x, y)) +
    coord_equal(xlim = range(test$x), ylim = range(test$y)) +
    geom_tile(data = test_habitat, fill = 'transparent', color = 'black') +
    geom_path(aes(group = day), alpha = 0.4)
  
  # make a histogram of the total number of encounters per track
  test %>%
    group_by(day) %>%
    summarize(encounters = sum(new_cell)) %>%
    ggplot() +
    geom_histogram(aes(encounters), bins = ceiling(log2(N_DAYS)) + 1,
                   color = 'black')
  
  test %>%
    group_by(day) %>%
    summarize(encounters = sum(new_cell)) %>%
    pull(encounters) %>%
    range()
}
