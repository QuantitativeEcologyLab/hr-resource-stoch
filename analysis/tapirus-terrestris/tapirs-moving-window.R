library('dplyr')   # for data wrangling
library('purrr')   # for functional programming
library('tidyr')   # for data wrangling
library('ggplot2') # for fancy plots
library('cowplot') # for fancy multi-panel plots
library('ctmm')    # for movement modeling
source('functions/window_hr.R') # function to calculate HRs and create figures
theme_set(theme_bw()) # change default theme

# import tapir data from https://github.com/StefanoMezzini/tapirs
tapirs <- readRDS('../tapirs/models/tapirs-final.rds') %>%
  select(region, name, name.short, data) %>%
  filter(tapirs, region == 'cerrado')

# check ranges of tracking data
tapirs <-
  tapirs %>%
  mutate(start = map_chr(data, \(d) min(d$timestamp) %>% as.character()),
         end = map_chr(data, \(d) max(d$timestamp) %>% as.character()),
         days = map_dbl(data, \(d) diff(range(d$timestamp))),
         clicks = map_dbl(data, nrow),
         daily_clicks = clicks / days)%>%
  arrange(days, clicks, desc(daily_clicks))

if(FALSE) {
  View(tapirs)
  
  # plot the tapir data
  tapirs %>%
    mutate(data = map(data, data.frame)) %>% # convert to df for unnesting
    unnest(data) %>%
    slice(seq(1, n(), by = 10)) %>%
    ggplot() +
    facet_wrap(~ name, scales = 'free') +
    geom_path(aes(longitude, latitude), alpha = 0.05) +
    geom_point(aes(longitude, latitude, color = timestamp), alpha = 0.1) +
    scale_color_viridis_c('Date') +
    theme(legend.position = c(x = 0.7, y = 0.1),
          legend.direction = 'horizontal',
          legend.key.width = unit(0.5, 'in'))
}

# find some good examples ----
# Cerrado ----
ce_proj <- '+proj=utm +zone=22 +datum=NAD83 +units=m'

cerrado <- filter(tapirs, region == 'cerrado')
select(cerrado, name, days, clicks, daily_clicks)

cerrado %>%
  filter(name != 'CE_08_COLOMBINA') %>% # remove tapir that fails 
  # filter(name.short %in% # filter to only a few good examples
  #          c('ANNA', 'DONALINA', 'KURUKA', 'LOU', 'ZEFA', 'ZACA', 'TITI', 'ZIGGY')) %>%
  pull(data) %>% # extract data column
  map(function(.d) {
    projection(.d) <- ce_proj # add projection based on mean long and lat
    cat('Running', .d@info$identity, '\n')
    .d <- window_hr(.d,
                    window = 15 %#% 'day', # 1 week of data for sufficient sample size
                    dt = 15 %#% 'day', # move window over by a single day each time
                    fig_path = 'figures/tapirs',
                    rds_path = 'models/tapirs',
                    cores = 1) # cannot parallelize on Windows
  })

# join all the Cerrado data into a single tibble
cerrado <-
  list.files(path = 'models/tapirs',
             pattern = '*-window-15-days-dt-15-days.rds',
             full.names = TRUE) %>%
  tibble(name = .,
         d = map(., \(x) readRDS(x))) %>%
  mutate(name = substr(name,
                       nchar('models/tapirs/CE_xx_x'),
                       nchar(name) -
                         nchar('-window-15-days-dt-15-days.rds'))) %>%
  unnest(cols = d) %>%
  select(name, t_center, hr_est_95, hr_lwr_95, hr_upr_95) %>%
  filter(! is.na(hr_est_95)) %>%
  mutate(t_center = as.POSIXct(t_center)) %>%
  group_by(name) %>%
  # add weights based on CI width
  mutate(weight = 1 / (hr_upr_95 - hr_lwr_95),
         weight = weight / mean(weight)) %>% # keep sum(weight) == n()
  ungroup()

# ensure weights are correct
sum(cerrado$weight) == nrow(cerrado)

saveRDS(cerrado, 'data/cerrado-hrs.rds')

# Pantanal ----
pantanal <- filter(tapirs, region == 'pantanal')
select(pantanal, name, days, clicks, daily_clicks) %>%
  arrange(desc(days))

# fails: PA_26_SACHIN, PA_16_DORA 
pantanal %>%
  filter(name %in% # filter to only a few good examples
           c('PA_08_BENJAMIN', 'PA_10_VIVEK', 'PA_30_NELSAO',
             'PA_31_JORDANO', 'PA_23_RICK', 'PA_24_LESLIE', 'PA_12_FELIPPE',
             'PA_21_CAIO')) %>%
  pull(data) %>% # extract data column
  map(function(.d) {
    cat('Running', .d@info$identity, '\n')
    .d <- window_hr(.d,
                    window = 7 %#% 'day', # 1 week of data for sufficient sample size
                    dt = 1 %#% 'day', # move window over by a single day each time
                    fig_path = 'figures/tapirs',
                    rds_path = 'models/tapirs',
                    cores = 1) # cannot parallelize on Windows
  })
