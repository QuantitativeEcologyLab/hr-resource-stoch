library('ctmm')      # for tracking data and movement modeling
library('dplyr')     # for data wrangling
library('ggplot2')   # for fancy plots
library('mgcv')      # for GAMs
library('gratia')    # for ggplot-based GAM graphics
library('lubridate') # for smoother date wrangling

source('functions/betals.r') # beta location-scale fam 
source('analysis/figures/default-figure-styling.R') # for NDVI color palette

# import 95% home range estimates with no CIs (optional)
ud <- readRDS('../tapirs/models/tapirs-akdes.rds') %>%
  filter(name.short == 'ANNA') %>%
  pull(akde) %>%
  purrr::map_dfr(\(x) {
    SpatialPolygonsDataFrame.UD(x, level.UD = 0.95, level = 0) %>%
      sf::st_as_sf() %>%
      sf::st_transform(crs = '+proj=longlat')
    })

# import NDVI data
d <- readRDS('data/ndvi-rasters/tapir-anna/tapir-anna-data.rds') %>%
  mutate(dec_date = decimal_date(date)) %>%
  # remove unusual raster (likely due to a flood)
  filter(date != '2017-12-19')

ggplot(d, aes(long, lat, fill = ndvi)) +
  facet_wrap(~ round(dec_date, 2)) +
  geom_raster() +
  geom_sf(dat = ud, inherit.aes = FALSE, alpha = 0.3) +
  scale_fill_gradientn('NDVI', colours = ndvi_pal, limits = c(-1, 1))

m_ndvi <-
  gam(list(
    # mean predictor
    ndvi ~ # not scaling because range is in (0, 1)
      s(long, lat, bs = 'ds', k = 50) + # mean over space
      s(dec_date, bs = 'tp', k = 10), # mean over time
    # scale predictor (sigma2 = mu * (1 - mu) * scale)
    ~
      s(long, lat, bs = 'ds', k = 30) +
      s(dec_date, bs = 'tp', k = 10)),
    family = betals(),
    data = d,
    method = 'REML')

draw(m_ndvi, n = 250, rug = FALSE) # plot smooths
saveRDS(m_ndvi, file = 'models/tapirs/CE_31_ANNA-mgcv-ndvi-betals.rds')
