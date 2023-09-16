library('ctmm')      # for tracking data and movement modeling
library('dplyr')     # for data wrangling
library('ggplot2')   # for fancy plots
library('mgcv')      # for GAMs
library('lubridate') # for smoother date wrangling
library('sf')        # for spatial objects

source('functions/betals.r') # beta location-scale fam 
source('analysis/figures/default-figure-styling.R') # for NDVI color palette

# import 95% home range estimates with no CIs (optional)
uds <- readRDS('../tapirs/models/tapirs-akdes.rds') %>%
  filter(region == 'cerrado') %>%
  pull(akde) %>%
  purrr::map_dfr(\(x) {
    SpatialPolygonsDataFrame.UD(x, level.UD = 0.95, level = 0) %>%
      sf::st_as_sf() %>%
      sf::st_transform(crs = '+proj=longlat') %>%
      st_union() %>%
      st_as_sf()
  })
uds_merged <- uds %>%
  st_union() %>%
  st_buffer(3e3)
plot(uds_merged, col = 'grey')
plot(uds, add = TRUE, col = 'forestgreen')

# import NDVI data
d <- readRDS('data/ndvi-rasters/tapirs/tapirs-data.rds') %>%
  mutate(dec_date = decimal_date(date))

d %>%
  filter(dec_date %in% unique(dec_date)[1:6]) %>%
  ggplot(aes(long, lat, fill = ndvi)) +
  facet_wrap(~ round(dec_date, 2)) +
  geom_raster() +
  geom_sf(dat = uds, inherit.aes = FALSE, alpha = 0.3) +
  scale_fill_gradientn('NDVI', colours = ndvi_pal, limits = c(-1, 1))

range(d$ndvi)
n_distinct(d$dec_date)

# run a test model
m_test <-
  bam(
    (ndvi + 1) / 2 ~ # not scaling because range is in (0, 1)
      s(long, lat, bs = 'ds', k = 20) + # mean over space
      s(dec_date, bs = 'tp', k = 10), # need high k to account for animal adapting
    family = betar(link = 'logit'),
    data = d,
    method = 'fREML',
    discrete = TRUE,
    control = gam.control(trace = TRUE))
plot(m_test, pages = 1, n = 400, scheme = 3)
summary(m_test)

m_ndvi <-
  gam(list(
    # mean predictor
    ndvi ~ # not scaling because range is in (0, 1)
      s(long, lat, bs = 'ds', k = 50) + # mean over space
      s(dec_date, bs = 'tp', k = 10), # need high k to account for animal adapting
    # scale predictor (sigma2 = mu * (1 - mu) * scale)
    ~
      s(long, lat, bs = 'ds', k = 25) +
      s(dec_date, bs = 'tp', k = 10)),
    family = betals(),
    data = d,
    method = 'REML',
    control = gam.control(trace = TRUE))

saveRDS(m_ndvi, file = 'models/tapirs/tapirs-ndvi-betals.rds')

if(FALSE) {
  plot(m_ndvi, pages = 1, scheme = 3) # plot smooths
}
