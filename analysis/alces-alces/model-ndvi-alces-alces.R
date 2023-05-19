library('ctmm')      # for tracking data and movement modeling
library('dplyr')     # for data wrangling
library('ggplot2')   # for fancy plots
library('mgcv')      # for GAMs
library('lubridate') # for smoother date wrangling

# for NDVI color palette
source('analysis/figures/barga-italy-2023/default-figure-styling-barga-2023.R')

# import NDVI data
d <- readRDS('data/ndvi-rasters/alces-alces-nw/alces-alces-ndvi.rds') %>%
  mutate(dec_date = decimal_date(date),
         year = year(date),
         doy = yday(date))

#' fit gaussian GAM with `bam()` to check estimated degrees of freedom
summarize(d,
          years = n_distinct(year),
          days = n_distinct(doy),
          locations = n() / years / days)

# takes a while to initialize (~20 seconds), but fits very fast (~10 seconds)
m_ndvi_0 <-
  bam(ndvi ~
        s(long, lat, bs = 'ds', k = 50) + # mean over space
        s(year, bs = 'tp', k = 9) +       # mean between years
        s(doy, bs = 'tp', k = 10),        # mean within years
      family = gaussian(),
      data = d,
      method = 'fREML',
      discrete = TRUE,
      control = gam.control(nthreads = 8, trace = TRUE))

# compare edf to k for each smooth
summary(m_ndvi_0)

layout(matrix(c(1, 1, 2, 3), ncol = 2))
plot(m_ndvi_0, scheme = 3)
layout(1)

m_ndvi <-
  gam(list(
    # mean predictor
    ndvi ~
      s(long, lat, bs = 'ds', k = 50) +         # mean over space
      s(year, bs = 'cr', k = 9) +                # mean between years
      s(doy, bs = 'cc', k = 10) +                # mean within years
      ti(year, doy, bs = c('cr', 'cc'), k = 5) + # s(doy) changes between years
      ti(long, lat, doy, bs = c('ds', 'cc'), k = c(10, 5), d = c(2, 1)),
    # precision (1/standard deviation) predictor
    ~
      s(long, lat, bs = 'ds', k = 25) +
      s(year, bs = 'cr', k = 4) +
      s(doy, bs = 'cc', k = 10) +
      ti(year, doy, bs = c('cr', 'cc'), k = 5) +
      ti(long, lat, doy, bs = c('ds', 'cc'), k = c(10, 5), d = c(2, 1))),
    family = gaulss(),
    data = d,
    method = 'REML',
    # optimizer = 'bfgs', # checking if this is faster
    knots = list(doy = c(0.5, 366.5)),
    control = gam.control(nthreads = 8, trace = TRUE))

if(FALSE) {
  summary(m_ndvi)
  plot(m_ndvi, pages = 1, scheme = 3, scale = 0, n = 250) # plot smooths
}

saveRDS(m_ndvi, file = 'models/alces-alces-nw-mgcv-ndvi-gaulss.rds')
