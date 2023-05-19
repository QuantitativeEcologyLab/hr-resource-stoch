library('ctmm')      # for tracking data and movement modeling
library('tidyr')     # for data wrangling
library('dplyr')     # for data wrangling
library('ggplot2')   # for fancy plots
library('mgcv')      # for GAMs
library('lubridate') # for smoother date wrangling

source('analysis/figures/default-figure-styling.R') # for NDVI color palette
source('functions/window_hr.R') # function to calculate HRs and create figures

# import movement model and telemetry
load('C:/Users/mezzinis/Dropbox/Empirical_Results/Aepyceros_melampus/Fits/Fits_1757.rda')
load('C:/Users/mezzinis/Dropbox/Uni/tels/Aepyceros_melampus.Rda')
tel <- Aepyceros_melampus %>%
  filter(individual.local.identifier == '1757') %>%
  as.telemetry()
plot(tel)
diff(range(tel$timestamp))

# import 95% home range estimates with no CIs (optional)
ud <- akde(tel, FIT)
plot(tel)
plot(ud, add = TRUE)

# fit moving window
if(! file.exists('models/aepyceros-melampus-1757-window-7-days-dt-7-days.rds')){
  tel@info$identity <- 'aepyceros-melampus-1757'
  window_hr(tel,
            window = 7 %#% 'day', # 1 week of data for sufficient sample size
            dt = 1 %#% 'day', # move window over by a single day each time
            movement_model = FIT,
            fig_path = 'figures',
            rds_path = 'models',
            cores = 1)
}

mw <- readRDS('models/aepyceros-melampus-1757-window-7-days-dt-7-days.rds')
hist(mw$hr_est_95)
abline(v = summary(FIT, units = FALSE)$CI[1, 2] / 1e3^2, col = 'red')
plot(hr_est_95 ~ t_center, mw)

ud_poly <- SpatialPolygonsDataFrame.UD(ud) %>%
  sf::st_as_sf() %>%
  sf::st_transform('+proj=longlat')

# import NDVI data
d <- readRDS('data/ndvi-rasters/Aepyceros_melampus/1757/Aepyceros_melampus-1757-ndvi.rds') %>%
  mutate(dec_date = decimal_date(date))

ggplot(d, aes(long, lat, fill = ndvi)) +
  facet_wrap(~ round(dec_date, 2), nrow = 3) +
  geom_raster() +
  geom_sf(data = ud_poly, inherit.aes = FALSE, alpha = 0.3) +
  scale_fill_gradientn('NDVI', colours = ndvi_pal, limits = c(-1, 1))

if(file.exists('models/aepyceros-melampus/1757-ndvi-gaulss.rds')) {
  m_ndvi <- readRDS(file = 'models/aepyceros-melampus/1757-ndvi-gaulss.rds')
} else {
  m_ndvi <-
    gam(list(
      # mean predictor
      ndvi ~
        s(long, lat, bs = 'ds', k = 20) + # mean over space
        s(dec_date, bs = 'tp', k = 20), # need high k to account for animal adapting
      # precision (1/standard deviation) predictor
      ~
        s(long, lat, bs = 'ds', k = 10) +
        s(dec_date, bs = 'tp', k = 10)),
      family = gaulss(b = 0.0001), # using minimum standard deviation of 0.0001
      data = d,
      method = 'REML',
      control = gam.control(nthreads = 4, trace = TRUE))
  
  saveRDS(m_ndvi, file = 'models/aepyceros-melampus/1757-ndvi-gaulss.rds')
  
  if(FALSE) {
    plot(m_ndvi, pages = 1, scheme = 3, scale = 0, n = 250) # plot smooths
  }
}
