# see https://rspatialdata.github.io/vegetation.html for example code and more info
#' install `rgeoboundaries` with `remotes::install_github('wmgeolab/rgeoboundaries')`
#' install `MODIStsp` with `remotes::install_github('ropensci/MODIStsp')`
#' `MODIStsp` version >= 2.0.7 from github fixes login issue due to 
library('dplyr')    # for data wrangling
library('sf')       # for spatial features
library('MODIStsp') # for downloading NDVI rasters
library('purrr')  # for functional programming
library('tidyr')  # for data wrangling
library('terra') # to import and save rasters
source('earthdata-login-info.R') # import personal login credentials for EarthData

# use a bounding box from tracking data ----
library('ctmm') # for movement data and models
load('C:/Users/mezzinis/Dropbox/Uni/tels/Aepyceros_melampus.Rda')
tel <- Aepyceros_melampus %>%
  filter(individual.local.identifier == '1757') %>%
  as.telemetry()

# import movement model and fit UD
load('H:/Empirical_Results/Aepyceros_melampus/Fits/Fits_1757.rda')
ud <- akde(data = tel, CTMM = FIT, weights = FALSE)
plot(ud, level.UD = 0.995)
plot(tel, add = TRUE, cex = 1, error = FALSE) # UD is fairly wide

bbox <-
  SpatialPolygonsDataFrame.UD(ud, level.UD = 0.9995, level = 0) %>%
  st_as_sf() %>%
  st_transform(crs = '+proj=longlat') %>%
  st_bbox()

# download NDVI
MODIStsp(gui = FALSE, # do not use the browser GUI, only run in R
         out_folder = 'data/ndvi-rasters/Aepyceros_melampus/1757/',
         selprod = 'Vegetation Indexes_16Days_250m (M*D13Q1)', # can't specify Terra here
         prod_version = '061', # 2022 raster version
         bandsel = 'NDVI', # Normalized Difference Vegetation Index layer only
         sensor = 'Terra', # only terrestrial values, ignore main bodies of water
         user = USERNAME, # your Earthdata username (for urs.earthdata.nasa.gov/home)
         password = PASSWORD, # your Earthdata password
         start_date = format(min(tel$timestamp) - 16, '%Y.%m.%d'), # 1 raster before min
         end_date = format(max(tel$timestamp) + 16, '%Y.%m.%d'), # 1 raster after end
         spatmeth = 'bbox', # use a bounding box for the extent
         bbox = bbox, # spatial file for raster extent
         out_projsel = 'User Defined', # use specified projection instead of default
         output_proj = '+proj=longlat', # download unprojected raster
         resampling = 'bilinear', # method for resampling raster if changing projection
         delete_hdf = TRUE, # delete HDF files after download is complete
         scale_val = TRUE, # convert from integers to floats within [-1, 1]
         out_format = 'GTiff', # output format: 'ENVI' (.hdr) or 'GTiff' (.tif)
         n_retries = 10, # number of times to try again if download fails before aborting
         verbose = TRUE, # print processing messages
         parallel = 4) # use TRUE for automatic number of cores (max 8), or specify

# check rasters ----
ud_poly <- SpatialPolygonsDataFrame.UD(ud, level.UD = 0.9995, level = 0) %>%
  st_as_sf() %>%
  st_transform('+proj=longlat')

rasters <-
  list.files(path = 'data/ndvi-rasters/Aepyceros_melampus/1757/VI_16Days_250m_v61/NDVI',
             pattern = '.tif', full.names = TRUE) %>%
  rast()

names(rasters) <-
  list.files(path = 'data/ndvi-rasters/Aepyceros_melampus/1757/VI_16Days_250m_v61/NDVI/',
             pattern = '.tif', full.names = FALSE) %>%
  substr(start = nchar('MOD13Q1_NDVI_X'), stop = nchar('MOD13Q1_NDVI_XXXX_XXX')) %>%
  as.Date(format = '%Y_%j')

if(FALSE) {
  raster::plot(rasters[[1]])
  plot(ud_poly, add = TRUE, color = 'transparent')
  tel %>%
    SpatialPoints.telemetry() %>%
    st_as_sf() %>%
    st_transform(crs = '+proj=longlat') %>%
    plot(add = TRUE)
  
  rasters %>%
    mask(ud_poly) %>%
    crop(ud_poly) %>%
    raster::plot()
}

# save NDVI data as an rds file of a tibble
rasters %>%
  as.data.frame(xy = TRUE) %>%
  pivot_longer(-c(x, y)) %>%
  transmute(long = x,
            lat = y,
            date = as.Date(name),
            ndvi = value) %>%
  saveRDS('data/ndvi-rasters/Aepyceros_melampus/1757/Aepyceros_melampus-1757-ndvi.rds')
