# see https://rspatialdata.github.io/vegetation.html for example code and more info
#' install `rgeoboundaries` with `remotes::install_github('wmgeolab/rgeoboundaries')`
#' install `MODIStsp` with `remotes::install_github('ropensci/MODIStsp')`
#' `MODIStsp` version >= 2.0.7 from github fixes login issue due to 
library('dplyr')    # for data wrangling
library('sf')       # for spatial features
library('MODIStsp') # for downloading NDVI rasters
library('purrr')    # for functional programming
library('tidyr')    # for data wrangling
library('terra')    # to import and save rasters
source('earthdata-login-info.R') # import personal login credentials for EarthData

# use a bounding box from tracking data ----
library('ctmm') # for movement data and models
tapir <- filter(readRDS('../tapirs/models/tapirs-final.rds'), name.short == 'ANNA')
tel <- tapir$data[[1]]
ud <- tapir$akde[[1]]

bbox <-
  SpatialPolygonsDataFrame.UD(ud, level.UD = 0.9995, level = 0) %>%
  st_as_sf() %>%
  st_transform(crs = '+proj=longlat') %>%
  st_bbox()

# download NDVI
MODIStsp(gui = FALSE, # do not use the browser GUI, only run in R
         out_folder = 'data/ndvi-rasters/tapir-anna', # '<folder>/VI_16Days_250m_v6/NDVI'
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
         parallel = 10) # use TRUE for automatic number of cores (max 8), or specify

# check rasters ----
ud_poly <- SpatialPolygonsDataFrame.UD(ud, level.UD = 0.9995, level = 0) %>%
  st_as_sf() %>%
  st_transform('+proj=longlat')

rasters <-
  list.files(path = 'data/ndvi-rasters/tapir-anna/VI_16Days_250m_v61/NDVI/',
             pattern = '.tif', full.names = TRUE) %>%
  rast()

if(FALSE) {
  rasters %>%
    mask(ud_poly) %>%
    crop(ud_poly) %>%
    terra::plot()
}

# save NDVI data as an rds file of a tibble
#' `rasterToPoints()` on a `stack()` of the rasters fails because there's too much data
rasters %>%
  as.data.frame(xy = TRUE) %>%
  pivot_longer(-c(x, y)) %>%
  transmute(long = x,
            lat = y,
            date = substr(name,
                          start = nchar('MOD13Q1_NDVI_x'),
                          stop = nchar(name)) %>%
              as.Date(format = '%Y_%j'),
            ndvi = value) %>%
  saveRDS('data/ndvi-rasters/tapir-anna/tapir-anna-data.rds')
