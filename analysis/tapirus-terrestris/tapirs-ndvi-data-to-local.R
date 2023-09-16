# see https://rspatialdata.github.io/vegetation.html for example code and more info
#' install `rgeoboundaries` with `remotes::install_github('wmgeolab/rgeoboundaries')`
#' install `MODIStsp` with `remotes::install_github('ropensci/MODIStsp')`
#' `MODIStsp` version >= 2.0.7 from github fixes login issue due to 
library('ctmm') # for movement data and models
library('dplyr')    # for data wrangling
library('sf')       # for spatial features
library('MODIStsp') # for downloading NDVI rasters
library('purrr')    # for functional programming
library('tidyr')    # for data wrangling
library('terra')    # to import and save rasters
source('earthdata-login-info.R') # import personal login credentials for EarthData

# use a bounding box from tracking data ----
tapirs <- filter(readRDS('../tapirs/models/tapirs-final.rds'),
                 region == 'cerrado')

# extract all telemetries and convert them to a single, large data frame
tels <- map_dfr(tapirs$data, data.frame)

# import 95% home range estimates with no CIs (optional)
uds <-
  map_dfr(tapirs$akde, \(x) {
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
plot(uds, col = 'forestgreen', add = TRUE)

bbox <- st_bbox(uds_merged)

plot(st_as_sf(st_as_sfc(bbox)))
plot(uds_merged, col = 'grey', add = TRUE)
plot(uds, col = 'forestgreen', add = TRUE)

# download NDVI
MODIStsp(gui = FALSE, # do not use the browser GUI, only run in R
         reprocess =TRUE,
         out_folder = 'data/ndvi-rasters/tapirs', # '<folder>/VI_16Days_250m_v6/NDVI'
         selprod = 'Vegetation Indexes_16Days_250m (M*D13Q1)', # can't specify Terra here
         prod_version = '061', # 2022 raster version
         bandsel = 'NDVI', # Normalized Difference Vegetation Index layer only
         sensor = 'Terra', # only terrestrial values, ignore main bodies of water
         user = .USERNAME, # your Earthdata username (for urs.earthdata.nasa.gov/home)
         password = .PASSWORD, # your Earthdata password
         start_date = format(min(tels$timestamp) - 16 %#% 'days', # 1 raster before min
                             '%Y.%m.%d'),
         end_date = format(max(tels$timestamp) + 16 %#% 'days', # 1 raster after end
                           '%Y.%m.%d'),
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
rasters <-
  list.files(path = 'data/ndvi-rasters/tapirs/VI_16Days_250m_v61/NDVI/',
             pattern = '.tif', full.names = TRUE) %>%
  rast()

if(FALSE) {
  terra::plot(rasters)
}

# save NDVI data as an rds file of a tibble
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
  saveRDS('data/ndvi-rasters/tapirs/tapirs-data.rds')
