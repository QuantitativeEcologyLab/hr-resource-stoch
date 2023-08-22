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
load('C:/Users/mezzinis/Dropbox/Uni/tels/Alces_alces.Rda')

# filter to north-western individuals to avoid having a spatial model with a
# large spatial gap
Alces_alces %>%
  group_by(individual.local.identifier) %>%
  summarize(location.lat = mean(location.lat),
            location.long = mean(location.long)) %>%
  plot(location.lat ~ location.long, ., col = 'red', pch = 19)

nw <- Alces_alces %>%
  group_by(individual.local.identifier) %>%
  summarize(nw = all(location.long < -109)) %>%
  pull(nw) %>%
  which()

Alces_alces <- filter(Alces_alces, location.long < -109)

fits <- lapply(files, \(x) {
  load(x)
  return(FIT)
})

hrs <- map_dbl(fits, \(x) summary(x, units = FALSE)$CI['area (square meters)', 'est'])
tp <- map_dbl(fits, \(x) summary(x, units = FALSE)$CI['τ[position] (seconds)', 'est'])

# more migratory individuals in the NW clusters
layout(matrix(1:4, ncol = 2))
hist(hrs / (1 %#% 'km^2'), xlim = c(0, 7000))
hist(hrs[nw] / (1 %#% 'km^2'), xlim = c(0, 7000))
hist(tp / (1 %#% 'months'), breaks = 20, xlim = c(0, 15))
hist(tp[nw] / (1 %#% 'months'), breaks = 20, xlim = c(0, 15))
layout(1)

Alces_alces %>%
  group_by(individual.local.identifier) %>%
  summarise(nrows = n(),
            days = round((max(timestamp, na.rm = TRUE) - min(timestamp, na.rm = TRUE))),
            obs_per_day = nrows / as.numeric(days)) %>%
  arrange(desc(days))

Alces_alces_sp <- as.telemetry(Alces_alces) %>%
  SpatialPoints.telemetry() %>%
  st_as_sf() %>%
  st_transform('+proj=longlat')

bbox <- Alces_alces_sp %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_buffer(dist = 5e3) %>% # buffer 5 km
  st_bbox() %>%
  st_as_sfc()

plot(bbox)
plot(Alces_alces_sp, add = TRUE, col = 'red', pch = 19)

# download NDVI
MODIStsp(gui = FALSE, # do not use the browser GUI, only run in R
         out_folder = 'data/ndvi-rasters/alces-alces-nw/', # '<folder>/VI_16Days_1Km_v6/NDVI'
         selprod = 'Vegetation_Indexes_16Days_1Km (M*D13A2)', # can't specify Terra here
         prod_version = '061', # 2022 raster version
         bandsel = 'NDVI', # Normalized Difference Vegetation Index layer only
         sensor = 'Terra', # only terrestrial values, ignore main bodies of water
         user = .USERNAME, # your Earthdata username (for urs.earthdata.nasa.gov/home)
         password = .PASSWORD, # your Earthdata password
         start_date = format(min(Alces_alces$timestamp) - 16, '%Y.%m.%d'), # 1 raster before min
         end_date = format(max(Alces_alces$timestamp) + 16, '%Y.%m.%d'), # 1 raster after end
         spatmeth = 'bbox', # use a bounding box for the extent
         bbox = st_bbox(bbox), # spatial file for raster extent
         out_projsel = 'User Defined', # use specified projection instead of default
         output_proj = '+proj=longlat', # download unprojected raster
         resampling = 'bilinear', # method for resampling raster if changing projection
         delete_hdf = TRUE, # delete HDF files after download is complete
         scale_val = TRUE, # convert from integers to floats within [-1, 1]
         out_format = 'GTiff', # output format: 'ENVI' (.hdr) or 'GTiff' (.tif)
         n_retries = 10, # number of times to try again if download fails before aborting
         verbose = TRUE, # print processing messages
         parallel = TRUE) # use TRUE for automatic number of cores (max 8), or specify

# check rasters ----
rasters <-
  list.files(path = 'data/ndvi-rasters/alces-alces-nw/VI_16Days_1Km_v61/NDVI/',
             pattern = '.tif', full.names = TRUE) %>%
  rast()

names(rasters) <-
  list.files(path = 'data/ndvi-rasters/alces-alces-nw/VI_16Days_1Km_v61/NDVI/',
             pattern = '.tif', full.names = FALSE) %>%
  substr(start = nchar('MOD13Q1_NDVI_X'), stop = nchar('MOD13Q1_NDVI_XXXX_XXX')) %>%
  as.Date(format = '%Y_%j')

if(FALSE) {
  terra::plot(rasters[[1]])
  plot(Alces_alces_sp, add = TRUE)
  terra::plot(rasters)
  
  for(i in seq(1, length(names(rasters)), by = 16)) {
    terra::plot(rasters[[i:(i + 15)]])
    readline(prompt = 'Press <Enter> to view next set of plots')
  }
}

# save NDVI data as an rds file of a tibble
rasters %>%
  as.data.frame(xy = TRUE) %>%
  pivot_longer(-c(x, y)) %>%
  transmute(long = x,
            lat = y,
            date = as.Date(name),
            ndvi = value) #%>%
  saveRDS('data/ndvi-rasters/alces-alces-nw/alces-alces-ndvi.rds')
