library('ctmm')   # for movement modeling
library('dplyr')  # for data wrangling (%>%, mutate())
library('terra')  # to work with raster data

# label each movement as whether a new patch was visited
label_visits <- function(.tel, .habitat) {
  .tel %>%
    data.frame() %>% # convert to a data.frame for easy plotting
    dplyr::mutate(cell_id =
                    terra::cellFromXY(object = .habitat,
                                      xy = SpatialPoints.telemetry(.tel) %>%
                                        as.data.frame() %>%
                                        select(x, y)) %>%
                    suppressWarnings(),
                  # check if animal moved to a new cell, first cell is a zero
                  new_cell = c(1, diff(cell_id)),
                  new_cell = new_cell != 0) # convert to TRUE/FALSE
}
