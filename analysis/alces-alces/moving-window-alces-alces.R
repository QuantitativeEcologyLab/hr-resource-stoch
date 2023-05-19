library('dplyr')   # for data wrangling
library('purrr')   # for functional programming
library('tidyr')   # for data wrangling
library('ggplot2') # for fancy plots
library('cowplot') # for fancy multi-panel plots
library('ctmm')    # for movement modeling
source('functions/window_hr.R') # function to calculate HRs and create figures
theme_set(theme_bw()) # change default theme

load('C:/Users/mezzinis/Dropbox/Uni/tels/Alces_alces.Rda')

tel <- filter(Alces_alces, location.long < -109) %>% # only NW data
  as.telemetry()

#' `Error in if (ANY) { : missing value where TRUE/FALSE needed` occasionally
#' occurs. Not sure why the error happens, but it can be ignored. Simply restart
#' the `for` loop
#' single-day HR has too little data
#' *117* (3):
#' *Error* in `mutate()`:
#' ! Problem while computing `ess = map_dbl(model, function(.m) summary(.m)$DOF["area"])`.
#' Caused by error in `summary(.m)$DOF`:
#' ! $ operator is invalid for atomic vectors`
for(ID in names(tel)) {
  message(paste0('Running moose ', ID, '...'))
  load(
    paste0('C:/Users/mezzinis/Dropbox/Empirical_Results/Alces_alces/Fits/Fits_',
           ID, '.rda'))
  window_hr(tel = tel[[ID]],
            window = 30 %#% 'day', # make sure sample size is high enough!
            dt = 5 %#% 'day',     # amount by which window moves over each time
            fig_path = 'figures/alces_alces',
            rds_path = 'models/alces_alces',
            cores = 1,
            movement_model = FIT)
}

# check output
if(FALSE) {
  d <- readRDS('models/alces_alces/10-window-3-days-dt-3-days.rds')
  head(select(d, 1:13))
  plot(hr_est_95 ~ date, d, type = 'b', pch = 20)
}
