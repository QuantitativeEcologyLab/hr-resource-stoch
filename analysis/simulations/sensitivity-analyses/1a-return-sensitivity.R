library('ctmm')    # for generating movement models
library('dplyr')   # for data wrangling (%>%, mutate(), slice(), filter(), etc.)
library('purrr')   # for functional programming (map_*(), etc.)
library('tidyr')   # for data wrangling (pivot_*(), etc.)
library('ggplot2') # for fancy plots
source('analysis/simulations/movement-model.R') # for consistency between scripts
source('functions/add_returns.R') # to add return trips to (0, 0)
source('functions/get_hr.R') # to extract gaussian HR estimate from spatial variance
theme_set(theme_bw() + # change default ggplot theme
            theme(legend.position = 'none', strip.background = element_blank(),
                  strip.text = element_blank()))

DELTA_T <- 30 # sample every 30 s
SAMPLES <- seq(0, 60 * 60 * 24, by = DELTA_T) # 24 h sampled every DELTA_T

newd <- tibble(hour = 1:24, # 24 hours
               t = as.POSIXct(60 * 60 * hour, origin = '1970-01-01'),
               longitude = 0, # always at (0, 0) at the beginning of each hour
               latitude = 0,
               animal = 'test animal') %>%
  as.telemetry(projection = PROJECTION)
head(newd)
sims_loop <- simulate(model, data = newd, dt = DELTA_T, seed = 1, complete = TRUE)
sims_out <- sims_loop[sims_loop$t %% (60 * 60) < (60 * 30), ] # only first half hours
sims_out <- sims_out[sims_out$t < (60 * 60 * 24), ] # remove last data point

# check if `sims_loop` is continuous while `sims_out` is intermittent
plot(rep(1, nrow(sims_loop)) ~ sims_loop$t, pch = '.', ylim = c(-1, 2))
points(rep(0, nrow(sims_out)) ~ sims_out$t, pch = '.')

fit_ctmm <- function(.data) {
  tictoc::tic() # start the "stopwatch"
  theta <- ctmm.guess(.data, interactive = FALSE)
  model <- ctmm.fit(data = .data, CTMM = theta)
  tictoc::toc() # print the elapsed time
  return(model)
}
                              # DELTA_t     : 5 min, 1 min,  30 s
m_loop <- fit_ctmm(sims_loop) # Fitting time:  16 s,  65 s, 168 s
m_out <- fit_ctmm(sims_out)   # Fitting time:   8 s,  32 s,  86 s

summary(m_loop)
summary(m_out)
