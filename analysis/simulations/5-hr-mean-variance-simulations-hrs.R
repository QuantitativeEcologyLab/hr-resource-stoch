library('ctmm')  # for continuous-time movement modeling
library('dplyr') # for data wrangling (e.g., %>%)
library('purrr') # for functional programming (e.g., map(), map_dbl())
library('tidyr') # for data wrangling (e.g., nested tibbles)
source('analysis/mean-variance-trends-panel-data.R') # trends mu and s2
source('functions/get_hr.R') # for gaussian HR estimate from spatial variance

N <- nrow(d55)
files <- list.files('analysis/sockeye-scripts/temp-files/',
                    pattern = 'days-row*')

# check how many models have fit
length(files)

# ensure all models fit
length(files) == N

# check which models didn't fit
which(! (substr(files,
                nchar('days-row-x'),
                nchar(files) - nchar('.rds')) %in% 1:N))

days_models <-
  map_dfr(list.files('analysis/sockeye-scripts/temp-files/',
                     pattern = 'days-row-*'),
          \(x) readRDS(paste0('analysis/sockeye-scripts/temp-files/', x)))

days_hrs <-
  mutate(days_models,
         # mean exploration length
         t_expl = map_dbl(tel, \(x) as.numeric(max(x$timestamp))),
         # position variance
         pos_var = map_dbl(model, \(x) ctmm:::area.covm(x$sigma)),
         # calculate the HR as a function of the positional variance
         hr_50 = get_hr(pos_var, quantile = 0.50),
         hr_95 = get_hr(pos_var, quantile = 0.95)) %>%
  select(-c(tel, theta, model))

saveRDS(days_hrs, 'simulations/days-hrs.rds')
