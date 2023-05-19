library('dplyr')     # for data wrangling
library('purrr')     # for functional programming
library('tidyr')     # for data wrangling
library('ctmm')      # for movement modeling
library('mgcv')      # for empirical Bayesian modeling
library('lubridate') # for smoother date wrangling
library('ggplot2')   # for fancy figures
library('cowplot')   # for fancy multi-panel plots

# for common theme and color palette
source('analysis/figures/barga-italy-2023/default-figure-styling-barga-2023.R')

theme_set(theme_get() + theme(panel.grid = element_blank()))

# function to make predictions from the model
ndvi_preds <- function(.data) {
  .data <- .data %>%
    data.frame() %>% # convert telemetry data to a data.frame
    mutate(year = year(timestamp),
           doy = yday(timestamp))
  
  predict.gam(m_ndvi, newdata = .data, type = 'response', se.fit = FALSE) %>%
    data.frame() %>%
    tibble() %>%
    transmute(mu = X1,
              sigma2 = (1/X2)^2) %>%
    return()
}

# import moose moving window data
load('C:/Users/mezzinis/Dropbox/Uni/tels/Alces_alces.Rda')
Alces_alces <- filter(Alces_alces, location.long < -109) # only NW data

# import NDVI data
d <- readRDS('data/ndvi-rasters/alces-alces-nw/alces-alces-ndvi.rds') %>%
  mutate(dec_date = decimal_date(date),
         year = year(date),
         doy = yday(date))

# extract the animal's telemetry
tel <- filter(Alces_alces, individual.local.identifier == '117') %>%
  data.frame() %>%
  rename(long = location.long, lat = location.lat) %>%
  mutate(dec_date = decimal_date(timestamp),
         year = year(timestamp),
         doy = yday(timestamp))

# subset the data to the area of the animal
d_0 <- filter(d,
              long > min(tel$long) - 0.05,
              long < max(tel$long) + 0.05,
              lat > min(tel$lat) - 0.05,
              lat < max(tel$lat) + 0.05,
              dec_date > min(tel$dec_date) - 0.1,
              dec_date < max(tel$dec_date) + 0.1)

if(FALSE) {
  d_0 %>%
    filter(date == min(date)) %>%
    ggplot() +
    geom_raster(aes(long, lat, fill = ndvi)) +
    geom_point(aes(long, lat), data = tel)
}

if(! file.exists('models/alces_alces/117-ndvi-gaulss.rds')) {
  # fit a location-scale model to the subset of the data
  DATES <- floor(n_distinct(d_0$dec_date))
  
  cat('Fitting GAM...\n')
  m_ndvi <-
    gam(list(
      # mean predictor
      ndvi ~
        s(long, lat, bs = 'ds', k = 15) + # mean over space
        s(dec_date, bs = 'cr', k = DATES) + # change over time
        ti(long, lat, dec_date, bs = c('ds', 'cr'), d = c(2, 1), k = c(30, 10)),
      # precision (1/standard deviation) predictor
      ~
        s(long, lat, bs = 'ds', k = 15) + # mean over space
        s(dec_date, bs = 'cr', k = round(DATES/4)) + # change over time
        ti(long, lat, dec_date, bs = c('ds', 'cr'), d = c(2, 1), k = c(10, 10))),
      family = gaulss(),
      data = d_0,
      method = 'REML',
      knots = list(doy = c(0.5, 366.5)),
      control = gam.control(nthreads = 8, trace = TRUE))
  
  saveRDS(m_ndvi, 'models/alces_alces/117-ndvi-gaulss.rds')
  
  if(FALSE) {
    summary(m_ndvi)
    plot(m_ndvi, pages = 1, scheme = 3, scale = 0) # plot smooths
  }
} else {
  m_ndvi <- readRDS('models/alces_alces/117-ndvi-gaulss.rds')
}

tel <- bind_cols(tel, ndvi_preds(tel))

moose <-
  readRDS('models/alces_alces/117-window-30-days-dt-5-days.rds') %>%
  filter(! is.na(hr_est_95)) %>%
  mutate(ess = map_dbl(model, \(.m) summary(.m)$DOF['area']),
         weight = ess / mean(ess), # weights based on effective sample size
         preds = map(dataset, \(.d) filter(tel, timestamp %in% .d$timestamp)),
         mu = map_dbl(preds, \(.d) mean(.d$mu)),
         sigma2 = map_dbl(preds, \(.d) mean(.d$sigma2))) %>%
  select(t_center, mu, sigma2, hr_lwr_95, hr_est_95, hr_upr_95, ess, weight) %>%
  pivot_longer(c(hr_lwr_95, hr_est_95, hr_upr_95), names_to = c('.value', 'quantile'),
               names_pattern = '(.+)_(.+)') %>%
  mutate(t_center = as.POSIXct(t_center, origin = origin),
         quantile = paste0(quantile, '%'),
         quantile = factor(quantile))

# create the figure
date_labs <- range(tel$timestamp) %>% as.Date()
# YLIMS <- c(0, 50)

# not using a conditional model because E(R) and V(R) are highly correlated
ggplot(moose) +
  geom_point(aes(mu, sigma2), alpha = 0.5) +
  labs(x = '\U1D707(t)', y = '\U1D70E\U00B2(t)')

grobs <-
  lapply(
    list(
      # mean, variance and HR
      ggplot(moose, aes(t_center, mu)) +
        geom_line(color = pal[1], linewidth = 2) +
        labs(x = NULL, y = 'Resource abundance, \U1D53C(\U1D445)'),
      ggplot(moose, aes(t_center, sigma2)) +
        geom_line(color = pal[2], linewidth = 2) +
        labs(x = NULL, y = 'Resource stochasticity, \U1D54D(\U1D445)'),
      ggplot(moose) +
        geom_line(aes(t_center, hr_est), color = pal[3], linewidth = 2) +
        labs(x = NULL, y = expression('30-day'~home~range~size~(km^2)~'   '))),
    as_grob) # convert to grid graphical objects (grobs)

# align left margins of all plots
aligned_widths <- align_margin(lapply(grobs, function(x) {x$widths}), 'first')

# Setting the dimensions of plots to the aligned dimensions
for (i in seq_along(grobs)) {
  grobs[[i]]$widths <- aligned_widths[[i]]
}

# Draw aligned plots
p_1 <- plot_grid(plotlist = grobs, ncol = 1, labels = c('a.', 'b.', 'c.'),
                 label_size = 22)

p_2 <-
  plot_grid(
    # mean
    ggplot(moose) +
      # coord_cartesian(ylim = YLIMS) +
      geom_point(aes(mu, hr_est), alpha = 0.5, color = pal[3]) +
      geom_smooth(aes(mu, hr_est), color = pal[1], se = FALSE, method = 'gam',
                  formula = y ~ x, linewidth = 2,
                  method.args = list(family = 'Gamma')) +
      scale_x_continuous('Resource abundance') +
      scale_y_continuous(expression('30-day'~home~range~size~(km^2)),
                         expand = c(0, 0)) +
      scale_size_manual(values = c(1, 0.5)),
    
    # variance
    ggplot(moose) +
      # coord_cartesian(ylim = YLIMS) +
      geom_point(aes(sigma2, hr_est), alpha = 0.5, color = pal[3]) +
      geom_smooth(aes(sigma2, hr_est), color = pal[2], se = FALSE,
                  method = 'gam', formula = y ~ x, linewidth = 2,
                  method.args = list(family = 'Gamma')) +
      scale_x_continuous('Resource stochasticity') +
      scale_y_continuous(expression('30-day'~home~range~size~(km^2)),
                         expand = c(0, 0)),
    ncol = 1)
p_2

p_2 <-
  moose %>%
  pivot_longer(c(mu, sigma2), names_to = 'param', values_to = 'value') %>%
  mutate(param = if_else(param == 'mu', 'Resource abundance',
                         'Resource stochasticity')) %>%
  ggplot() +
  # coord_cartesian(ylim = YLIMS) +
  facet_wrap(~ param, scales = 'free_x', strip.position = 'bottom', ncol = 1) +
  geom_point(aes(value, hr_est), alpha = 0.5, color = pal[3]) +
  geom_smooth(aes(value, hr_est, color = param), se = FALSE,
              method = 'gam', formula = y ~ x, linewidth = 2,
              method.args = list(family = 'Gamma'), show.legend = FALSE) +
  xlab(NULL) +
  scale_y_continuous(expression('30-day'~home~range~size~(km^2)),
                     expand = c(0, 0)) +
  scale_color_manual(values = pal) +
  theme(strip.placement = 'outside', strip.background = element_blank(),
        strip.text = element_text(size = 22))

p <- plot_grid(p_1, p_2, ncol = 2, labels = c('', 'd.'))

ggsave('figures/2023-GRS-GRC-barga-italy/alces-alces/moose-117-example.png',
       plot = p, height = 40, width = 40, units = 'cm', dpi = 600,
       bg = 'white')

m <- gam(hr_est ~ mu + sigma2 + mu:sigma2,
         family = Gamma('log'),
         data = moose,
         method = 'REML')

# make predictions from the model
preds <-
  tibble(
    mu = seq(min(moose$mu), max(moose$mu), length.out = 400),
    sigma2 = seq(min(moose$sigma2), max(moose$sigma2), length.out = 400)) %>%
  bind_cols(.,
            hr_mu = predict(m, newdata = ., type = 'response',
                            terms = c('(Intercept), mu')) %>%
              as.numeric(),
            hr_sigma2 = predict(m, newdata = ., type = 'response',
                                terms = c('(Intercept)', 'sigma2')) %>%
              as.numeric()) %>%
  pivot_longer(c(hr_mu, hr_sigma2), names_to = 'pred', names_pattern = 'hr_(.+)',
               values_to = 'hr') %>%
  pivot_longer(c(mu, sigma2), names_to = 'pred2', values_to = 'x') %>%
  filter(pred == pred2) %>%
  mutate(param = if_else(pred == 'mu', 'Resource abundance',
                         'Resource stochasticity'))

moose %>%
  pivot_longer(c(mu, sigma2), names_to = 'param', values_to = 'value') %>%
  mutate(param = if_else(param == 'mu', 'Resource abundance',
                         'Resource stochasticity')) %>%
  ggplot() +
  # coord_cartesian(ylim = c(0, 1e4)) +
  facet_wrap(~ param, scales = 'free_x', strip.position = 'bottom', ncol = 1) +
  geom_point(aes(value, hr_est), alpha = 0.5, color = pal[3]) +
  geom_line(aes(x, hr, color = param), preds, linewidth = 2,
            show.legend = FALSE) +
  xlab(NULL) +
  scale_y_continuous(expression('30-day'~home~range~size~(km^2)),
                     expand = c(0, 0)) +
  scale_color_manual(values = pal) +
  theme(strip.placement = 'outside', strip.background = element_blank(),
        strip.text = element_text(size = 22))


plot(m, all.terms = TRUE, pages = 1)
