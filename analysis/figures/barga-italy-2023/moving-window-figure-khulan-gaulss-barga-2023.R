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

theme_set(theme_get() +
            theme(panel.grid = element_blank()))

# import location-scale model for predictions of mean and variance in NDVI
m_ndvi <- readRDS('models/equus-hemionus-hemionus-4-mgcv-ndvi-gaulss-100-reduction.rds')

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

# import khulan moving window data
load('data/equus-hemionus-hemionus/Equus_hemionus_hemionus.Rda')
tel <- filter(Equus_hemionus_hemionus, Collar_ID == '1') %>%
  as.telemetry()
load('data/equus-hemionus-hemionus/Fits/Fits_4.rda')
m_ouf <- FIT

tel@info['identity'] <- 'equus-hemonius-hemonius-4'

rm(Equus_hemionus_hemionus, FIT)

tel <- data.frame(tel) %>%
  rename(long = longitude, lat = latitude) %>%
  mutate(dec_date = decimal_date(timestamp))
tel <- bind_cols(tel, ndvi_preds(tel))

khulan <-
  readRDS('models/equus-hemonius-hemonius-4-window-30-days-dt-1-days.rds') %>%
  filter(! is.na(hr_est_95)) %>%
  mutate(preds = map(dataset, \(.d) filter(tel, timestamp %in% .d$timestamp)),
         mu = map_dbl(preds, \(.d) mean(.d$mu)),
         sigma2 = map_dbl(preds, \(.d) mean(.d$sigma2))) %>%
  select(t_center, mu, sigma2, hr_lwr_95, hr_est_95, hr_upr_95) %>%
  pivot_longer(c(hr_lwr_95, hr_est_95, hr_upr_95), names_to = c('.value', 'quantile'),
               names_pattern = '(.+)_(.+)') %>%
  mutate(t_center = as.POSIXct(t_center, origin = origin),
         quantile = paste0(quantile, '%'),
         quantile = factor(quantile))

# create the figure
date_labs <- range(tel$timestamp) %>% as.Date()
# YLIMS <- c(0, 13)

# not using a conditional model because E(R) and V(R) are highly correlated
ggplot(khulan) +
  geom_point(aes(mu, sigma2)) +
  geom_smooth(aes(mu, sigma2), method = 'gam') +
  labs(x = '\U1D707(t)', y = '\U1D70E\U00B2(t)')

grobs <-
  lapply(
    list(
      # mean, variance and HR
      ggplot(khulan, aes(t_center, mu)) +
        geom_line(color = pal[1], linewidth = 2) +
        labs(x = NULL, y = 'Resource abundance, \U1D53C(\U1D445)'),
      ggplot(khulan, aes(t_center, sigma2)) +
        geom_line(color = pal[2], linewidth = 2) +
        labs(x = NULL, y = 'Resource stochasticity, \U1D54D(\U1D445)'),
      ggplot(khulan) +
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
    ggplot(khulan) +
      # coord_cartesian(ylim = YLIMS) +
      geom_point(aes(mu, hr_est), alpha = 0.3, color = pal[3]) +
      geom_smooth(aes(mu, hr_est), color = pal[1], se = FALSE, method = 'gam',
                  formula = y ~ x, method.args = list(family = "Gamma"),
                  linewidth = 2) +
      scale_x_continuous('Resource abundance') +
      scale_y_continuous(expression('30-day'~home~range~size~(km^2)),
                         expand = c(0, 0)) +
      scale_size_manual(values = c(1, 0.5)),
    
    # variance
    ggplot(khulan) +
      # coord_cartesian(ylim = YLIMS) +
      geom_point(aes(sigma2, hr_est), alpha = 0.3, color = pal[3]) +
      geom_smooth(aes(sigma2, hr_est), color = pal[2], se = FALSE,
                  method = 'gam', formula = y ~ x, linewidth = 2,
                  method.args = list(family = "Gamma")) +
      scale_x_continuous('Resource stochasticity') +
      scale_y_continuous(expression('30-day'~home~range~size~(km^2)),
                         expand = c(0, 0)),
    ncol = 1)
p_2

p_2 <-
  khulan %>%
  pivot_longer(c(mu, sigma2), names_to = 'param', values_to = 'value') %>%
  mutate(param = if_else(param == 'mu', 'Resource abundance',
                         'Resource stochasticity')) %>%
  ggplot() +
  # coord_cartesian(ylim = YLIMS) +
  facet_wrap(~ param, scales = 'free_x', strip.position = 'bottom', ncol = 1) +
  geom_point(aes(value, hr_est), alpha = 0.5, color = pal[3]) +
  geom_smooth(aes(value, hr_est, color = param), se = FALSE,
              method = 'gam', formula = y ~ x, linewidth = 2,
              method.args = list(family = "Gamma"), show.legend = FALSE) +
  xlab(NULL) +
  scale_y_continuous(expression('30-day'~home~range~size~(km^2)),
                     expand = c(0, 0)) +
  scale_color_manual(values = pal) +
  theme(strip.placement = 'outside', strip.background = element_blank(),
        strip.text = element_text(size = 22)); p_2

p <- plot_grid(p_1, p_2, ncol = 2, labels = c('', 'd.'))

ggsave('figures/2023-GRS-GRC-barga-italy/khulan-example.png', plot = p,
       height = 40, width = 40, units = 'cm', dpi = 600, bg = 'white')

