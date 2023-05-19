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
m_ndvi <- readRDS('models/CE_31_ANNA-mgcv-ndvi-gaulss.rds')

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

# import tapir moving window data
tapir <- readRDS('../tapirs/models/tapirs-final.rds') %>%
  filter(name.short == 'ANNA')
tel <- tapir$data[[1]]

if(FALSE) {
  m_ouf <- tapir$model[[1]]
  track <-
    predict(m_ouf, data = tel, dt = 100, complete = TRUE) %>%
    data.frame() %>%
    select(timestamp, longitude, latitude) %>%
    rename(long = longitude, lat = latitude) %>%
    mutate(dec_date = decimal_date(timestamp))
  
  saveRDS(track, 'models/CE_31_ANNA-track-dt-100.rds')
} else {
  track <- readRDS('models/CE_31_ANNA-track-dt-100.rds')
}

tel <- data.frame(tel) %>%
  rename(long = longitude, lat = latitude) %>%
  mutate(dec_date = decimal_date(timestamp))
tel <- bind_cols(tel, ndvi_preds(tel))

tapir <-
  readRDS('models/CE_31_ANNA-window-7-days-dt-1-days.rds') %>%
  # mutate(ess = map_dbl(model, \(.m) summary(.m)$DOF['area']),
  #        weight = ess / mean(ess)) %>%
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
YLIMS <- c(0, 13)

# not using a conditional model because E(R) and V(R) are highly correlated
ggplot(tapir) +
  geom_point(aes(mu, sigma2), alpha = 0.1) +
  labs(x = '\U1D707(t)', y = '\U1D70E\U00B2(t)')

grobs <-
  lapply(
    list(
      # mean, variance and HR
      ggplot(tapir, aes(t_center, mu)) +
        geom_line(color = pal[1], linewidth = 2) +
        labs(x = NULL, y = 'Resource abundance, \U1D53C(\U1D445)'),
      ggplot(tapir, aes(t_center, sigma2)) +
        geom_line(color = pal[2], linewidth = 2) +
        labs(x = NULL, y = 'Resource stochasticity, \U1D54D(\U1D445)'),
      ggplot(tapir) +
        geom_line(aes(t_center, hr_est), color = pal[3], linewidth = 2) +
        labs(x = NULL, y = expression('7-day'~home~range~size~(km^2)~'   '))),
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

# regression plots ----
# cannot use a GLM or GAM because mu and sigma2 are correlated
plot(sigma2 ~ mu, tapir)

p_2 <-
  tapir %>%
  pivot_longer(c(mu, sigma2), names_to = 'param', values_to = 'value') %>%
  mutate(param = if_else(param == 'mu', 'Resource abundance',
                         'Resource stochasticity')) %>%
  ggplot() +
  coord_cartesian(ylim = YLIMS) +
  facet_wrap(~ param, scales = 'free_x', strip.position = 'bottom', ncol = 1) +
  geom_point(aes(value, hr_est), alpha = 0.5, color = pal[3]) +
  geom_smooth(aes(value, hr_est, color = param), se = FALSE,
              method = 'gam', formula = y ~ x, linewidth = 2,
              method.args = list(family = 'Gamma'), show.legend = FALSE) +
  xlab(NULL) +
  scale_y_continuous(expression('7-day'~home~range~size~(km^2)),
                     expand = c(0, 0)) +
  scale_color_manual(values = pal) +
  theme(strip.placement = 'outside', strip.background = element_blank(),
        strip.text = element_text(size = 22)); p_2

p <- plot_grid(p_1, p_2, ncol = 2, labels = c('', 'd.'))

ggsave('figures/2023-GRS-GRC-barga-italy/tapir-example.png', plot = p,
       scale = 4/3, height = 30, width = 40, units = 'cm', dpi = 600,
       bg = 'white')
