library('dplyr')     # for data wrangling
library('purrr')     # for functional programming
library('tidyr')     # for data wrangling
library('ctmm')      # for movement modeling
library('mgcv')      # for GAMs
library('scam')      # for shape-constrained GAMs
library('lubridate') # for smoother date wrangling
library('ggplot2')   # for fancy figures
library('cowplot')   # for fancy multi-panel plots

# for common theme and color palette
source('analysis/figures/barga-italy-2023/default-figure-styling-barga-2023.R')

theme_set(theme_get() +
            theme(panel.grid = element_blank()))

# import location-scale model for predictions of mean and variance in NDVI
m_ndvi <- readRDS('models/aepyceros-melampus/1757-ndvi-gaulss.rds')

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

# import impala moving window data
load('C:/Users/mezzinis/Dropbox/Empirical_Results/Aepyceros_melampus/Fits/Fits_1757.rda')
load('C:/Users/mezzinis/Dropbox/Uni/tels/Aepyceros_melampus.Rda')
tel <- Aepyceros_melampus %>%
  filter(individual.local.identifier == '1757') %>%
  as.telemetry()
plot(tel)
diff(range(tel$timestamp))

tel <- data.frame(tel) %>%
  rename(long = longitude, lat = latitude) %>%
  mutate(dec_date = decimal_date(timestamp))
tel <- bind_cols(tel, ndvi_preds(tel))

impala <-
  readRDS('models/aepyceros-melampus-1757-window-7-days-dt-1-days.rds') %>%
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

hist(impala$hr_est)
# drop outlier HRs
impala <- filter(impala, hr_est <= 7.5)
hist(impala$hr_est)

# create the figure
date_labs <- range(tel$timestamp) %>% as.Date()
# YLIMS <- c(0, 5)

# not using a conditional model because E(R) and V(R) are highly correlated
ggplot(impala) +
  geom_point(aes(mu, sigma2)) +
  labs(x = '\U1D707(t)', y = '\U1D70E\U00B2(t)')

grobs <-
  lapply(
    list(
      # mean, variance and HR
      ggplot(impala, aes(t_center, mu)) +
        geom_line(color = pal[1], linewidth = 2) +
        labs(x = NULL, y = 'Resource abundance, \U1D53C(\U1D445)'),
      ggplot(impala, aes(t_center, sigma2)) +
        geom_line(color = pal[2], linewidth = 2) +
        labs(x = NULL, y = 'Resource stochasticity, \U1D54D(\U1D445)'),
      ggplot(impala) +
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

# fit regression ----
m <- gam(hr_est ~ mu + sigma2,
         family = Gamma('log'),
         data = impala)
plot(m, pages = 1, all.terms = TRUE)
layout(matrix(1:4, ncol = 2))
gam.check(m)
layout(1)

preds <-
  tibble(
    mu = seq(min(impala$mu), max(impala$mu), length.out = 400),
    sigma2 = seq(min(impala$sigma2), max(impala$sigma2), length.out = 400)) %>%
  bind_cols(.,
            hr_mu = predict(m, newdata = ., type = 'response',
                            terms = c('(Intercept)', 'mu')) %>%
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

p_2 <-
  impala %>%
  pivot_longer(c(mu, sigma2), names_to = 'param', values_to = 'value') %>%
  mutate(param = if_else(param == 'mu', 'Resource abundance',
                         'Resource stochasticity')) %>%
  ggplot() +
  # coord_cartesian(ylim = YLIMS) +
  facet_wrap(~ param, scales = 'free_x', strip.position = 'bottom', ncol = 1) +
  geom_point(aes(value, hr_est), alpha = 0.5, color = pal[3]) +
  geom_line(aes(x, hr, color = param), preds, linewidth = 2,
            show.legend = FALSE) +
  xlab(NULL) +
  scale_y_continuous(expression('7-day'~home~range~size~(km^2))) +
  scale_color_manual(values = pal) +
  theme(strip.placement = 'outside', strip.background = element_blank(),
        strip.text = element_text(size = 22))

p <- plot_grid(p_1, p_2, ncol = 2, labels = c('', 'd.'))

ggsave('figures/2023-GRS-GRC-barga-italy/impala-example.png', plot = p,
       scale = 4/3, height = 30, width = 40, units = 'cm', dpi = 600,
       bg = 'white')
