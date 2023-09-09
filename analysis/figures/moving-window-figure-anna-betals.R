library('dplyr')     # for data wrangling
library('purrr')     # for functional programming
library('tidyr')     # for data wrangling
library('ctmm')      # for movement modeling
library('mgcv')      # for empirical Bayesian modeling
library('scam')      # for shape-constrained splines
library('lubridate') # for smoother date wrangling
library('ggplot2')   # for fancy figures
library('cowplot')   # for fancy multi-panel plots

source('analysis/figures/default-figure-styling.R') # for theme & palettes

# axis lables
e_r <- 'Resource abundance, \U1D707(t)'
v_r <- 'Resource stochasticity, \U1D70E\U00B2(t)'
hr_lab <- expression('7-day'~home~range~size~(km^2)~'   ')

# import location-scale model for predictions of mean and variance in NDVI
m_ndvi <- readRDS('models/tapirs/CE_31_ANNA-mgcv-ndvi-betals.rds')

# function to make predictions from the model for each telemetry
ndvi_preds <- function(.data) {
  .data <- data.frame(.data) # convert telemetry data to a data.frame
  
  predict(m_ndvi, newdata = .data, type = 'response', se.fit = FALSE) %>%
    data.frame() %>%
    tibble() %>%
    transmute(mu = X1,
              sigma2 = X1 * (1 - X1) * X2) %>% # didn't scale NDVI to [0, 1]
    return()
}

# import tapir data (from https://doi.org/10.1186/s40462-022-00313-w)
tapir <- readRDS('../tapirs/models/tapirs-final.rds') %>%
  filter(name.short == 'ANNA')

tel <- data.frame(tapir$data[[1]]) %>%
  rename(long = longitude, lat = latitude) %>%
  mutate(dec_date = decimal_date(timestamp)) %>%
  bind_cols(., ndvi_preds(.)) %>% 
  as_tibble()

tapir <-
  readRDS('models/tapirs/CE_31_ANNA-window-7-days-dt-1-days.rds') %>%
  mutate(est = map(dataset, \(.d) filter(tel, timestamp %in% .d$timestamp)),
         mu = map_dbl(est, \(.d) mean(.d$mu)),
         sigma2 = map_dbl(est, \(.d) mean(.d$sigma2))) %>%
  select(t_center, mu, sigma2, hr_lwr_95, hr_est_95, hr_upr_95) %>%
  pivot_longer(c(hr_lwr_95, hr_est_95, hr_upr_95), names_to = c('.value', 'quantile'),
               names_pattern = '(.+)_(.+)') %>%
  mutate(t_center = as.POSIXct(t_center),
         quantile = paste0(quantile, '%'))

# create the figure
date_labs <- range(tel$timestamp) %>% as.Date()
YLIMS <- c(0, 13)

# E(R) and V(R) are highly correlated
ggplot(tapir) +
  geom_point(aes(mu, sigma2), alpha = 0.1) +
  labs(x = '\U1D707(t)', y = '\U1D70E\U00B2(t)')

l_grobs <-
  lapply(
    list(
      # mean, variance and HR
      ggplot(tapir, aes(t_center, mu)) +
        geom_line(color = pal[1], linewidth = 2) +
        labs(x = NULL, y = e_r),
      ggplot(tapir, aes(t_center, sigma2)) +
        geom_line(color = pal[2], linewidth = 2) +
        labs(x = NULL, y = v_r),
      ggplot(tapir) +
        geom_line(aes(t_center, hr_est), color = pal[3], linewidth = 2) +
        labs(x = NULL, y = hr_lab)),
    as_grob) # convert to grid graphical objects (grobs)

# align left margins of all plots
aligned_widths <- align_margin(map(l_grobs, function(x) {x$widths}), 'first')

# Setting the dimensions of plots to the aligned dimensions
for (i in seq_along(l_grobs)) {
  l_grobs[[i]]$widths <- aligned_widths[[i]]
}

# Draw aligned plots
p_left <- plot_grid(plotlist = l_grobs, ncol = 1, labels = c('a.', 'b.', 'c.'),
                    label_size = 22, label_x = -0.01, label_y = 1.025)

# regression plots ----
# using shape-constrained splines to avoid artifacts
m <- scam(hr_est ~ s(mu, bs = 'mpd', k = 4) + s(sigma2, bs = 'mpi', k = 4),
          family = Gamma('log'),
          data = tapir)
plot(m, pages = 1, scheme = 2, trans = exp, ylim = c(log(0), log(4)))

# removing data from d and splitting into two panels
preds <- tibble(mu = gratia:::seq_min_max(tapir$mu, n = 250),
                sigma2 = gratia:::seq_min_max(tapir$sigma2, n = 250)) %>%
  bind_cols(predict(m, newdata = ., terms = 's(mu)',
                    type = 'link', se.fit = TRUE) %>%
              as.data.frame() %>%
              transmute(hr_mu_est = exp(fit),
                        hr_mu_lwr = exp(fit - 1.96 * se.fit),
                        hr_mu_upr = exp(fit + 1.96 * se.fit)),
            predict(m, newdata = ., terms = 's(sigma2)',
                    type = 'link', se.fit = TRUE) %>%
              as.data.frame() %>%
              transmute(hr_sigma2_est = exp(fit),
                        hr_sigma2_lwr = exp(fit - 1.96 * se.fit),
                        hr_sigma2_upr = exp(fit + 1.96 * se.fit)))

p_d <- ggplot() +
  coord_cartesian(ylim = c(0, 12.5)) +
  geom_point(aes(mu, hr_est), tapir, color = pal[3], alpha = 0.5) +
  geom_ribbon(aes(mu, ymin = hr_mu_lwr, ymax = hr_mu_upr), preds,
              fill = pal[1], alpha = 0.3) +
  geom_line(aes(mu, hr_mu_est), preds, color = pal[1], linewidth = 2) +
  xlab(e_r) +
  scale_y_continuous(hr_lab)

p_e <- ggplot() +
  coord_cartesian(ylim = c(0, 12.5)) +
  geom_point(aes(sigma2, hr_est), tapir, color = pal[3], alpha = 0.5) +
  geom_ribbon(aes(sigma2, ymin = hr_sigma2_lwr, ymax = hr_sigma2_upr), preds,
              fill = pal[2], alpha = 0.3) +
  geom_line(aes(sigma2, hr_sigma2_est), preds, color = pal[2], linewidth = 2) +
  xlab(v_r) +
  scale_y_continuous(hr_lab)

p_f <- ggplot(tapir) +
  geom_point(aes(mu, sigma2), alpha = 0.5) +
  labs(x = e_r, y = v_r)

r_grobs <- map(list(p_d, p_e, p_f), as_grob)

# align right margins of all plots
aligned_widths <- align_margin(map(r_grobs, function(x) {x$widths}), 'first')

# Setting the dimensions of plots to the aligned dimensions
for (i in seq_along(r_grobs)) {
  r_grobs[[i]]$widths <- aligned_widths[[i]]
}

p_right <- plot_grid(plotlist = r_grobs, ncol = 1, labels = c('d.', 'e.', 'f.'),
                     label_size = 22, label_x = -0.01, label_y = 1.025)

# add space to avoid cutting off x axis for (e.)
p <- plot_grid(p_left, p_right, NULL, nrow = 1, rel_widths = c(1, 1, 0.01))

ggsave('figures/tapir-example-with-data.png', plot = p, height = 11,
       width = 14, units = 'in', dpi = 600, bg = 'white')
