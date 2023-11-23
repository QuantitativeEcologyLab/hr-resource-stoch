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
e_r <- bquote(paste(bold('Resource abundance, '), '\U1D6CD', bold('('),
                    bolditalic('t'), bold(')')))
v_r <- bquote(paste(bold('Resource stochasticity, '), '\U1D6D4\U00B2',
                    bold('('), bolditalic('t'), bold(')')))
hr_lab <- bquote(paste(bold('7-day home range size (km'), '\U00B2',
                       bold(')')))

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

if(! file.exists('data/anna-hr-ndvi-data.rds')) {
  tapir <-
    readRDS('models/tapirs/CE_31_ANNA-window-7-days-dt-1-days.rds') %>%
    mutate(est = map(dataset, \(.d) filter(tel, timestamp %in% .d$timestamp)),
           mu = map_dbl(est, \(.d) mean(.d$mu)),
           sigma2 = map_dbl(est, \(.d) mean(.d$sigma2))) %>%
    select(date, mu, sigma2, hr_lwr_95, hr_est_95, hr_upr_95) %>%
    mutate(weight = 1 / (hr_upr_95 - hr_lwr_95), # weights based on CI width
           weight = weight / mean(weight)) # standardize weights
  sum(tapir$weight) == nrow(tapir) # ensure weights are right
  saveRDS(tapir, 'data/anna-hr-ndvi-data.rds')
} else {
  tapir <- readRDS('data/anna-hr-ndvi-data.rds')
}

# create the figure
date_labs <- range(tel$timestamp) %>% as.Date()

# E(R) and V(R) are highly correlated
ggplot(tapir) +
  geom_point(aes(mu, sigma2, alpha = weight)) +
  labs(x = '\U1D707(t)', y = '\U1D70E\U00B2(t)') +
  theme(text = element_text(face = 'plain'))

l_grobs <-
  lapply(
    list(
      # mean, variance and HR
      ggplot(tapir, aes(date, mu)) +
        geom_line(color = pal[1], linewidth = 2) +
        labs(x = NULL, y = e_r),
      ggplot(tapir, aes(date, sigma2)) +
        geom_line(color = pal[2], linewidth = 2) +
        labs(x = NULL, y = v_r),
      ggplot(tapir) +
        geom_line(aes(date, hr_est_95), color = pal[3], linewidth = 2) +
        labs(x = NULL, y = hr_lab)),
    as_grob) # convert to grid graphical objects (grobs)

# align left margins of all plots
aligned_widths <- align_margin(map(l_grobs, function(x) {x$widths}), 'first')

# Setting the dimensions of plots to the aligned dimensions
for (i in seq_along(l_grobs)) {
  l_grobs[[i]]$widths <- aligned_widths[[i]]
}

# Draw aligned plots
p_left <- plot_grid(plotlist = l_grobs, ncol = 1, labels = 'AUTO',
                    label_size = 22, label_x = -0.01, label_y = 1.025)

# regression plots ----
m <- gam(hr_est_95 ~ s(mu, bs = 'ts', k = 5) + s(sigma2, bs = 'ts', k = 5),
         family = Gamma('log'),
         weights = weight, # weights based on CI size
         data = tapir,
         method = 'REML')
# trends should be monotonic...
plot(m, pages = 1, scheme = 1)

# using shape-constrained splines to avoid artifacts due to correlation
m <- scam(hr_est_95 ~
            s(mu, bs = 'mpd', k = 5) +
            s(sigma2, bs = 'mpi', k = 5),
          family = Gamma('log'),
          weights = weight,
          data = tapir)
plot(m, pages = 1, scheme = 1)

m0 <- scam(hr_est_95 ~ s(mu, bs = 'mpd', k = 5),
           family = Gamma('log'),
           weights = weight,
           data = tapir)
AIC(m, m0) # model that accounts for both mu and sigma2 fits much better

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
  geom_point(aes(mu, hr_est_95, alpha = weight), tapir, color = pal[3]) +
  geom_ribbon(aes(mu, ymin = hr_mu_lwr, ymax = hr_mu_upr), preds,
              fill = pal[1], alpha = 0.3) +
  geom_line(aes(mu, hr_mu_est), preds, color = pal[1], linewidth = 2) +
  labs(x = e_r, y = hr_lab) +
  scale_alpha_continuous('Weight', range = c(0.3, 1)) +
  theme(legend.position = c(0.93, 0.86),
        legend.box.background = element_rect(color = '#00000080'),
        legend.background = element_blank()) +
  guides(alpha = guide_legend(override.aes = list(color = 'black')))

p_e <- ggplot() +
  coord_cartesian(ylim = c(0, 12.5)) +
  geom_point(aes(sigma2, hr_est_95, alpha = weight), tapir, color = pal[3]) +
  geom_ribbon(aes(sigma2, ymin = hr_sigma2_lwr, ymax = hr_sigma2_upr), preds,
              fill = pal[2], alpha = 0.3) +
  geom_line(aes(sigma2, hr_sigma2_est), preds, color = pal[2], linewidth = 2) +
  labs(x = v_r, y = hr_lab) +
  scale_alpha_continuous(range = c(0.3, 1))

p_f <- ggplot(tapir) +
  geom_point(aes(mu, sigma2, alpha = weight)) +
  labs(x = e_r, y = v_r) +
  scale_alpha_continuous(range = c(0.3, 1))

r_grobs <- map(list(p_d, p_e, p_f), as_grob)

# align right margins of all plots
aligned_widths <- align_margin(map(r_grobs, function(x) {x$widths}), 'first')

# Setting the dimensions of plots to the aligned dimensions
for (i in seq_along(r_grobs)) {
  r_grobs[[i]]$widths <- aligned_widths[[i]]
}

p_right <-
  plot_grid(plotlist = r_grobs, ncol = 1, labels = c('D', 'E', 'F'),
            label_size = 22, label_x = -0.01, label_y = 1.025)

# add space to avoid cutting off x axis for (e.)
p <- plot_grid(p_left, p_right, NULL, nrow = 1, rel_widths = c(1, 1, 0.01))
p

ggsave('figures/tapir-example-with-data.png', plot = p, height = 12,
       width = 14, units = 'in', dpi = 600, bg = 'white')
