library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('mgcv')    # for empirical Bayes GAMs
library('ggplot2') # for fancy plots
library('cowplot') # for fancy multi-panel plots

# for common theme and color palette
source('analysis/figures/barga-italy-2023/default-figure-styling-barga-2023.R')
theme_set(theme_get() + theme(legend.position = 'none'))

e_r <- 'Resource abundance, \U1D53C(\U1D445)'
v_r <- 'Resource stochasticity, \U1D54D(\U1D445)'

sims <- readRDS('simulations/days-hrs.rds') %>%
  filter(! mean %in% c('constant', 'erratic'),
         ! variance %in% c('constant', 'erratic')) %>%
  tibble() %>%
  pivot_longer(c(hr_50, hr_95), names_to = 'quantile', values_to = 'hr') %>%
  filter(quantile == 'hr_95')

m <- gam(
  list(
    # predictor for the mean
    hr ~ mu + sigma2, # interactions are flat
    # predictor for the scale (variance = scale * mean)
    ~ mu + sigma2), # interactions are flat
  family = gammals(),
  data = sims)

plot(m, all.terms = TRUE, pages = 1)

#' *NOTE*: using a GLM instead of a GAM because sampling of the space is
#' autocorrelated and not uniform
ggplot(sims, aes(mu, sigma2)) +
  geom_point(alpha = 0.3)

# estimate predictions ----
# effect of mu
newd_mu <- expand_grid(mu = seq(min(sims$mu), max(sims$mu),
                                length.out = 400),
                       sigma2 = mean(sims$sigma2),
                       quantile = unique(sims$quantile))

preds_mu <-
  bind_cols(
    newd_mu,
    predict(m, newdata = newd_mu, type = 'response', se.fit = TRUE) %>%
      data.frame() %>%
      transmute(hr = fit.1,
                var_hr = fit.1 * exp(fit.2)))

# effect of sigma2
newd_sigma2 <-
  expand_grid(mu = mean(sims$mu),
              sigma2 = seq(min(sims$sigma2), max(sims$sigma2),
                           length.out = 400),
              quantile = unique(sims$quantile))

preds_sigma2 <-
  bind_cols(
    newd_sigma2,
    predict(m, newdata = newd_sigma2, type = 'response', se.fit = TRUE) %>%
      data.frame() %>%
      transmute(hr = fit.1,
                var_hr = fit.1 * exp(fit.2)))

# final figure
preds <- bind_rows(mutate(preds_mu, x = e_r) %>%
                     rename(value = mu) %>%
                     select(-sigma2),
                   mutate(preds_sigma2, x = v_r) %>%
                     rename(value = sigma2) %>%
                     select(-mu))

sims_l <- sims %>%
  mutate(e = hr - predict(m, newdata = sims, type = 'response')[, 1]) %>%
  pivot_longer(c(mu, sigma2)) %>%
  mutate(x = if_else(name == 'mu', e_r, v_r))

ggplot() +
  facet_grid(. ~ x, scales = 'free', switch = 'x') +
  geom_point(aes(value, hr), sims_l, alpha = 0.5, color = pal[3]) +
  geom_line(aes(value, hr, group = quantile, color = x), preds,
            linewidth = 2) +
  scale_color_manual(values = pal) +
  scale_x_continuous(NULL, breaks = NULL) + 
  scale_y_continuous('Home range size', breaks = NULL) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 22))

ggsave('figures/2023-GRS-GRC-barga-italy/simulation-regression-plots.png',
       width = 48, height = 17.5, dpi = 'print', bg = 'white', units = 'cm')
