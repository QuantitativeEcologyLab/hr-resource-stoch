library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('mgcv')    # for empirical Bayes GAMs
library('ggplot2') # for fancy plots
library('cowplot') # for fancy multi-panel plots
source('analysis/figures/default-figure-styling.R') # for color palette
theme_set(theme_bw() + theme(legend.position = 'none'))

e_r <- 'Resource abundance, E(\U1D445)'
v_r <- 'Resource stochasticity, Var(\U1D445)'

sims <- readRDS('simulations/days-hrs.rds') %>%
  tibble() %>%
  pivot_longer(c(hr_50, hr_95), names_to = 'quantile', values_to = 'hr') %>%
  mutate(quantile = factor(quantile))

# sampling is not uniform, but not highly problematic
# removing the data with constant mu or sigma2 has no effect on the model
ggplot(sims, aes(mu, sigma2)) +
  geom_point(alpha = 0.3)

m_1 <- gam(
  list(
    # predictor for the mean
    hr ~ quantile + mu + sigma2,
    # predictor for the scale (variance = scale * mean)
    ~ quantile + mu + sigma2),
  family = gammals(),
  data = sims)

m_2 <- gam(
  list(
    hr ~ quantile + s(mu) + s(sigma2),
    ~ quantile + s(mu) + s(sigma2)),
  family = gammals(),
  data = sims)

m_3 <- gam(
  list(
    hr ~ quantile + s(mu) + s(sigma2) + ti(mu, sigma2, k = 5),
    ~ quantile + s(mu) + s(sigma2) + ti(mu, sigma2, k = 5)),
  family = gammals(),
  data = sims)

AIC(m_1, m_2, m_3) # smooth model with interactions is best
BIC(m_1, m_2, m_3) # smooth model with interactions is best

m <- m_3

plot(m, all.terms = TRUE, pages = 1, scheme = 3, scale = 0)

# estimate predictions ----
# effect of mu
newd_mu <-
  expand_grid(mu = seq(min(sims$mu), max(sims$mu), length.out = 400),
              sigma2 = mean(sims$sigma2),
              quantile = unique(sims$quantile))

preds_mu <-
  bind_cols( # bind new data and predictions
    newd_mu,
    predict(m, newdata = newd_mu, type = 'link', se.fit = TRUE) %>%
      data.frame() %>% # convert list to data frame
      # 95% CIs assuming Gaussian credible intervals on the link scale
      transmute(hr = exp(fit.1),
                hr_lwr = exp(fit.1 - se.fit.1 * 1.96),
                hr_upr = exp(fit.1 + se.fit.1 * 1.96)))


# effect of sigma2
newd_sigma2 <-
  expand_grid(mu = mean(sims$mu),
              sigma2 = seq(min(sims$sigma2), max(sims$sigma2),
                           length.out = 400),
              quantile = unique(sims$quantile))

preds_sigma2 <-
  bind_cols(
    newd_sigma2,
    predict(m, newdata = newd_sigma2, type = 'link', se.fit = TRUE) %>%
      data.frame() %>%
      transmute(hr = exp(fit.1),
                hr_lwr = exp(fit.1 - se.fit.1 * 1.96),
                hr_upr = exp(fit.1 + se.fit.1 * 1.96)))

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
  facet_grid(. ~ x, scales = 'free', switch = 'x') + # label below x axis
  geom_point(aes(value, hr), filter(sims_l, quantile == 'hr_95'),
             alpha = 0.1, color = pal[3]) +
  geom_ribbon(aes(value, ymin = hr_lwr, ymax = hr_upr, fill = x),
              filter(preds, quantile == 'hr_95'), alpha = 0.5) +
  geom_line(aes(value, hr, color = x), linewidth = 1,
            filter(preds, quantile == 'hr_95')) +
  scale_color_manual(values = pal) +
  scale_x_continuous(NULL, breaks = NULL) + 
  scale_y_continuous('Space-use requirements, \U1D43B', breaks = NULL) +
  # to make facet strip look like the x axis
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_text(size = 12))

ggsave('figures/simulation-regression-plots.png',
       width = 6, height = 3, dpi = 'print', bg = 'white', scale = 1.5)

# effect of sqrt(sigma2) is nonlinear since it is on the same scale as mu
sd_r <- "Standard deviation in \U1D445, \U221A Var(\U1D445)"

preds_sd <- bind_rows(preds,
                      filter(preds, x == v_r) %>%
                        mutate(value = sqrt(value),
                               x = sd_r)) %>%
  filter(quantile == 'hr_95')

sims_sd <- sims %>%
  filter(quantile == 'hr_95') %>%
  mutate(sigma = sqrt(sigma2)) %>%
  pivot_longer(c(mu, sigma2, sigma)) %>%
  mutate(x = case_when(name == 'mu' ~ e_r,
                       name == 'sigma2' ~ v_r,
                       name == 'sigma' ~ sd_r))

ggplot() +
  facet_grid(. ~ x, scales = 'free', switch = 'x') +
  geom_point(aes(value, hr), sims_sd, alpha = 0.1, color = pal[3]) +
  geom_line(aes(value, hr, group = quantile, color = x), preds_sd,
            linewidth = 1) +
  scale_color_manual(values = pal[-(3:6)]) +
  scale_x_continuous(NULL, breaks = NULL) + 
  scale_y_continuous('Home range size, \U1D43B', breaks = NULL) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12))
