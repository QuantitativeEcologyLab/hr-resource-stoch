library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('mgcv')    # for empirical Bayes GAMs
library('ggplot2') # for fancy plots
library('cowplot') # for fancy multi-panel plots
source('analysis/figures/default-figure-styling.R') # for color palette
theme_set(theme_get() +
            theme(legend.position = 'none',
                  strip.background = element_blank(),
                  strip.text = element_text(size = 14, face = 'bold')))

e_r <- 'paste(bold("Resource abundance, E("), bolditalic("R"), bold(")"))'
v_r <- 'paste(bold("Resource stochasticity, Var("), bolditalic("R"), bold(")"))'

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

# check if CV = sigma/mu is better than sigma
m_4 <- gam(
  list(
    hr ~ quantile + s(mu) + s(cv) + ti(mu, cv, k = 5),
    ~ quantile + s(mu) + s(cv) + ti(mu, cv, k = 5)),
  family = gammals(),
  data = mutate(sims, cv = sqrt(sigma2) / mu))

# ti term of the scale parameter looks over-fit and too uncertain
plot(m_4, all.terms = TRUE, pages = 1, scheme = 3, scale = 0)

AIC(m_1, m_2, m_3, m_4) # smooth model with interactions is best
BIC(m_1, m_2, m_3, m_4) # smooth model with interactions is best

summary(m_3)
summary(m_4)

# choose the best model
m <- m_3

plot(m, all.terms = TRUE, pages = 1, scheme = 3, scale = 0)

# estimate predictions ----
# effect of mu
newd_mu <-
  expand_grid(mu = seq(min(sims$mu), max(sims$mu), length.out = 400),
              sigma2 = mean(sims$sigma2),
              quantile = 'hr_95')

preds_mu <-
  bind_cols( # bind new data and predictions
    newd_mu,
    predict(m, newdata = newd_mu, type = 'link', se.fit = TRUE) %>%
      data.frame() %>% # convert list to data frame
      # 95% CIs assuming Gaussian credible intervals on the link scale
      transmute(hr = exp(fit.1),
                hr_lwr = exp(fit.1 - se.fit.1 * 1.96),
                hr_upr = exp(fit.1 + se.fit.1 * 1.96)))

p_mu <-
  ggplot() +
  geom_point(aes(mu, hr), filter(sims, quantile == 'hr_95'),
             alpha = 0.1, color = pal[3]) +
  geom_ribbon(aes(mu, ymin = hr_lwr, ymax = hr_upr), fill = pal[1],
              preds_mu, alpha = 0.5) +
  geom_line(aes(mu, hr), color = pal[1], linewidth = 1, preds_mu) +
  scale_x_continuous(
    bquote(paste(bold('Resource abundance, E('),
                 bolditalic('R'), bold(')'))), breaks = NULL) +
  scale_y_continuous(bquote(paste(bold('Space-use requirements, '),
                                  bolditalic('H'))), breaks = NULL)

# effect of sigma2
newd_sigma2 <-
  expand_grid(mu = mean(sims$mu),
              sigma2 = seq(min(sims$sigma2), max(sims$sigma2),
                           length.out = 400),
              quantile = 'hr_95')

preds_sigma2 <-
  bind_cols(
    newd_sigma2,
    predict(m, newdata = newd_sigma2, type = 'link', se.fit = TRUE) %>%
      data.frame() %>%
      transmute(hr = exp(fit.1),
                hr_lwr = exp(fit.1 - se.fit.1 * 1.96),
                hr_upr = exp(fit.1 + se.fit.1 * 1.96)))

p_sigma2 <-
  ggplot() +
  geom_point(aes(sigma2, hr), filter(sims, quantile == 'hr_95'),
             alpha = 0.1, color = pal[3]) +
  geom_ribbon(aes(sigma2, ymin = hr_lwr, ymax = hr_upr), fill = pal[2],
              preds_sigma2, alpha = 0.5) +
  geom_line(aes(sigma2, hr), color = pal[2], linewidth = 1,
            preds_sigma2) +
  scale_x_continuous(
    bquote(paste(bold('Resource stochasticity, Var('),
                 bolditalic('R'), bold(')'))), breaks = NULL) +
  scale_y_continuous(bquote(paste(bold('Space-use requirements, '),
                                  bolditalic('H'))), breaks = NULL)

# full interaction plot
newd_mu_sigma2 <-
  expand_grid(mu = seq(min(sims$mu), max(sims$mu), length.out = 400),
              sigma2 = seq(min(sims$sigma2), max(sims$sigma2),
                           length.out = 400),
              quantile = 'hr_95')

preds_mu_sigma2 <-
  bind_cols(newd_mu_sigma2,
            predict(m, newdata = newd_mu_sigma2, type = 'response',
                    se.fit = TRUE, terms = c('ti(mu,sigma2)',
                                             'ti(mu,sigma2).1')) %>%
              as.data.frame()) %>%
  rename(hr = fit.1)

p_mu_sigma2 <-
  ggplot(preds_mu_sigma2) +
  geom_raster(aes(mu, sigma2, fill = hr)) +
  geom_contour(aes(mu, sigma2, z = hr), color = 'black') +
  scale_x_continuous(
    bquote(paste(bold('Resource abundance, E('),
                 bolditalic('R'), bold(')'))), breaks = NULL,
    expand = c(0, 0)) +
  scale_y_continuous(
    bquote(paste(bold('Resource stochasticity, Var('),
                 bolditalic('R'), bold(')'))), breaks = NULL,
    expand = c(0, 0)) +
  scale_fill_gradient(bquote(paste(atop(bold('Space-use'),
                                        paste(bold('requirements, '),
                                              bolditalic('H  '))))),
                      breaks = range(preds_mu_sigma2$hr),
                      labels = c('Low', 'High'),
                      low = 'grey90', high = pal[3]) +
  theme(legend.position = 'top')

# final figure
cowplot::plot_grid(p_mu, p_sigma2, p_mu_sigma2, labels = 'AUTO',
                   nrow = 1)

ggsave('figures/simulation-regression-plots.png',
       width = 9, height = 3, dpi = 'print', bg = 'white', scale = 1.5)

# effect of sqrt(sigma2) is nonlinear since it is on the same scale as mu
sd_r <- 'paste(bold("Standard deviation in "), bolditalic("R"))'

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
  facet_grid(. ~ x, scales = 'free', switch = 'x',
             labeller = label_parsed) + # label below x axis
  geom_point(aes(value, hr), sims_sd, alpha = 0.1, color = pal[3]) +
  geom_line(aes(value, hr, group = quantile, color = x), preds_sd,
            linewidth = 1) +
  scale_color_manual(values = pal[-(3:6)]) +
  scale_x_continuous(NULL, breaks = NULL) + 
  scale_y_continuous(bquote(paste(bold('Space-use requirements, '),
                                  bolditalic('H'))), breaks = NULL) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 12))
