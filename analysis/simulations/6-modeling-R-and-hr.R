library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('mgcv')    # for empirical Bayes GAMs
library('gratia')  # for empirical Bayes GAMs
library('ggplot2') # for fancy plots
library('cowplot') # for fancy multi-panel plots
source('analysis/figures/default-figure-styling.R') # for color palette
source('functions/gammals-variance-simulation-and-derivatives.R')
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
  data = sims,
  method = 'REML')
plot(m_1, pages = 1, all.terms = TRUE) # acceptable but likely too rigid

m_2 <- gam(
  list(
    hr ~ quantile + s(mu) + s(sigma2),
    ~ quantile + s(mu) + s(sigma2)),
  family = gammals(),
  data = sims,
  method = 'REML')
plot(m_2, pages = 1, all.terms = TRUE) # scale terms are inappropriate

m_3 <- gam(
  list(
    hr ~ quantile + s(mu, k = 5) + s(sigma2, k = 5) + ti(mu, sigma2, k = 5),
    ~ quantile + s(mu, k = 5) + s(sigma2, k = 5) + ti(mu, sigma2, k = 5)),
  family = gammals(),
  data = sims,
  method = 'REML')
plot(m_3, pages = 1, all.terms = TRUE, scheme = 3) # much better

# check if CV = sigma/mu is better than sigma
m_4 <- gam(
  list(
    hr ~ quantile + s(mu) + s(cv) + ti(mu, cv, k = 5),
    ~ quantile + s(mu) + s(cv) + ti(mu, cv, k = 5)),
  family = gammals(),
  data = mutate(sims, cv = sqrt(sigma2) / mu),
  method = 'REML')
# ti term of the scale parameter looks over-fit
plot(m_4, all.terms = TRUE, pages = 1, scheme = 3, scale = 0)

AIC(m_1, m_2, m_3, m_4) # smooth model with interactions is best
BIC(m_1, m_2, m_3, m_4) # smooth model with interactions is best

summary(m_3)
summary(m_4)

# choose the best model
m <- m_3

plot(m, all.terms = TRUE, pages = 1, scheme = 3, scale = 0)

# filter to 95% HR for simplicity
sims <- filter(sims, quantile == 'hr_95')

# predictions for calculating (squared) residuals
sims$hr_est_95 <- predict(m, newdata = sims, type = 'response',
                          se.fit = FALSE)[, 1]

# effects of mu ----
newd_mu <-
  expand_grid(mu = seq(min(sims$mu), max(sims$mu), length.out = 400),
              sigma2 = mean(sims$sigma2),
              quantile = 'hr_95')

preds_mu <-
  left_join(gammals_mean(m, newd_mu, nsims = 1e4, unconditional = TRUE),
            gammals_var(m, newd_mu, nsims = 1e4, unconditional = TRUE),
            by = c('row', 'mu', 'sigma2', 'quantile', 'simulation')) %>%
  group_by(row, mu, sigma2, quantile) %>%
  summarise(hr_mean = mean(mean),
            hr_mean_975 = quantile(mean, 0.025),
            hr_mean_025 = quantile(mean, 0.975),
            hr_var = mean(variance),
            hr_var_975 = quantile(variance, 0.025),
            hr_var_025 = quantile(variance, 0.975),
            .groups = 'drop')

p_mu_mean <-
  ggplot() +
  geom_point(aes(mu, hr), sims, alpha = 0.1, color = pal[3]) +
  geom_ribbon(aes(mu, ymin = hr_mean_025, ymax = hr_mean_975), fill = pal[1],
              preds_mu, alpha = 0.3) +
  geom_line(aes(mu, hr_mean), color = pal[1], linewidth = 1, preds_mu) +
  scale_x_continuous(
    bquote(paste(bold('Resource abundance, E('),
                 bolditalic('R'), bold(')'))), breaks = NULL) +
  scale_y_continuous(bquote(paste(bold('Space-use requirements, '),
                                  bolditalic('H'))), breaks = NULL)

p_mu_var <-
  ggplot() +
  # geom_point(aes(mu, (hr - hr_est_95)^2), sims, alpha = 0.1,
  #            color = pal[3]) +
  geom_ribbon(aes(mu, ymin = hr_var_025, ymax = hr_var_975), fill = pal[1],
              preds_mu, alpha = 0.3) +
  geom_line(aes(mu, hr_var), color = pal[1], linewidth = 1, preds_mu) +
  scale_x_continuous(
    bquote(paste(bold('Resource abundance, E('),
                 bolditalic('R'), bold(')'))), breaks = NULL) +
  scale_y_continuous(bquote(paste(
    bold('Variance in space-use requirements'))), breaks = NULL)

# effects of sigma2 ----
newd_sigma2 <-
  expand_grid(mu = mean(sims$mu),
              sigma2 = seq(min(sims$sigma2), max(sims$sigma2),
                           length.out = 400),
              quantile = 'hr_95')

preds_sigma2 <-
  left_join(gammals_mean(m, newd_sigma2, nsims = 1e4, unconditional = TRUE),
            gammals_var(m, newd_sigma2, nsims = 1e4, unconditional = TRUE),
            by = c('row', 'mu', 'sigma2', 'quantile', 'simulation')) %>%
  group_by(row, mu, sigma2, quantile) %>%
  summarise(hr_mean = mean(mean),
            hr_mean_975 = quantile(mean, 0.025),
            hr_mean_025 = quantile(mean, 0.975),
            hr_var = mean(variance),
            hr_var_975 = quantile(variance, 0.025),
            hr_var_025 = quantile(variance, 0.975),
            .groups = 'drop')

p_sigma2_mean <-
  ggplot() +
  geom_point(aes(sigma2, hr), sims, alpha = 0.1, color = pal[3]) +
  geom_ribbon(aes(sigma2, ymin = hr_mean_025, ymax = hr_mean_975, group = mu),
              fill = pal[2], preds_sigma2, alpha = 0.3) +
  geom_line(aes(sigma2, hr_mean, lty = factor(mu)),
            color = pal[2], linewidth = 1,
            preds_sigma2) +
  scale_x_continuous(
    bquote(paste(bold('Resource stochasticity, Var('),
                 bolditalic('R'), bold(')'))), breaks = NULL) +
  scale_y_continuous(bquote(paste(bold('Space-use requirements, '),
                                  bolditalic('H'))), breaks = NULL)

p_sigma2_var <-
  ggplot() +
  # geom_point(aes(sigma2, (hr - hr_est_95)^2), sims, alpha = 0.1,
  #            color = pal[3]) +
  geom_ribbon(aes(sigma2, ymin = hr_var_025, ymax = hr_var_975),
              fill = pal[2], preds_sigma2, alpha = 0.3) +
  geom_line(aes(sigma2, hr_var), color = pal[2], linewidth = 1,
            preds_sigma2) +
  scale_x_continuous(
    bquote(paste(bold('Resource stochasticity, Var('),
                 bolditalic('R'), bold(')'))), breaks = NULL) +
  scale_y_continuous(bquote(paste(
    bold('Variance in space-use requirements'))), breaks = NULL)

# full plots with interactions ----
newd_full <-
  expand_grid(mu = seq(min(sims$mu), max(sims$mu), length.out = 100),
              sigma2 = seq(min(sims$sigma2), max(sims$sigma2),
                           length.out = 100),
              quantile = 'hr_95')

preds_full <-
  left_join(gammals_mean(m, newd_full, nsims = 1e4, unconditional = TRUE),
            gammals_var(m, newd_full, nsims = 1e4, unconditional = TRUE),
            by = c('row', 'mu', 'sigma2', 'quantile', 'simulation')) %>%
  group_by(row, mu, sigma2, quantile) %>%
  summarise(hr_mean = mean(mean),
            hr_var = mean(variance),
            .groups = 'drop')

# mean
fill_mean_lab <- bquote(atop(bold('Mean space-use'),
                             paste(bold('requirements, E('),
                                   bolditalic('H'), bold(')  '))))

preds_full %>%
  filter(sigma2 %in% sort(unique(sigma2))[c(10, 30, 50, 70, 90)]) %>%
  ggplot() +
  geom_line(aes(mu, hr_mean, group = sigma2), color = pal[1]) +
  scale_x_continuous( bquote(paste(bold('Resource abundance, E('),
                                   bolditalic('R'), bold(')'))),
                      breaks = NULL) +
  scale_y_continuous(fill_mean_lab, breaks = NULL)

p_full_mean <-
  ggplot(preds_full) +
  geom_raster(aes(mu, sigma2, fill = hr_mean)) +
  geom_contour(aes(mu, sigma2, z = hr_mean), color = 'black') +
  scale_x_continuous(
    bquote(paste(bold('Resource abundance, E('),
                 bolditalic('R'), bold(')'))), breaks = NULL,
    expand = c(0, 0)) +
  scale_y_continuous(
    bquote(paste(bold('Resource stochasticity, Var('),
                 bolditalic('R'), bold(')'))), breaks = NULL,
    expand = c(0, 0)) +
  scale_fill_gradient(fill_mean_lab, low = 'grey90', high = pal[3], 
                       breaks = range(preds_full$hr_mean),
                       labels = c('low', 'high')) +
  theme(legend.position = 'top')

# variance
fill_var_lab <- bquote(atop(bold('Variance in space-use'),
                            paste(bold('requirements, Var('),
                                  bolditalic('H'), bold(')  '))))

preds_full %>%
  filter(mu %in% sort(unique(mu))[c(10, 30, 50, 70, 90)]) %>%
  ggplot() +
  geom_line(aes(sigma2, hr_mean, group = mu), color = pal[2]) +
  scale_x_continuous(bquote(paste(bold('Resource stochasticity, Var('),
                                  bolditalic('R'), bold(')'))),
                     breaks = NULL) +
  scale_y_continuous(fill_mean_lab, breaks = NULL)

p_full_var <-
  ggplot(preds_full) +
  geom_raster(aes(mu, sigma2, fill = hr_var)) +
  geom_contour(aes(mu, sigma2, z = hr_var), color = 'black') +
  scale_x_continuous(
    bquote(paste(bold('Resource abundance, E('),
                 bolditalic('R'), bold(')'))), breaks = NULL,
    expand = c(0, 0)) +
  scale_y_continuous(
    bquote(paste(bold('Resource stochasticity, Var('),
                 bolditalic('R'), bold(')'))), breaks = NULL,
    expand = c(0, 0)) +
  scale_fill_gradient(fill_var_lab, low = 'grey90', high = pal[7], 
                      breaks = range(preds_full$hr_var),
                      labels = c('low', 'high')) +
  theme(legend.position = 'top')

# final figures ----
# mean
cowplot::plot_grid(p_mu_mean, p_sigma2_mean, p_full_mean, labels = 'AUTO',
                   nrow = 1)

ggsave('figures/simulation-regression-mean.png',
       width = 9, height = 3, dpi = 'print', bg = 'white', scale = 1.5)

# variance
cowplot::plot_grid(p_mu_var, p_sigma2_var, p_full_var, labels = 'AUTO',
                   nrow = 1)

ggsave('figures/simulation-regression-var.png',
       width = 9, height = 3, dpi = 'print', bg = 'white', scale = 1.5)

# mean and variance together
cowplot::plot_grid(p_mu_mean, p_sigma2_mean, p_full_mean,
                   p_mu_var, p_sigma2_var, p_full_var,
                   labels = 'AUTO', nrow = 2)

ggsave('figures/simulation-regression-mean-and-var.png',
       width = 9, height = 6, dpi = 'print', bg = 'white', scale = 1.5)

# effect of sigma2 for a few values of mu ----
newd_mu_2 <-
  expand_grid(mu = seq(min(sims$mu), max(sims$mu), length.out = 200),
              sigma2 = round(quantile(sims$sigma2)),
              quantile = 'hr_95')

preds_mu_2 <-
  left_join(gammals_mean(m, newd_mu_2, nsims = 1e4, unconditional = TRUE),
            gammals_var(m, newd_mu_2, nsims = 1e4, unconditional = TRUE),
            by = c('row', 'mu', 'sigma2', 'quantile', 'simulation')) %>%
  group_by(row, mu, sigma2, quantile) %>%
  summarise(hr_mean = mean(mean),
            hr_mean_lwr = quantile(mean, 0.025),
            hr_mean_upr = quantile(mean, 0.975),
            hr_var = mean(variance),
            hr_var_lwr = quantile(variance, 0.025),
            hr_var_upr = quantile(variance, 0.975),
            .groups = 'drop')

ggplot(preds_mu_2) +
  geom_ribbon(aes(mu, ymin = hr_mean_lwr, ymax = hr_mean_upr,
                  fill = sigma2, group = sigma2), alpha = 0.2) +
  geom_line(aes(mu, hr_mean, color = sigma2, group = sigma2)) +
  scale_x_continuous(
    bquote(paste(bold('Resource abundance, E('),
                 bolditalic('R'), bold(')'))), breaks = NULL,
    expand = c(0, 0)) +
  scale_y_continuous(
    bquote(paste(bold('Resource stochasticity, Var('),
                 bolditalic('R'), bold(')'))), breaks = NULL,
    expand = c(0, 0)) +
  scale_color_viridis_c(bquote(paste(bold('Resource stochasticity, Var('),
                                     bolditalic('R'), bold(')'))),
                        limits = range(sims$sigma2), direction = -1,
                        breaks = range(sims$sigma2),
                        labels = c('low', 'high'),
                        aesthetics = c('color', 'fill')) +
  theme(legend.position = 'top')

# effect of SD(R) on E(H) ----
ggplot() +
  geom_point(aes(sqrt(sigma2), hr), sims, alpha = 0.1, color = pal[3]) +
  geom_ribbon(aes(sqrt(sigma2), ymin = hr_mean_025, ymax = hr_mean_975),
              fill = pal[2], preds_sigma2, alpha = 0.3) +
  geom_line(aes(sqrt(sigma2), hr_mean), color = pal[2], linewidth = 1,
            preds_sigma2) +
  scale_x_continuous(bquote(paste(bold('Standard deviation in '), bolditalic('R')))) +
  scale_y_continuous(bquote(paste(bold('Space-use requirements, '),
                                  bolditalic('H'))), breaks = NULL)

ggplot() +
  # geom_point(aes(sqrt(sigma2), (hr - hr_est_95)^2), sims, alpha = 0.1,
  #            color = pal[3]) +
  geom_ribbon(aes(sqrt(sigma2), ymin = hr_var_025, ymax = hr_var_975),
              fill = pal[2], preds_sigma2, alpha = 0.3) +
  geom_line(aes(sqrt(sigma2), hr_var), color = pal[2], linewidth = 1,
            preds_sigma2) +
  scale_x_continuous(
    bquote(paste(bold('Resource stochasticity, Var('),
                 bolditalic('R'), bold(')'))), breaks = NULL) +
  scale_y_continuous(bquote(paste(
    bold('Variance in space-use requirements, Var('),
    bolditalic('H'), bold(')'))), breaks = NULL)
