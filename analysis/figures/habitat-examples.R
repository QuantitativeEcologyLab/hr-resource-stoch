# plot(1, xlab = '\U1D53C') # run this before attaching ragg to plot unicode
# if notation doesn't show, see https://github.com/r-lib/ragg/issues/118
library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('purrr')     # for functional programming
library('ctmm')      # for movement modeling
library('lubridate') # for working with data more smoothly
library('ggplot2')   # for fancy plots
library('latex2exp') # for LaTeX in plots
library('cowplot')   # for fancy plots in grids
library('ragg')      # needed for custom alpha with coord_cartesian
# see https://github.com/tidyverse/ggplot2/issues/4029

source('functions/rgamma2.R') # function to generate rgamma() from mean and variance
source('functions/qgamma2.R') # function to generate qgamma() from mean and variance
source('analysis/figures/default-figure-styling.R')
theme_set(theme_get() + theme(legend.position = 'top', panel.grid = element_blank()))

e_r <- '\U1D53C(\U1D445)'
var_r <- '\U1D54D(\U1D445)'

# resource abundance palette
LOW <- '#744700'
MID <- '#d9bb94'
HIGH <- 'darkgreen'

# variable definition ----
# functions for trends in mean and variance
mean_fcn <- function(.t) {
  sinpi(.t / 4.5 - 0.75) * 20 + 35
}

var_fcn <- function(.t) {
  (sinpi(.t / 4.5 - 0.75) * 15 + 17.5)
}

params <- tibble(t = seq(1, 10, length.out = 250),
                 mu = mean_fcn(t),
                 sigma2 = var_fcn(t),
                 sigma2_constant = mean(sigma2))

lim_food <- c(0,
              qgamma2(p = 0.01, mu = max(params$mu), max(params$sigma2)) %>%
                ceiling())

lim_sigma2 <- c(floor(min(params$sigma2)), ceiling(max(params$sigma2)))

# changing mean, constant variance ----
# samples
d_1 <- expand_grid(t = 1:10,
                   patch = letters[t]) %>%
  mutate(mu = mean_fcn(t),
         sigma2 = mean(var_fcn(t)), # constant values
         food = rgamma2(mu = mu, sigma2 = sigma2, N = n()),
         e = food - mu)

# probability density of R
r_density <- tibble(r = seq(0, 10, length.out = 400),
                    dens = dgamma(r, shape = 3, scale = 1))

# a. raster of abundance with events
p_1a <-
  ggplot(d_1, aes(t, patch)) +
  geom_raster(aes(fill = food)) +
  geom_vline(aes(xintercept = t + 0.5)) +
  geom_point(color = 'black', size = 3) +
  geom_point(aes(color = sigma2), size = 2) +
  scale_x_continuous('Time', expand = c(0, 0), breaks = NULL) +
  scale_y_discrete('Patch', expand = c(0, 0)) +
  scale_fill_gradient2(e_r, low = LOW, mid = MID, high = HIGH,
                       limits = lim_food, breaks = lim_food, labels = c('Low', 'High'),
                       midpoint = quantile(lim_food, 0.5)) +
  scale_color_viridis_c(var_r, direction = -1, option = 'D', end = 1,
                        begin = 0.3, limits = lim_sigma2, breaks = lim_sigma2,
                        labels = c('Low', 'High')) +
  guides(color = guide_colorbar(order = 0), fill = guide_colorbar(order = 1))

# b. variable definition
# b-1: Gamma distribution
plot_expr <- function(string) {
  ggplot() +
    annotate('text', 0, 0, label = string, size = 8) +
    theme_void() +
    theme(text = element_text(family = 'Symbola'))
}

# probability density of R
r_pdf <-
  ggplot(r_density, aes(r, dens)) +
  geom_area(fill = pal[3], alpha = 0.5, color = pal[3]) +
  scale_x_continuous('\U1D445', breaks = NULL) +
  scale_y_continuous('P(\U1D445 = \U1D45F)', breaks = NULL)

## b-2, b-3: mean and variance trends
p_mean <-
  ggplot() +
  geom_point(aes(t, food), d_1, alpha = 0.25, show.legend = FALSE) +
  geom_line(aes(t, mu, color = mu), params, lwd = 1, show.legend = FALSE) +
  scale_x_continuous('Time', breaks = NULL) +
  scale_y_continuous('\U1D707(t)', limits = lim_food, breaks = lim_food,
                     labels = c('Low', 'High')) +
  scale_color_gradient2(low = LOW, mid = MID, high = HIGH, limits = lim_food,
                        midpoint = quantile(lim_food, 0.5))

p_var <-
  ggplot() +
  geom_line(aes(t, sigma2_constant, color = sigma2_constant), params, lwd = 1,
            show.legend = FALSE) +
  scale_x_continuous('Time', breaks = NULL) +
  scale_y_continuous('\U1D70E\U00B2', limits = lim_sigma2,
                     breaks = lim_sigma2, labels = c('Low', 'High')) +
  scale_color_viridis_c(direction = -1, option = 'D', end = 1,
                        begin = 0.3, limits = lim_sigma2)

p_1b <-
  plot_grid(
    plot_expr('\U1D445 \U007E \U0393(\U1D707(t), \U1D70E\U00B2)'), r_pdf,
    plot_expr('\U1D53C(\U1D445) = \U1D707(t)'), p_mean,
    plot_expr('\U1D54D(\U1D445) = \U1D70E\U00B2'), p_var,
    ncol = 2); p_1b

p_1 <- plot_grid(p_1a,
                 plot_grid(p_1b,
                           nrow = 1, labels = c('b.', '')),
                 labels = c('a.', ''), ncol = 1)

ggsave('figures/habitat-examples-constant-variance.png',
       plot = p_1, width = 6, height = 8, bg ='white')

# varying mean, varying variance ----
# samples
set.seed(4) # to ensure consistent samples within raster fill scale of R
d_2 <- expand_grid(t = 1:10,
                   patch = letters[t]) %>%
  mutate(mu = mean_fcn(t),
         sigma2 = var_fcn(t), # changing values
         food = rgamma2(mu = mu, sigma2 = sigma2, N = n()),
         e = food - mu)

# a. raster of abundance with events
p_2a <-
  ggplot(d_2, aes(t, patch)) +
  geom_raster(aes(fill = food)) +
  geom_vline(aes(xintercept = t + 0.5)) +
  geom_point(color = 'black', size = 3) +
  geom_point(aes(color = sigma2), size = 2) +
  scale_x_continuous('Time', expand = c(0, 0), breaks = NULL) +
  scale_y_discrete('Patch', expand = c(0, 0)) +
  scale_color_viridis_c(var_r, direction = -1, option = 'D', end = 1,
                        begin = 0.3, limits = lim_sigma2, breaks = lim_sigma2,
                        labels = c('Low', 'High')) +
  scale_fill_gradient2(e_r, low = LOW, mid = MID, high = HIGH,
                       limits = lim_food, breaks = lim_food, labels = c('Low', 'High'),
                       midpoint = quantile(lim_food, 0.5)) +
  guides(color = guide_colorbar(order = 0), fill = guide_colorbar(order = 1))

# b. variable definition
## b-2, b-3: mean and variance trends
p_mean <-
  ggplot() +
  geom_point(aes(t, food), d_2, alpha = 0.25, show.legend = FALSE) +
  geom_line(aes(t, mu, color = mu), params, lwd = 1, show.legend = FALSE) +
  scale_x_continuous('Time', breaks = NULL) +
  scale_y_continuous('\U1D707(t)', limits = lim_food, breaks = lim_food,
                     labels = c('Low', 'High')) +
  scale_color_gradient2(low = LOW, mid = MID, high = HIGH, limits = lim_food,
                        midpoint = quantile(lim_food, 0.5))

p_var <-
  ggplot() +
  geom_line(aes(t, sigma2, color = sigma2), params, lwd = 1,
            show.legend = FALSE) +
  scale_x_continuous('Time', breaks = NULL) +
  scale_y_continuous('\U1D70E\U00B2(t)', limits = lim_sigma2,
                     breaks = lim_sigma2, labels = c('Low', 'High')) +
  scale_color_viridis_c(direction = -1, option = 'D', end = 1,
                        begin = 0.3, limits = lim_sigma2)

p_2b <-
  plot_grid(
    plot_expr('\U1D445 \U007E \U0393(\U1D707(t), \U1D70E\U000B2(t))'), r_pdf,
    plot_expr('\U1D53C(\U1D445) = \U1D707(t)'), p_mean,
    plot_expr('\U1D54D(\U1D445) = \U1D70E\U000B2(t)'), p_var,
    ncol = 2)

p_2 <-
  plot_grid(
    p_2a,
    plot_grid(p_2b,
              nrow = 1, labels = c('b.', NA)),
    labels = c('a.', NA), ncol = 1)

ggsave('figures/habitat-examples-changing-variance.png',
       plot = p_2, width = 6, height = 8, bg ='white')
