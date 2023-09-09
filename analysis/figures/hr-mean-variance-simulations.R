library('ggplot2') # for fancy plots
library('cowplot') # for multi-panel fancy plots
source('analysis/figures/default-figure-styling.R') # defaults for figures
source('analysis/mean-variance-trends-panel-data.R') # tibble of mu, var

# import simulation data
days <- readRDS('simulations/days-hrs.rds')

# main simulation figure of change in home range (5x5 plot) ----
e_r <- 'Resource abundance, E(\U1D445|t) = \U1D707(\U1D461)'
v_r <- 'Resource stochasticity, Var(\U1D445|t) = \U1D70E\U00B2(\U1D461)'
hr_lab <- 'Space-use requirements, \U1D43B' # label for y axis
p_sim <-
  ggplot(days) +
  facet_grid(variance ~ mean, drop = FALSE) + # facet by trends in mu and var
  # add lines for 95% HR
  geom_line(aes(animal, hr_95), color = pal[3], linewidth = 1) +
  scale_x_continuous('Time', breaks = NULL) +
  scale_y_continuous(hr_lab, breaks = NULL, limits = c(0, NA)) +
  # remove facet names since they appear in the marginal panels
  theme(strip.background = element_blank(), strip.text = element_blank())

# trends in the mean
p_mean <-
  ggplot(d55, aes(animal, mu)) +
  facet_grid(. ~ mean) +
  geom_line(color = pal[1], linewidth = 1) +
  scale_x_continuous(e_r, breaks = NULL, position = 'top') +
  scale_y_continuous(hr_lab, breaks = NULL) +
  theme(strip.background = element_blank(),
        axis.title.y = element_text(colour = 'transparent'))

# trends in the variance
p_variance <-
  ggplot(d55, aes(animal, sigma2)) +
  facet_grid(variance ~ ., switch = 'y') +
  geom_line(color = pal[2], linewidth = 1) +
  scale_x_continuous('Time', breaks = NULL) +
  scale_y_continuous(v_r, breaks = NULL) +
  theme(strip.background = element_blank(),
        axis.title.x = element_text(color = 'transparent'))

# create a single plot
p_full <- plot_grid(NULL, p_mean, p_variance, p_sim, rel_widths = c(1, 4.5),
                    rel_heights = c(1, 3.5), nrow = 2)
p_full

# save the plot as a png
ggsave('figures/mean-variance-5-by-5-hr-sims.png', p_full,
       width = 10, height = 5.625, units = 'in', dpi = 'print', bg = 'white')
