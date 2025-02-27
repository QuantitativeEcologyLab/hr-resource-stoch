library('ggplot2') # for fancy plots
library('cowplot') # for multi-panel fancy plots
source('analysis/figures/default-figure-styling.R') # defaults for figures
source('analysis/simulations/mean-variance-trends-panel-data.R')

# import simulation data
days <- readRDS('simulations/days-hrs.rds')

# main simulation figure of change in home range (5x5 plot) ----
e_r <- bquote(paste(bold('Resource abundance, E('), bolditalic('R'),
                    bold(') = '), '\U1D6CD', bold('('), bolditalic('t'), bold(')')))
v_r <- bquote(paste(bold('Resource stochasticity, Var('), bolditalic('R'),
                    bold(') = '), '\U1D6D4\U00B2', bold('('),
                    bolditalic('t'), bold(')')))
hr_lab <- bquote(paste(bold('Home-range size, '), bolditalic('H')))

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
p_full <- plot_grid(NULL, p_mean, p_variance, p_sim, rel_widths = c(1, 4),
                    rel_heights = c(1, 3.1), nrow = 2)
p_full

# save the plot as a png
ggsave('figures/mean-variance-5-by-5-hr-sims.png', p_full,
       width = FULL_WIDTH, height = 100, units = 'mm', scale = 1.5,
       dpi = 600, bg = 'white')
