library('ggplot2') # for fancy plots
library('cowplot') # for multi-panel fancy plots
source('analysis/mean-variance-trends-panel-data.R') # tibble of mu, var

# for common theme and color palette
source('analysis/figures/barga-italy-2023/default-figure-styling-barga-2023.R')

# import simulation data
days <- readRDS('simulations/days-hrs.rds')

# main simulation figure of change in home range (5x5 plot) ----
e_r <- 'Resource abundance, \U1D53C(\U1D445)' # E(R) with blackboard bold
v_r <- 'Resource stochasticity, \U1D54D(\U1D445)' # V(R) w blackboard bold
hr_lab <- 'Home range size' # label for y axis
p_sim <-
  ggplot(days) +
  facet_grid(variance ~ mean, drop = FALSE) + # facet by trends in mu and var
  # add lines for 50% and 95% HRs
  geom_line(aes(animal, hr_95), color = pal[3], linewidth = 2) +
  scale_x_continuous('Time', breaks = NULL) +
  scale_y_continuous(hr_lab, breaks = NULL, limits = c(0, NA)) +
  # remove facet names since they appear in the marginal panels
  theme(strip.background = element_blank(), strip.text = element_blank())

# trends in the mean
p_mean <-
  ggplot(d55, aes(animal, mu)) +
  facet_grid(. ~ mean) +
  geom_line(color = pal[1], linewidth = 2) +
  scale_x_continuous(e_r, breaks = NULL, position = 'top') +
  scale_y_continuous(hr_lab, breaks = NULL) +
  theme(strip.background = element_blank(),
        axis.title.y = element_text(colour = 'transparent'))

# trends in the variance
p_variance <-
  ggplot(d55, aes(animal, sigma2)) +
  facet_grid(variance ~ ., switch = 'y') +
  geom_line(color = pal[2], linewidth = 2) +
  scale_x_continuous('Time', breaks = NULL) +
  scale_y_continuous(v_r, breaks = NULL) +
  theme(strip.background = element_blank(),
        axis.title.x = element_text(color = 'transparent'))

# create a single plot
p_full <- plot_grid(NULL, p_mean, p_variance, p_sim, rel_widths = c(1, 4.5),
                    rel_heights = c(1, 3.5), nrow = 2)
p_full

# save the plot as a png
ggsave('figures/2023-GRS-GRC-barga-italy/mean-variance-5-by-5-hr-sims.png',
       p_full, width = 72, height = 37, units = 'cm', dpi = 600,
       bg = 'white')
