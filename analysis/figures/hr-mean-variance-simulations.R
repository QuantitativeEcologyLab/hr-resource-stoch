library('ggplot2') # for fancy plots
library('cowplot') # for multi-panel fancy plots
source('analysis/figures/default-figure-styling.R') # defaults for figures
source('analysis/mean-variance-trends-panel-data.R') # tibble of mu, var

# import simulation data
days <- readRDS('simulations/days-hrs.rds')

# main simulation figure of change in home range (5x5 plot) ----
e_r <- 'Resource abundance, \U1D53C(\U1D445|t) = \U1D707(\U1D461)'
v_r <- 'Resource stochasticity, \U1D54D(\U1D445|t) = \U1D70E\U00B2(\U1D461)'
hr_lab <- 'Home range size, \U1D43B' # label for y axis
p_sim <-
  ggplot(days) +
  facet_grid(variance ~ mean, drop = FALSE) + # facet by trends in mu and var
  # area between 50% and 95% HRs
  # geom_ribbon(aes(animal, ymin = hr_50, ymax = hr_95), fill = pal[3],
  #             alpha = 0.2) +
  # add lines for 50% and 95% HRs
  # geom_line(aes(animal, hr_50), color = pal[3]) +
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

# estimate fraction of animals not satisfied if only mean is accounted for ----
days_suff <-
  left_join(days,
            filter(days, variance == 'constant') %>%
              select(mean, animal, hr_50, hr_95) %>%
              rename(hr_50_x = hr_50, hr_95_x = hr_95),
            by = c('mean', 'animal')) %>%
  # see if the animal has enough space based on the 95% quantile
  mutate(enough_50 = hr_50 <= hr_95_x,
         enough_95 = hr_95 <= hr_95_x,
         additional = if_else(enough_95, hr_95, hr_95_x))

p_sim_suff <-
  ggplot(days_suff) +
  facet_grid(variance ~ mean, drop = FALSE) + # facet by trends in mu and var
  # area between 50% and 95% HRs
  geom_ribbon(aes(animal, ymin = hr_50_x, ymax = hr_95_x), fill = pal[3],
              alpha = 0.2) +
  # area missing to satisfy needs
  geom_ribbon(aes(animal, ymin = hr_95, ymax = additional), fill = pal[7],
              alpha = 0.2) +
  # add lines for 50% and 95% HRs
  geom_line(aes(animal, hr_50, color = enough_50, group = 1)) +
  geom_line(aes(animal, hr_95, color = enough_95, group = 1)) +
  scale_x_continuous('Time', breaks = NULL) +
  scale_y_continuous(hr_lab, breaks = NULL, limits = c(0, NA)) +
  scale_color_manual('Space is sufficient', values = pal[c(7, 3)],
                     labels = c('No', 'Yes'), aesthetics = c('color', 'fill')) +
  # remove facet names since they appear in the marginal panels
  theme(strip.background = element_blank(), strip.text = element_blank(),
        legend.position = 'top')

p_full_suff <-
  plot_grid(get_legend(p_sim_suff),
            plot_grid(NULL, p_mean, rel_widths = c(1, 4.4), nrow = 1),
            plot_grid(p_variance, p_sim_suff + theme(legend.position = 'none'),
                      rel_widths = c(1, 4.5), nrow = 1),
            rel_heights = c(0.2, 1, 4), nrow = 3)
p_full_suff
