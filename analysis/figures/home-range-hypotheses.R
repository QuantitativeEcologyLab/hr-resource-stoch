library('dplyr')   # for data wrangling
library('ggplot2') # for fancy plots
library('ggbrace') # for curly braces on fancy plots
library('cowplot') # for fancy multi-panel plots
source('analysis/figures/default-figure-styling.R') # theme & color palette

theme_set(theme_get() + theme(legend.position = 'none',
                              text = element_text(size = 22)))

hr_lab <- 'Range size' # label for y axis

# hypothesis for change in H over E(R) ----
d <-
  tibble(mu = seq(0, 1, length.out = 400),
         sigma2 = rev(mu),
         split = mu < 0.35,
         h_1 = if_else(mu < 0.2, 2, sinpi(0.5 - (mu - 0.2) / 2)^10*2) +0.2,
         h_2 = if_else(split, 2 * mu^2 - 7 * mu + 3.916873, h_1),
         h_3 = if_else(split, h_2 + 30 * (0.35 - mu)^2, h_2))

notes <- tibble(x = c(0.13, 0.32),
                y = c(2.05, 3.2),
                text = c('range-resident', 'nomadic or\ndispersing'))

p_mu <-
  ggplot() +
  coord_cartesian(ylim = c(0, 4)) +
  geom_ribbon(aes(mu, ymin = h_1, ymax = h_3), d, color = pal[1],
              alpha = 0.2, linewidth = 2, fill = pal[1]) +
  geom_text(aes(x = x, y = y, label = text), notes, size = 6,
            fontface = 'bold') +
  scale_x_continuous(bquote(paste(bold('Resource abundance, E('),
                                       bolditalic(R), bold(')'))),
                     breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(hr_lab, breaks = 0, expand = c(0, 0)); p_mu

# hypothesis for change in H over V(R) ----
p_s2 <-
  ggplot() +
  coord_cartesian(ylim = c(0, 4)) +
  geom_ribbon(aes(mu, ymin = h_1, ymax = h_3), d, color = pal[2],
              alpha = 0.2, linewidth = 2, fill = pal[2]) +
  geom_text(aes(x = x, y = y, label = text), notes, size = 6,
            fontface = 'bold') +
  scale_x_reverse(bquote(paste(bold(Resource~'stochasticity, Var('),
                                    bolditalic(R),bold(')'))),
                  breaks = NULL,
                  expand = c(0, 0)) +
  scale_y_continuous(hr_lab, breaks = 0, expand = c(0, 0)); p_s2

# save the plots ----
# individual plots
ggsave('figures/mean-abundance-hr-hypotheses.png', plot = p_mu,
       width = HALF_WIDTH, height = HALF_WIDTH, scale = 2.5, units = 'mm',
       dpi = 600, bg = 'white')
ggsave('figures/variance-abundance-hr-hypotheses.png', plot = p_s2,
       width = HALF_WIDTH, height = HALF_WIDTH, scale = 2.5, units = 'mm',
       dpi = 600, bg = 'white')

# figure the GRS and GRC poster
p_grid <- plot_grid(p_mu, NULL, p_s2, labels = c('a.', '', 'b.'),
                    rel_widths = c(1, 0.5, 1), label_size = 22, nrow = 1)
ggsave('figures/2023-GRS-GRC-barga-italy/grid-hr-hypotheses.png',
       plot = p_grid, width = 60.25, height = 20, units = 'cm', dpi = 600,
       bg = 'white')

# single figure for manuscript
p_grid <- plot_grid(p_mu, p_s2, labels = c('A', 'B'),
                    label_size = 22, nrow = 1)
ggsave('figures/hr-hypotheses.png', plot = p_grid, scale = 2.5,
       width = FULL_WIDTH, height = HALF_WIDTH, units = 'mm', dpi = 600,
       bg = 'white')
