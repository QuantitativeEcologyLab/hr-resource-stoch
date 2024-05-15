library('ctmm')    # for movement modeling
library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('ggplot2') # for fancy plots
library('ggbrace') # for curly braces on fancy plots
library('cowplot') # for fancy multi-panel plots
source('analysis/figures/default-figure-styling.R') # theme & color palette

theme_set(theme_get() + theme(legend.position = 'none',
                              text = element_text(size = 22)))

hr_lab <- 'Space-use requirements' # label for y axis

# hypothesis for change in H over E(R) ----
d <-
  tibble(mu = seq(0, 1, length.out = 400),
         sigma2 = rev(mu),
         split = mu < 0.35,
         h_1 = if_else(mu < 0.2, 2, sinpi(0.5 - (mu - 0.2) / 2)^10*2) +0.2,
         h_2 = if_else(split, 2 * mu^2 - 7 * mu + 3.916873, h_1),
         h_3 = if_else(split, h_2 + 30 * (0.35 - mu)^2, h_2))

notes <- tibble(x = c(0.13, 0.32, 0.19, 0.68),
                y = c(2.05, 3.2, 1.35, 2.05),
                text = c('range-resident', 'nomadic or\ndispersing',
                         'evolutionary timescale', 'ecological timescale'))

p_mu <-
  ggplot() +
  coord_cartesian(ylim = c(0, 4)) +
  geom_ribbon(aes(mu, ymin = h_1, ymax = h_3), d, color = pal[1],
              alpha = 0.2, linewidth = 2, fill = pal[1]) +
  stat_brace(aes(x = c(0, 0.345), y = c(1.45, 1.7)), rotate = 180,
             outside = FALSE) +
  stat_brace(aes(x = c(0.36, 1), y = c(1.95, 1.7)), rotate = 0,
             outside = FALSE) +
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
  stat_brace(aes(x = c(0, 0.345), y = c(1.45, 1.7)), rotate = 180,
             outside = FALSE) +
  stat_brace(aes(x = c(0.36, 1), y = c(1.95, 1.7)), rotate = 0,
             outside = FALSE) +
  geom_text(aes(x = x, y = y, label = text), notes, size = 6,
            fontface = 'bold') +
  scale_x_reverse(bquote(paste(bold(Resource~'stochasticity, Var('),
                                    bolditalic(R),bold(')'))),
                  breaks = NULL,
                  expand = c(0, 0)) +
  scale_y_continuous(hr_lab, breaks = 0, expand = c(0, 0)); p_s2

# save the plots ----
ggsave('figures/mean-abundance-hr-hypotheses.png', plot = p_mu,
       width = 8, height = 8, scale = 1, dpi = 300, bg = 'white')
ggsave('figures/variance-abundance-hr-hypotheses.png', plot = p_s2,
       width = 8, height = 8, scale = 1, dpi = 300, bg = 'white')

# figure the GRS and GRC poster
p_grid <- plot_grid(p_mu, NULL, p_s2, labels = c('a.', '', 'b.'),
                    rel_widths = c(1, 0.5, 1), label_size = 22, nrow = 1)
ggsave('figures/2023-GRS-GRC-barga-italy/grid-hr-hypotheses.png',
       plot = p_grid, width = 60.25, height = 20, units = 'cm', dpi = 600,
       bg = 'white')
