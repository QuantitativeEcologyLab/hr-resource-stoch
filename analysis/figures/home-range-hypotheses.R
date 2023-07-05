library('ctmm')    # for movement modeling
library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('ggplot2') # for fancy plots
library('ggbrace') # for curly braces on fancy plots
library('cowplot') # for fancy multi-panel plots
source('analysis/figures/default-figure-styling.R') # for common theme and color palette

theme_set(theme_get() + theme(legend.position = 'none',
                              text = element_text(size = 22)))

hr_lab <- 'Spatial needs' # label for y axis

# hypothesis for change in H over E(R) ----
d <-
  tibble(mu = seq(0, 1, length.out = 400),
         sigma2 = rev(mu),
         split = mu < 0.35,
         h_1 = if_else(mu < 0.2, 2, sinpi(0.5 - (mu - 0.2) / 2)^10 * 2) + 0.2,
         h_2 = if_else(split, 2 * mu^2 - 7 * mu + 3.916873, h_1),
         h_3 = if_else(split, h_2 + 30 * (0.35 - mu)^2, h_2))

notes <- tibble(x = c(0.13, 0.3, 0.175, 0.68),
                y = c(2.05, 3.2, 1.35, 2.05),
                text = c('range-resident', 'nomadic or\ndispersing',
                         'evolutionary timescale', 'ecological timescale'))

p_mu <-
  ggplot() +
  coord_cartesian(ylim = c(0, 4)) +
  geom_ribbon(aes(mu, ymin = h_1, ymax = h_3), d, color = pal[1], alpha = 0.2,
              linewidth = 2, fill = pal[1]) +
  geom_brace(aes(x = c(0, 0.345), y = c(1.45, 1.7)), rotate = 180) +
  geom_brace(aes(x = c(0.36, 1), y = c(1.95, 1.7)), rotate = 0) +
  geom_text(aes(x = x, y = y, label = text), notes, size = 6) +
  scale_x_continuous('Resource abundance, \U1D53C(\U1D445)', breaks = NULL,
                     expand = c(0, 0)) +
  scale_y_continuous(hr_lab, breaks = 0, expand = c(0, 0)); p_mu

# hypothesis for change in H over V(R) ----
p_s2 <-
  ggplot() +
  coord_cartesian(ylim = c(0, 4)) +
  geom_ribbon(aes(mu, ymin = h_1, ymax = h_3), d, color = pal[2], alpha = 0.2,
              linewidth = 2, fill = pal[2]) +
  geom_brace(aes(x = c(0, 0.345), y = c(1.45, 1.7)), rotate = 180) +
  geom_brace(aes(x = c(0.36, 1), y = c(1.95, 1.7)), rotate = 0) +
  geom_text(aes(x = x, y = y, label = text), notes, size = 6) +
  scale_x_reverse('Resource stochasticity, \U1D54D(\U1D445)', breaks = NULL,
                  expand = c(0, 0)) +
  scale_y_continuous(hr_lab, breaks = 0, expand = c(0, 0)); p_s2

# save the plots ----
ggsave('figures/mean-abundance-hr-hypotheses.png', plot = p_mu, width = 8,
       height = 8, scale = 1, dpi = 300, bg = 'white')
ggsave('figures/variance-abundance-hr-hypotheses.png', plot = p_s2, width = 8,
       height = 8, scale = 1, dpi = 300, bg = 'white')

p_grid <- plot_grid(p_mu, NULL, p_s2, labels = c('a.', '', 'b.'),
                    rel_widths = c(1, 0.5, 1), label_size = 22, nrow = 1)
ggsave('figures/2023-GRS-GRC-barga-italy/grid-hr-hypotheses.png', plot = p_grid,
       width = 60.25, height = 20, units = 'cm', dpi = 600, bg = 'white')

# add simulated animal movement ----
p_track <-
  simulate(ctmm(tau = c(Inf, 1), sigma = 1), t = seq(1, 400, by = 0.05), seed = 5) %>%
  SpatialPoints.telemetry() %>%
  data.frame() %>%
  ggplot(aes(x, y)) +
  coord_equal() +
  geom_path() +
  scale_x_continuous(NULL, breaks = NULL) +
  scale_y_continuous(NULL, breaks = NULL) +
  theme(plot.background = element_blank())

p_e <-
  ggdraw(p_mu) +
  draw_plot(p_track, x = 1, y = 1.0335, width = 0.45, height = 0.45, hjust = 1, vjust = 1)

p_v <-
  ggdraw(p_s2) +
  draw_plot(p_track, x = 0.096, y = 1.0335, width = 0.45, height = 0.45, hjust = 0,
            vjust = 1)

# condense the two figures into a single plot
notes_both <- tibble(x = c(0.13, 0.3, 0.82, 0.52),
                     y = c(2.05, 3.5, 0.88, 2.90),
                     text = c('range-resident',
                              'nomadic or\ndispersing',
                              'ecological\ntimescale',
                              'evolutionary\ntimescale'),
                     timescale = grepl('timescale', text))
# segments <- tibble(x1 = c(0.50, 0.50),
#                    x2 = c(0.25, 0.75),
#                    y1 = c(3.0, 3.0),
#                    y2 = c(2.7, 2.7))
# 
p_both <-
  ggplot() +
  coord_cartesian(ylim = c(0, 4), expand = FALSE) +
  # geom_ribbon(aes(mu, ymin = h_1, ymax = h_3), d, color = pal[1], fill = pal[1],
  #             alpha = 0.2, linewidth = 1) +
  # geom_ribbon(aes(sigma2, ymin = h_1, ymax = h_3), d, color = pal[2],
  #             fill = pal[2], alpha = 0.2, linewidth = 1) +
  # geom_brace(aes(x = 0.67 + c(0, 0.05), y = c(0, 1.75)), rotate = 90) +
  # geom_brace(aes(x = 0.6 + c(0, 0.05), y = c(1.8, 4)), rotate = 270) +
  geom_ribbon(aes(value, ymin = h_1, ymax = h_3, color = name, fill = name),
              pivot_longer(d, cols = c(mu, sigma2)),
              alpha = 0.2, linewidth = 1) +
  geom_text(aes(x = x, y = y, label = text), filter(notes_both, ! timescale),
            fontface = 'bold') +
  scale_x_continuous(NULL, breaks = NULL) +
  scale_y_continuous('Positional variance', breaks = 0) +
  scale_color_manual(NULL, values = pal, aesthetics = c('color', 'fill')) +
  theme(axis.title.y = element_text(color = pal[3]),
        legend.position = 'bottom', text = element_text(face = 'bold')); p_both

# p_both <-
#   ggplot() +
#   coord_cartesian(ylim = c(0, 4), expand = FALSE) +
#   geom_ribbon(aes(value, ymin = h_1, ymax = h_3, color = name, fill = name),
#               pivot_longer(d, cols = c(mu, sigma2)),
#               alpha = 0.2, linewidth = 1) +
#   geom_text(aes(x = 0.5, y = 1.6, label = 'range-resident'), data = NULL) +
#   geom_text(aes(x = 0.5, y = 3.5, label = 'nomadic or\ndispersing'), data = NULL) +
#   geom_segment(aes(x = 0.5, xend = 0.5 + c(-0.35, 0.35), y = 1.7, yend = 2.17),
#                arrow = grid::arrow()) +
#   geom_segment(aes(x = 0.5, xend = 0.5 + c(-0.275, 0.275), y = 3.35, yend = 3),
#                arrow = grid::arrow()) +
#   scale_x_continuous(NULL, breaks = NULL) +
#   scale_y_continuous('Positional variance', breaks = 0) +
#   scale_color_manual(NULL, values = pal, aesthetics = c('color', 'fill')) +
#   theme(axis.title.y = element_text(color = pal[3]), legend.position = 'bottom')

p_braces <-
  ggplot() +
  geom_brace(aes(x = c(0, 0.02), y = c(0, 1.75)), rotate = 90) +
  geom_brace(aes(x = c(0, 0.02), y = c(1.8, 4)), rotate = 90) +
  geom_text(aes(x = 0.04, y = y, label = text), filter(notes_both, timescale),
            fontface= 'bold') +
  scale_y_continuous(limits = c(-0.4, 4.05), expand = c(0, 0)) +
  xlim(x = c(0, 0.075)) +
  theme_void(); p_braces

p_final <-
  (p_both +
     theme(legend.position = 'none',
           axis.title.y = element_text(color = pal[3], size = 14))) %>%
  add_sub(label = 'Resource abundance', color = pal[1], fontface = 'bold',
          vpadding = grid::unit(1, 'line')) %>%
  add_sub(label = 'Resource stochasticity', color = pal[2],
          fontface = 'bold', vpadding = grid::unit(0, 'line')) %>%
  ggdraw() %>%
  plot_grid(p_braces, nrow = 1, rel_widths = c(4, 1)); p_final

ggsave('figures/mean-variance-hr-hypotheses.png', plot = p_final,
       width = 6, height = 5, scale = 1.5, dpi = 600, bg = 'white')

XLAB <- ggplot() +
  theme_void() +
  geom_text(aes(x = 0.25, y = 0), label = 'Resource abundance', color = pal[1]) +
  geom_text(aes(x = 0.5, y = 0), label = '/', color = 'black') +
  geom_text(aes(x = 0.75, y = 0), label = 'Resource stochasticity', color = pal[2]) +
  xlim(c(0, 1))

add_sub(p_both, 'Resource abundance', color = pal[1], x = 0.5) %>%
  add_sub('Resource stochasticity', color = pal[2], x = 0.5) %>%
  ggdraw()
