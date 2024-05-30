library('ggplot2') # for fancy figures
library('khroma')  # for color palettes

.check_colors <- FALSE

# set default ggplot theme
theme_set(theme_bw() +
            theme(legend.position = 'none',
                  text = element_text(size = 15, face = 'bold'),
                  panel.grid = element_blank()))

# custom color-blind-friendly palette
pal <- c('#F09D00', '#3388E1', '#1C7734', '#992266')

if(.check_colors) plot_scheme_colorblind(pal)

# custom NDVI color palette
create_ndvi_pal <- colorRampPalette(c('darkblue', 'dodgerblue', '#744700',
                                      '#d9bb94', 'darkgreen'))
ndvi_pal <- create_ndvi_pal(100)
if(.check_colors) plot_scheme_colorblind(ndvi_pal)

# mean and variance notations
e_r <- bquote(paste(bold('E('), bolditalic('R'), bold(')')))
var_r <- bquote(paste(bold('Var('), bolditalic('R'), bold(')')))
h <- bquote(paste(bold('Home-range size, '), bolditalic('H')))
