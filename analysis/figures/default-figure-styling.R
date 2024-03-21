library('ggplot2') # for fancy figures

# set default ggplot theme
theme_set(theme_bw() +
            theme(legend.position = 'none',
                  text = element_text(size = 15, face = 'bold'),
                  panel.grid = element_blank()))

# custom color-blind palette
pal <- c('#ff8c00', '#4477AA', '#009900', '#66CCEE',
         '#CCBB44', '#EE6677', '#AA3377', '#BBBBBB')

# custom NDVI color palette
create_ndvi_pal <- colorRampPalette(c('darkblue', 'dodgerblue', '#744700', '#d9bb94',
                                      'darkgreen'))
ndvi_pal <- create_ndvi_pal(100)

# mean and variance notations
e_r <- bquote(paste(bold('E('), bolditalic('R'), bold(')')))
var_r <- bquote(paste(bold('Var('), bolditalic('R'), bold(')')))
h <- bquote(paste(bold('Space-use requirements, '), bolditalic('H')))
