library('ggplot2') # for fancy figures
library('stringi') # for working with strings

# set default ggplot theme
# cannot bold text because unicode characters fail
theme_set(theme_bw() +
            theme(legend.position = 'none',
                  text = element_text(size = 15)))

# custom color-blind palette
pal <- c('#ff8c00', '#4477AA', '#009900', '#66CCEE',
         '#CCBB44', '#EE6677', '#AA3377', '#BBBBBB')

# custom NDVI color palette
create_ndvi_pal <- colorRampPalette(c('darkblue', 'dodgerblue', '#744700', '#d9bb94',
                                      'darkgreen'))
ndvi_pal <- create_ndvi_pal(100)

# mean and variance notations
e_r <- expression(E(italic(R)))
var_r <- expression(Var(italic(R)))
h <- 'Space-use requirements'
