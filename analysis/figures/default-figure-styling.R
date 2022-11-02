library('ggplot2') # for fancy figures
library('stringi') # for working with strings

# set default ggplot theme
theme_set(theme_bw() + theme(text = element_text(face = 'bold')))

# custom color-blind palette
pal <- c('#ff8c00', '#4477AA', '#009900', '#66CCEE',
         '#CCBB44', '#EE6677', '#AA3377', '#BBBBBB')

# custom NDVI color palette
create_ndvi_pal <- colorRampPalette(c('darkblue', 'dodgerblue', '#744700', '#d9bb94',
                                      'darkgreen'))
ndvi_pal <- create_ndvi_pal(100)
