library('ctmm')  # for continuous-time movement modeling
library('dplyr') # for data wrangling

pull_parameter <- function(.object, .parameter, .quantile,
                           .value = c('low', 'est', 'high')) {
  CIs <- summary(.object, units = FALSE, level.UD = .quantile)$CI
  
  if(length(.value) > 1) stop('Please choose a single value from {low, est, high}')
  
  # take the model summary and extract the confidence intervals
  CIs %>%
    as_tibble() %>% # convert to a tibble (i.e., fancy data.frame)
    filter(grepl(pattern = .parameter, rownames(CIs))) %>%
    pull(get(.value))
}
