#' function to get home range estimate starting from positional variance (`sigma`)
#' **NOTE:** only works when AKDE is isotropic and Gaussian!
#' @param .sigma positional variance estimated by a CTMM model
#' @param quantile quantile within the interval (0, 1) for which HR should be calculated
get_hr <- function(.sigma, quantile) -2 * log(1 - quantile) * pi * .sigma
