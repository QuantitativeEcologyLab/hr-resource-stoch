rbeta2 <- function(mu, sigma2, N = 1e4) {
  
  phi <- sigma2 / (mu * (1 - mu)) # scale parameter
  a <- mu*(1 / phi - 1) # shape_1
  b <- (1 - mu) * (1 / phi - 1) # shape 2
  
  cat(paste0('a = ', a, ', b = ', b, ', phi = ', phi, '\n'))
  rbeta(n = N, shape1 = a, shape2 = b)
}


if(FALSE) {
  eps <- .Machine$double.eps
  layout(matrix(1:4, ncol = 2, byrow = TRUE))
  for(.mu in c(0.5, 0.9)) {
    for(.sigma2 in c(1e-4, 0.05)) {
      hist(rbeta2(mu = .mu, sigma2 = .sigma2 - eps), xlab = 'R',
           main = paste0('\U03BC = ', .mu, ', \U03C3\U00B2 = ', .sigma2))
    }
  }
  layout(1)
}
