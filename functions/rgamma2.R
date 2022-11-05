# re-parameterize rgamma() as a function of mean and variance
rgamma2 <- function(mu, sigma2, N = n()) {
  # mean = k * theta
  # sigma^2 = k * theta^2
  rgamma(n = N,
         shape = mu^2 / sigma2, # (k * theta)^2 / (k * theta^2)
         scale = sigma2 / mu)   # (k * theta^2) / (k * theta)
}
