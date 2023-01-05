qgamma2 <- function(p, mu, sigma2, lower.tail = FALSE) {
  # mean = k * theta
  # sigma^2 = k * theta^2
  qgamma(p = p,
         shape = mu^2 / sigma2, # (k * theta)^2 / (k * theta^2)
         scale = sigma2 / mu,   # (k * theta^2) / (k * theta)
         lower.tail = lower.tail)
}
