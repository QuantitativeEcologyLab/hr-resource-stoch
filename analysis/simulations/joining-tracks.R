library('ctmm')    # for movement modeling
library('dplyr')   # for data wrangling (%>%, mutate(), slice(), etc.)
library('ggplot2') # for fancy plots
theme_set(theme_bw()) # change default ggplot theme

# OUF model ----
m <- ctmm(tau = c(10, 2), # tau_p can be Inf
          sigma = 5,
          mu = c(0,0)) # mean of the stationary process, does not work with OUF?

set.seed(1) # for consistent results
sim1 <- simulate(m,t= seq(1, 10, by = 0.005))
d <- sim1[c(1, nrow(sim1)), ]
sim2 <- simulate(m, data = d, t = seq(1, 10, by = 0.005))

# does not start and end at (0, 0)
ggplot() +
  coord_equal() +
  geom_path(aes(x, y), sim1, color = 'black') +
  geom_path(aes(x, y), sim2, color = 'red') +
  geom_point(aes(0, 0), col = 'blue')

# IOU model ----
DELTA_T <- 0.1 # time between measurements
SAMPLES <- seq(0, 15, by = DELTA_T) # sampling times
m <- ctmm(tau = c(Inf, 1), sigma = 2, mu = c(0, 0)) # infinitely diffusive model

sim1 <- simulate(m, t = SAMPLES)
d <- sim1[c(1, nrow(sim1)), ]
sim2 <- simulate(m, data = d, t = SAMPLES)

# starts and ends at (0, 0)
ggplot() +
  coord_equal() +
  geom_path(aes(x, y), sim1, color = 'black') +
  geom_path(aes(x, y), sim2, color = 'red') +
  geom_point(aes(0, 0), col = 'blue') +
  geom_point(aes(x, y), filter(data.frame(sim1), t == max(t)), col = 'darkorange')

# both tracks actually start from (0, 0), so neither is returning to (0, 0)
ggplot() +
  coord_equal() +
  geom_path(aes(x, y, color = t), sim1) +
  geom_path(aes(x, y, color = t), sim2) +
  geom_point(aes(0, 0), col = 'blue') +
  geom_point(aes(x, y), filter(data.frame(sim1), t == max(t)), col = 'darkorange')

# if d is a tibble, tracks won't match
d <- sim1[c(1, nrow(sim1)), ] %>% data.frame() %>% tibble()
sim2 <- simulate(m, data = d, t = SAMPLES)

# starts and ends at (0, 0)
ggplot() +
  coord_equal() +
  geom_path(aes(x, y), sim1, color = 'black') +
  geom_path(aes(x, y), sim2, color = 'red') +
  geom_point(aes(0, 0), col = 'blue') +
  geom_point(aes(x, y), filter(data.frame(sim1), t == max(t)), col = 'darkorange')

# can join at multiple points
d <- sim1[c(1, 75, nrow(sim1)), ]
sim2 <- simulate(m, data = d, t = SAMPLES)

ggplot() +
  coord_equal() +
  geom_path(aes(x, y), sim1, color = 'black') +
  geom_path(aes(x, y), sim2, color = 'red') +
  geom_point(aes(0, 0), col = 'blue') +
  geom_point(aes(x, y), sim1[75, ], col = 'darkgreen') +
  geom_point(aes(x, y), filter(data.frame(sim1), t == max(t)), col = 'darkorange')

# of course, the more points you add, the more similar the tracks will be
d <- sim1[1:15 * 10, ]
sim2 <- simulate(m, data = d, t = SAMPLES)

ggplot() +
  coord_equal() +
  geom_path(aes(x, y), sim1, color = 'black') +
  geom_path(aes(x, y), sim2, color = 'red') +
  geom_point(aes(x, y), sim1[(1:15) * 10, ])
