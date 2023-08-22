# Main scripts (to be run in the following order)

1. `analysis/simulations/tracks.R` (see checks below before running 2)
2. `analysis/simulations/hr-mean-variance-simulations-days.R`
3. `analysis/simulations/hr-mean-variance-simulations-days-summarized.R`
4. `analysis/simulations/hr-mean-variance-simulations-modeling.R`
5. `analysis/simulations/hr-mean-variance-simulations-hrs.R`

# Checks after 1:

1. `analysis/figures/simulations/return-sensitivity.R`
2. `analysis/figures/simulations/delta-t-sensitivity.R`
3. `analysis/simulations/hr-simulation-extreme-scenarios.R`

The `analysis/simulations/joining-tracks.R` script includes information on joining tracks at given points in space.
