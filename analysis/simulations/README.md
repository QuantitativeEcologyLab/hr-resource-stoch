# Main scripts (to be run in the following order)

1. `analysis/simulations/1-tracks.R` (see checks below before running 2)
2. `analysis/simulations/2-hr-mean-variance-simulations-days.R`
3. `analysis/simulations/3-hr-mean-variance-simulations-days-summarized.R`
4. `analysis/simulations/4-hr-mean-variance-simulations-modeling.R`
5. `analysis/simulations/5-hr-mean-variance-simulations-hrs.R`
6. `analysis/simulations/6-modeling-R-and-hr.R`

# Checks after 1:

1. `analysis/simulations/sensitivity-analyses/1a-return-sensitivity.R`
2. `analysis/simulations/sensitivity-analyses/1b-delta-t-sensitivity.R`
3. `analysis/simulations/sensitivity-analyses/1c-hr-simulation-extreme-scenarios.R`

The `analysis/simulations/joining-tracks.R` script includes information on joining tracks at given points in space and time.

The `analysis/simulations/movement-model` script includes the movement model to generate the tracks and the raster used to determine when the organism gathers resources.
