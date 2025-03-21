# [How resource abundance and resource stochasticity affect organisms' range sizes](https://doi.org/10.1186/s40462-025-00546-5)

Mezzini S., Fleming C.H., Medici E.P. & Noonan M.J. (2025). How resource abundance and resource stochasticity affect organisms' range sizes. Movement Ecology 13, 20. https://doi.org/10.1186/s40462-025-00546-5

This repo contains the data and code supporting the results reported in the above manuscript, with the exception of two simulated datasets that were greater than 100 MB and the tapir data. The simulated data can be produced by running the scripts in the `analysis` folder, while the tapir data is available at (https://github.com/StefanoMezzini/tapirs)[https://github.com/StefanoMezzini/tapirs].

* `analysis` contains all scripts used for analyses, including figure creation and model fitting. Scripts contain comments that explain the code, but they assume users have basic knowledge of `R` and that they will refer to help pages for packages and functions (which can be accessed with `?`, e.g.: `?data.frame`).
* `data` contains all datasets not publicly available elsewhere.
* `figures` contains all figures used in this project. The scripts used to create the figures can be found in the `analysis/figures` folder.
* `functions` contains all custom functions sourced in analysis scripts.
* `issues` contains all files related to issues with the project and repository. Please report any reproducibility issues as an issue on the repository.
* `models` contains all model saved as `.rds` files, except for models which are too large to store on GitHub. All models can be fit using the scripts in the `analysis/models` folder.
* `posters` contains all posters presented on this work.
* `simulations` contains all data from simulations (except for two that were > 100 MB)
* `writing` contains all files used to create the final manuscript, including appendices.
