<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/tobiste/tectonicr/workflows/R-CMD-check/badge.svg)](https://github.com/tobiste/tectonicr/actions)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7510800.svg)](https://doi.org/10.5281/zenodo.7510800)
<!-- badges: end -->

# tectonicr

**tectonicr** is a free and open-source **R** package for modeling and analyzing the direction of the maximum horizontal stress (SHmax) based on the empirical link between the direction of intraplate stress and the direction of the relative motion of neighboring plates. The following methods are available:

- **Theoretical direction of SHmax**: The predicted stress field adjacent to a plate boundary is calculated using the relative plate motion of the neighboring plates using the function `model_shmax()`. The goodness-of-fit can be statistically tested by `norm_chisq()`, `norm_rayleigh()`, or `confidence_interval()`.
- **Distance to plate boundary**: `distance_from_pb()` gives the distance between the stress data point and the plate boundary measured along the stress trajectories.
- **Visualization of the trajectories of the theoretical stress field** in terms of small circles, great circles, and lines of constant bearing. The `eulerpole_paths()` functions generates an  `sf` object containing spatial information that is suitable to plot with, for instance, `ggplot()`. 
- **Relative rotations from a given set of plate motion parameters**: `equivalent_rotation()` transfers a set of plate motion parameters into the relative plate motions among the given plates. 
- **Average direction of a set of SHmax data** using the (weighted) mean, quasi-median, and other parameters to statistically estimate pi-directional data. 
- **Spatial interpolation of of SHmax**: `PoR_stress2grid()` uses distance, method, and quality-weighted mean direction of stress data without being affected by angular distortions.

## Prerequisites

You must have R installed on your system (see http://r-project.org). To install **tectonicr** from Github, you also need the `remotes` package. This can be installed by typing the following code at the R command line prompt:

```
install.packages("remotes")
```

## Installation

The most recent development version of **tectonicr** is available from Github and can be installed on your system as follows:

```
remotes::install_github('tobiste/tectonicr')
library('tectonicr')
```

## Documentation
https://tobiste.github.io/tectonicr/articles/tectonicr.html

## Author
Tobias Stephan

## How to cite
When referencing this package, please cite the package DOI: [10.5281/zenodo.7510800](https://doi.org/10.5281/zenodo.7510800).


## Useful References
- <div class="csl-entry">Wdowinski, S. (1998). A theory of intraplate tectonics. <i>Journal of Geophysical Research: Solid Earth</i>, <i>103</i>(3), 5037–5059. http://dx.doi.org/10.1029/97JB03390</div>

- <div class="csl-entry">Heidbach, O., Rajabi, M., Reiter, K., Ziegler, M., &#38; Team, W. (2016). <i>World Stress Map Database Release 2016. V. 1.1</i>. GFZ Data Services. https://doi.org/10.5880/WSM.2016.001</div>

- <div class="csl-entry">Mardia, K. V., and Jupp, P. E. (Eds.). (1999). <i>Directional Statistics</i>. Hoboken, NJ, USA: John Wiley & Sons, Inc. https://doi.org/10.1002/9780470316979</div>

## License
GPL-3.0 License
