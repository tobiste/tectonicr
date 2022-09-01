<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/tobiste/tectonicr/workflows/R-CMD-check/badge.svg)](https://github.com/tobiste/tectonicr/actions)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6360893.svg)](https://doi.org/10.5281/zenodo.6360893)
[![Travis build status](https://travis-ci.com/tobiste/tectonicr.svg?branch=main)](https://travis-ci.com/tobiste/tectonicr)
<!-- badges: end -->

# tectonicr

**tectonicr** is a free and open-source **R** package for modeling and analyzing the direction of the maximum horizontal stress (SHmax) based on the empirical link between the orientation of intraplate stress and the direction of the relative motion of neighboring plates. The following methods are available:

- **Theoretical direction of SHmax**: The predicted stress field adjacent to a plate boundary is calculated using the relative plate motion of the  neighboring plates using the function `model_shmax()`. The deviation or misfit of the prediction to the observation can be obtained from the function `misfit_shmax()` and statistically evaluated by applying `norm_chisq()`.
- **Distance to plate boundary**: `distance_from_pb()` gives the distance between the stress data point and the plate boundary measured along the stress trajectories.
- **Visualization of the trajectories of the theoretical stress field** in terms of small circles, great circles, and lines of constant bearing. The `eulerpole_paths()` functions generates an  `sf` object containing spatial information that is suitable to plot with, for instance, `ggplot()`. 
- **Relative rotations from a given set of plate motion parameters**: `equivalent_rotation()` transfers a set of plate motion parameters into the relative plate motions among the given plates. 
- **Average direction of a set of SHmax data** using the (weighted) mean or median for pi-directional data. 
- **Spatial interpolation of of SHmax**: `stress2grid()` uses distance, method, and quality-weighted mean direction of stress data

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
When referencing this package, please cite the package DOI: [10.5281/zenodo.6360893](https://doi.org/10.5281/zenodo.6360893).


## Useful References
- <div class="csl-entry">Wdowinski, S. (1998). A theory of intraplate tectonics. <i>Journal of Geophysical Research: Solid Earth</i>, <i>103</i>(3), 5037–5059. http://dx.doi.org/10.1029/97JB03390</div>

- <div class="csl-entry">Heidbach, O., Reinecker, J., Tingay, M., Müller, B., Sperner, B., Fuchs, K., &#38; Wenzel, F. (2007). Plate boundary forces are not enough: Second- and third-order stress patterns highlighted in the World Stress Map database. <i>Tectonics</i>, <i>26</i>(6), n/a-n/a. https://doi.org/10.1029/2007TC002133</div>

- <div class="csl-entry">Heidbach, O., Rajabi, M., Reiter, K., Ziegler, M., &#38; Team, W. (2016). <i>World Stress Map Database Release 2016. V. 1.1</i>. GFZ Data Services. https://doi.org/10.5880/WSM.2016.001</div>

- <div class="csl-entry">Zoback, M. Lou, Zoback, M. D., Adams, J., Assumpção, M., Bell, S., Bergman, E. A., Blümling, P., Brereton, N. R., Denham, D., Ding, J., Fuchs, K., Gay, N., Gregersen, S., Gupta, H. K., Gvishiani, A., Jacob, K., Klein, R., Knoll, P., Magee, M., … Zhizhin, M. (1989). Global patterns of tectonic stress. <i>Nature</i>, <i>341</i>(6240), 291–298. https://doi.org/10.1038/341291a0</div>

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

## License
GPL-3.0 License
