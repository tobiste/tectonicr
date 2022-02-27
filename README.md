<!-- badges: start -->
[![R-CMD-check](https://github.com/tobiste/PlateTectonicStressR/workflows/R-CMD-check/badge.svg)](https://github.com/tobiste/PlateTectonicStressR/actions)
[![Codecov test coverage](https://codecov.io/gh/tobiste/PlateTectonicStressR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/tobiste/PlateTectonicStressR?branch=main)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5885309.svg)](https://doi.org/10.5281/zenodo.5885309)
<!-- badges: end -->


# PlateTectonicStressR

**PlateTectonicStressR** is a free and open-source **R** package for modeling and analyzing the direction of the maximum horizontal stress based on the empirical link between the orientation of intraplate stress and the direction of the relative motion of neighboring plates.

## Prerequisites

You must have R installed on your system (see http://r-project.org). Additionally, to install PlateTectonicStressR from Github, you also need the `remotes` package. This can be installed by typing the following code at the R command line prompt:

```
install.packages("remotes")
```

## Installation

The most recent development version of PlateTectonicStressR is available from Github and can be installed on your system as follows:

```
remotes::install_github('tobiste/PlateTectonicStressR')
library('PlateTectonicStressR')
```

## Documentation
https://tobiste.github.io/PlateTectonicStressR/articles/PlateTectonicStressR.html

## Author
Tobias Stephan

## How to cite
When referencing this package, please cite the package DOI ([10.5281/zenodo.5834182](https://doi.org/10.5281/zenodo.5834182)).


## Useful References
- <div class="csl-entry">Wdowinski, S. (1998). A theory of intraplate tectonics. <i>Journal of Geophysical Research: Solid Earth</i>, <i>103</i>(3), 5037–5059. http://dx.doi.org/10.1029/97JB03390</div>

- <div class="csl-entry">Heidbach, O., Reinecker, J., Tingay, M., Müller, B., Sperner, B., Fuchs, K., &#38; Wenzel, F. (2007). Plate boundary forces are not enough: Second- and third-order stress patterns highlighted in the World Stress Map database. <i>Tectonics</i>, <i>26</i>(6), n/a-n/a. https://doi.org/10.1029/2007TC002133</div>

- <div class="csl-entry">Heidbach, O., Rajabi, M., Reiter, K., Ziegler, M., &#38; Team, W. (2016). <i>World Stress Map Database Release 2016. V. 1.1</i>. GFZ Data Services. https://doi.org/10.5880/WSM.2016.001</div>

- <div class="csl-entry">Zoback, M. Lou, Zoback, M. D., Adams, J., Assumpção, M., Bell, S., Bergman, E. A., Blümling, P., Brereton, N. R., Denham, D., Ding, J., Fuchs, K., Gay, N., Gregersen, S., Gupta, H. K., Gvishiani, A., Jacob, K., Klein, R., Knoll, P., Magee, M., … Zhizhin, M. (1989). Global patterns of tectonic stress. <i>Nature</i>, <i>341</i>(6240), 291–298. https://doi.org/10.1038/341291a0</div>

## License
GPL-3.0 License
