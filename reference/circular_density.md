# Circular Kernel Density Estimation

Kernel density estimates for circular data from a given kernel (von
Mises or wrapped Cauchy distribution) and bandwidth

## Usage

``` r
circular_density(
  x,
  z = NULL,
  bw = NULL,
  weights = NULL,
  na.rm = TRUE,
  from = 0,
  to = 360,
  n = 512L,
  axial = TRUE,
  kappa = NULL,
  rho = NULL,
  kernel = c("vonmises", "wrappedcauchy"),
  adjust = 1,
  subdensity = FALSE
)
```

## Arguments

- x:

  numeric. A vector of angles (in degrees) from which the estimate is to
  be computed.

- z:

  numeric. Angles where the density is estimated. If `NULL` equally
  spaced angles are used according to the parameters `from`, `to` and
  `n`.

- bw, kappa, rho:

  numeric. Smoothing bandwidth expressed as the concentration parameter
  \\\kappa\\ for the von Mises distribution or \\\rho\\ for the wrapped
  Cauchy distribution. Small and large values for the von Mises and
  wrapped Cauchy distribution, respectively, gives smooth density lines.
  If not specified, parameter will be estimated using
  [`est.kappa()`](https://tobiste.github.io/tectonicr/reference/estimate-kappa.md)
  for the von Mises distribution, or set to \\p \exp(-1)\\ for the
  wrapped Cauchy distribution (where \\p = 2\\ when `axial=TRUE` and 1
  otherwise).

- weights:

  numeric. A vector of observation weights, of the same length as `x`,
  to give individual observations weight in the density estimate. Should
  sum to 1; a warning is issued if it doesn't (unless
  `subdensity = TRUE`). Defaults to equal weight `1/length(x)` per
  observation, matching
  [`stats::density()`](https://rdrr.io/r/stats/density.html).

- na.rm:

  logical; if `TRUE`, missing values are removed from `x`. If `FALSE`
  any missing values cause an error.

- from, to:

  the left and right-most points of the grid at which the density is to
  be estimated; the defaults are `cut * bw` outside of `range(x)`.

- n:

  integer. Number of equally spaced angles at which the density is to be
  estimated.

- axial:

  Logical. Whether data are uniaxial (`axial=FALSE`) or biaxial (`TRUE`,
  the default).

- kernel:

  character. The smoothing kernel to be used; one of `"vonmises"` (the
  default) or `"wrappedcauchy"` for the von Mises or the Wrapped Cauchy
  distribution.

- adjust:

  the bandwidth used is actually `adjust*bw`. This makes it easy to
  specify values like ‘half the default’ bandwidth.

- subdensity:

  logical. If `TRUE`, suppress the "sum(weights) != 1" warning, for when
  a deliberately partial (sub-)density is wanted, e.g. one group's
  contribution to a shared total. See
  [`stats::density()`](https://rdrr.io/r/stats/density.html).

## Value

Object of class `"density"`

## See also

[`stats::density()`](https://rdrr.io/r/stats/density.html),
[`dvm()`](https://tobiste.github.io/tectonicr/reference/vonmises.md) and
[`dwcauchy()`](https://tobiste.github.io/tectonicr/reference/wcauchy.md)

## Examples

``` r
w <-  weighting(san_andreas$unc)

# von Mises kernel density
circular_density(san_andreas$azi, kappa = 100)
#> 
#> Call:
#>  circular_density(x = san_andreas$azi, kappa = 100)
#> 
#> Data: san_andreas$azi (1126 obs.);   Bandwidth 'bw' = 100
#> 
#>        x             y          
#>  Min.   :  0   Min.   :0.01295  
#>  1st Qu.: 90   1st Qu.:0.04817  
#>  Median :180   Median :0.13341  
#>  Mean   :180   Mean   :0.31963  
#>  3rd Qu.:270   3rd Qu.:0.48201  
#>  Max.   :360   Max.   :1.11629  
circular_density(san_andreas$azi, weights = w, kappa = 100)
#> 
#> Call:
#>  circular_density(x = san_andreas$azi, weights = w, kappa = 100)
#> 
#> Data: san_andreas$azi (1126 obs.);   Bandwidth 'bw' = 100
#> 
#>        x             y          
#>  Min.   :  0   Min.   :0.01116  
#>  1st Qu.: 90   1st Qu.:0.04113  
#>  Median :180   Median :0.12319  
#>  Mean   :180   Mean   :0.31971  
#>  3rd Qu.:270   3rd Qu.:0.48292  
#>  Max.   :360   Max.   :1.13495  

# wrapped Cauchy kernel density
circular_density(san_andreas$azi, rho = 0.9, kernel = "wrappedcauchy")
#> 
#> Call:
#>  circular_density(x = san_andreas$azi, rho = 0.9, kernel = "wrappedcauchy")
#> 
#> Data: san_andreas$azi (1126 obs.);   Bandwidth 'bw' = 0.9
#> 
#>        x             y          
#>  Min.   :  0   Min.   :0.04110  
#>  1st Qu.: 90   1st Qu.:0.06602  
#>  Median :180   Median :0.16949  
#>  Mean   :180   Mean   :0.31945  
#>  3rd Qu.:270   3rd Qu.:0.49674  
#>  Max.   :360   Max.   :0.95665  
circular_density(san_andreas$azi, weights = w, rho = 0.9, kernel = "wrappedcauchy")
#> 
#> Call:
#>  circular_density(x = san_andreas$azi, weights = w, rho = 0.9,     kernel = "wrappedcauchy")
#> 
#> Data: san_andreas$azi (1126 obs.);   Bandwidth 'bw' = 0.9
#> 
#>        x             y          
#>  Min.   :  0   Min.   :0.03926  
#>  1st Qu.: 90   1st Qu.:0.05818  
#>  Median :180   Median :0.15915  
#>  Mean   :180   Mean   :0.31952  
#>  3rd Qu.:270   3rd Qu.:0.50049  
#>  Max.   :360   Max.   :0.97728  
```
