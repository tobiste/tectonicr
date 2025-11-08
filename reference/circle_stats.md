# Summary Statistics of Circular Data

Calculate the (weighted median) and standard deviation of orientation
data.

## Usage

``` r
circular_mean(x, w = NULL, axial = TRUE, na.rm = TRUE)

circular_var(x, w = NULL, axial = TRUE, na.rm = TRUE)

circular_sd(x, w = NULL, axial = TRUE, na.rm = TRUE)

circular_median(x, w = NULL, axial = TRUE, na.rm = TRUE)

circular_quantiles(x, w = NULL, axial = TRUE, na.rm = TRUE)

circular_IQR(x, w = NULL, axial = TRUE, na.rm = TRUE)
```

## Arguments

- x:

  numeric vector. Values in degrees.

- w:

  (optional) Weights. A vector of positive numbers and of the same
  length as `x`.

- axial:

  logical. Whether the data are axial, i.e. pi-periodical (`TRUE`, the
  default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`).

- na.rm:

  logical value indicating whether `NA` values in `x` should be stripped
  before the computation proceeds.

## Value

numeric vector

## Note

Weighting may be the reciprocal of the data uncertainties.

Weightings have no effect on quasi-median and quasi-quantiles if
`length(x) %% 2 != 1` and `length(x) %% 4 == 0`, respectively.

## References

Mardia, K.V. (1972). Statistics of Directional Data: Probability and
Mathematical Statistics. London: Academic Press.

Mardia, K.V., and Jupp, P.E (1999). Directional Statistics, Wiley Series
in Probability and Statistics. John Wiley & Sons, Inc., Hoboken, NJ,
USA. [doi:10.1002/9780470316979](https://doi.org/10.1002/9780470316979)

N.I. Fisher (1993) Statistical Analysis of Circular Data, Cambridge
University Press.

Ziegler, M. O.; Heidbach O. (2019). Manual of the Matlab Script
Stress2Grid v1.1. *WSM Technical Report* 19-02, GFZ German Research
Centre for Geosciences.
[doi:10.2312/wsm.2019.002](https://doi.org/10.2312/wsm.2019.002)

Heidbach, O., Tingay, M., Barth, A., Reinecker, J., Kurfess, D., &
Mueller, B. (2010). Global crustal stress pattern based on the World
Stress Map database release 2008. *Tectonophysics* **482**,
3\<U+2013\>15,
[doi:10.1016/j.tecto.2009.07.023](https://doi.org/10.1016/j.tecto.2009.07.023)

## Examples

``` r
set.seed(1)
x <- rvm(10, 0, 100) %% 180
unc <- stats::runif(100, 0, 10)
w <- weighting(unc)
circular_mean(x, w)
#> [1] NA
circular_var(x, w)
#> [1] NA
circular_sd(x, w)
#> [1] NA
circular_median(x, w)
#> [1] NA
circular_quantiles(x, w)
#> 25% 50% 75% 
#>  NA  NA  NA 
circular_IQR(x, w)
#> [1] NA

data("san_andreas")
w2 <- weighting(san_andreas$unc)
circular_mean(san_andreas$azi)
#> [1] 10.64134
circular_mean(san_andreas$azi, w2)
#> [1] 10.85382
circular_median(san_andreas$azi)
#> [1] 35.5
circular_median(san_andreas$azi, w2)
#> [1] 35.53572
circular_quantiles(san_andreas$azi)
#>   25%   50%   75% 
#>  15.0  35.5 160.0 
circular_quantiles(san_andreas$azi, w2)
#>       25%       50%       75% 
#>  15.00000  35.53572 160.00000 
circular_var(san_andreas$azi)
#> [1] 0.3177094
circular_var(san_andreas$azi, w2)
#> [1] 0.2927477

data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
circular_mean(sa.por$azi.PoR, w2)
#> [1] 140.8843
circular_median(sa.por$azi.PoR, w2)
#> [1] 136.8329
circular_var(sa.por$azi.PoR, w2)
#> [1] 0.2614357
circular_quantiles(sa.por$azi.PoR, w2)
#>      25%      50%      75% 
#> 124.8483 136.8329 150.2174 
```
