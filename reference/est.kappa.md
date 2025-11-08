# Concentration parameter of von Mises distribution

Computes the maximum likelihood estimate of \\\kappa\\, the
concentration parameter of a von Mises distribution, given a set of
angular measurements.

## Usage

``` r
est.kappa(x, w = NULL, bias = FALSE, axial = TRUE)
```

## Arguments

- x:

  numeric. angles in degrees

- w:

  numeric. weightings

- bias:

  logical parameter determining whether a bias correction is used in the
  computation of the MLE. Default for bias is `FALSE` for no bias
  correction.

- axial:

  logical. Whether the data are axial, i.e. pi-periodical (`TRUE`, the
  default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`).

## Value

numeric.

## Examples

``` r
set.seed(123)
est.kappa(rvm(100, 90, 10), w = weighting(runif(100, 0, 10)))
#> [1] 3.216802
```
