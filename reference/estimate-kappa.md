# Concentration parameter of von Mises distribution

Estimates the concentration parameter of a von Mises distribution, given
a set of angular measurements.

## Usage

``` r
est.kappa.MLE(x, w = NULL, bias = FALSE)

est.kappa(x, w = NULL, p = 2)
```

## Arguments

- x:

  numeric. Angles in degrees

- w:

  numeric. Weightings

- bias:

  logical parameter determining whether a bias correction is used in the
  computation of the MLE. Default for bias is `FALSE` for no bias
  correction.

- p:

  integer. Number of parameters in the data space: 2 for circle (the
  default), 3 for a sphere.

## Value

numeric. Concentration of a von Mises distribution

## Details

`est.kappa.MLE()` is the maximum likelihood estimate for MLE for
\\\kappa\\.

`est.kappa()` uses an approximation based on the empirical equation:
\$\$\kappa = \frac{\bar{R}(p-\bar{R}^2)}{1-\bar{R}^2} \$\$ where
\\\bar{R}\\ is the mean resultant length and \\p\\ is the dimensionality
of the data (2 for circular data).

## See also

[`vonmises()`](https://tobiste.github.io/tectonicr/reference/vonmises.md)

## Examples

``` r
set.seed(123)
x <- rvm(100, 90, 10)
w <- weighting(runif(100, 0, 10))

est.kappa(x, w)
#> [1] 10.69211

est.kappa.MLE(x, w)
#> [1] 10.33191
```
