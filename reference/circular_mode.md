# Circular Mode

MLE angle (maximum density) using a von Mises distribution kernel with
specified concentration.

## Usage

``` r
circular_mode(x, kappa = NULL, axial = TRUE, n = 512)
```

## Arguments

- x:

  numeric vector. Values in degrees.

- kappa:

  von Mises distribution concentration parameter. Will be estimated
  using
  [`est.kappa()`](https://tobiste.github.io/tectonicr/reference/est.kappa.md)
  if not provided.

- axial:

  logical. Whether the data are axial, i.e. pi-periodical (`TRUE`, the
  default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`).#' @param
  kappa

- n:

  the number of equally spaced points at which the density is to be
  estimated.

## Value

numeric

## References

N.I. Fisher (1993) Statistical Analysis of Circular Data, Cambridge
University Press.

## Examples

``` r
set.seed(1)
x <- rvm(10, 0, 100)
circular_mode(x, kappa = est.kappa(x))
#> [1] 358.591
```
