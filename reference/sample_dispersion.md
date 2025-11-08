# Sample circular dispersion

Alternative versions of variance, dispersion a distance (Mardia and
Jupp, 1999; pp. 19-20). These alternative dispersion has a minimum at
the sample median.

## Usage

``` r
sample_circular_variance(x, w = NULL, axial = TRUE)

sample_circular_distance(x, y, axial = TRUE, na.rm = TRUE)

sample_circular_dispersion(
  x,
  y = NULL,
  w = NULL,
  w.y = NULL,
  axial = TRUE,
  na.rm = TRUE
)
```

## Arguments

- x, y:

  vectors of numeric values in degrees. `length(y)` is either `1` or
  `length(x)`

- w, w.y:

  (optional) Weights. A vector of positive numbers and of the same
  length as `x`. `w.y` is the (optional) weight of `y`.

- axial:

  logical. Whether the data are axial, i.e. pi-periodical (`TRUE`, the
  default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`).

- na.rm:

  logical. Whether `NA` values in `x` should be stripped before the
  computation proceeds.

## References

N.I. Fisher (1993) Statistical Analysis of Circular Data, Cambridge
University Press.

Mardia, K.V., and Jupp, P.E (1999). Directional Statistics, Wiley Series
in Probability and Statistics. John Wiley & Sons, Inc., Hoboken, NJ,
USA. [doi:10.1002/9780470316979](https://doi.org/10.1002/9780470316979)

## Examples

``` r
a <- c(0, 2, 359, 6, 354)
sample_circular_distance(a, 10) # distance to single value
#> [1] 5.0 4.0 5.5 2.0 8.0

b <- a + 90
sample_circular_distance(a, b) # distance to multiple values
#> [1] 45 45 45 45 45

data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
sample_circular_variance(sa.por$azi.PoR)
#> [1] 0.6037799
sample_circular_dispersion(sa.por$azi.PoR, y = 135)
#> [1] 10.88795
sample_circular_dispersion(sa.por$azi.PoR, y = 135, w = weighting(san_andreas$unc))
#> [1] 10.41816
```
