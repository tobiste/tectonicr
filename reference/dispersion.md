# Circular Distance and Dispersion

Circular distance between two angles and circular dispersion of angles
about a specified angle.

## Usage

``` r
circular_distance(x, y, axial = TRUE, na.rm = TRUE)

circular_dispersion(
  x,
  y = NULL,
  w = NULL,
  w.y = NULL,
  axial = TRUE,
  na.rm = TRUE
)

circular_sd2(x, y, w = NULL, axial = TRUE, na.rm = TRUE)
```

## Arguments

- x, y:

  vectors of numeric values in degrees. `length(y)` is either `1` or
  `length(x)`

- axial:

  logical. Whether the data are axial, i.e. pi-periodical (`TRUE`, the
  default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`).

- na.rm:

  logical. Whether `NA` values in `x` should be stripped before the
  computation proceeds.

- w, w.y:

  (optional) Weights. A vector of positive numbers and of the same
  length as `x`. `w.y` is the (optional) weight of `y`.

## Value

`circular_distance` returns a numeric vector of positive numbers,
`circular_dispersion` and `circular_sd2()` return a positive number.

## Details

Circular dispersion is a measure for the spread of data like the
variance. Dispersion measures the spread about a given angles, whereas
the variance measures the spread about the mean (Mardia and Jupp, 1999).
When `y = NULL` the dispersion is identical to the variance.

Circular standard deviation in `circular_sd2()` is the transformed
dispersion instead of the variance as for
[`circular_sd()`](https://tobiste.github.io/tectonicr/reference/circle_stats.md).

## Note

If `y` is `NULL`, than the circular variance is returned.

## References

Mardia, K.V. (1972). Statistics of Directional Data: Probability and
Mathematical Statistics. London: Academic Press.

Mardia, K.V., and Jupp, P.E (1999). Directional Statistics, Wiley Series
in Probability and Statistics. John Wiley & Sons, Inc., Hoboken, NJ,
USA. [doi:10.1002/9780470316979](https://doi.org/10.1002/9780470316979)

## See also

[`circular_mean()`](https://tobiste.github.io/tectonicr/reference/circle_stats.md),
[`circular_var()`](https://tobiste.github.io/tectonicr/reference/circle_stats.md).

## Examples

``` r
a <- c(0, 2, 359, 6, 354)
circular_distance(a, 10) # distance to single value
#> [1] 0.030153690 0.019369152 0.036408073 0.004865966 0.075975952

b <- a + 90
circular_distance(a, b) # distance to multiple values
#> [1] 1 1 1 1 1

data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
circular_dispersion(sa.por$azi.PoR, y = 135)
#> [1] 0.1495228
circular_dispersion(sa.por$azi.PoR, y = 135, w = weighting(san_andreas$unc))
#> [1] 0.1384805
circular_sd2(sa.por$azi.PoR, y = 135, w = weighting(san_andreas$unc))
#> [1] NA
```
