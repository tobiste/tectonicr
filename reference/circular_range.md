# Circular Range

Length of the smallest arc which contains all the observations. The
circular range is based on the sample circular distance.

## Usage

``` r
circular_range(x, axial = TRUE, na.rm = TRUE)
```

## Arguments

- x:

  numeric vector. Values in degrees.

- axial:

  logical. Whether the data are axial, i.e. pi-periodical (`TRUE`, the
  default) or directional, i.e. \\2 \pi\\-periodical (`FALSE`).

- na.rm:

  logical value indicating whether `NA` values in `x` should be stripped
  before the computation proceeds.

## Value

numeric. angle in degrees

## References

Mardia, K.V., and Jupp, P.E (1999). Directional Statistics, Wiley Series
in Probability and Statistics. John Wiley & Sons, Inc., Hoboken, NJ,
USA. [doi:10.1002/9780470316979](https://doi.org/10.1002/9780470316979)

## See also

[`sample_circular_distance()`](https://tobiste.github.io/tectonicr/reference/sample_dispersion.md)

## Examples

``` r
roulette <- c(43, 45, 52, 61, 75, 88, 88, 279, 357)
circular_range(roulette, axial = FALSE)
#> [1] 169

data("san_andreas")
circular_range(san_andreas$azi)
#> [1] 173
```
