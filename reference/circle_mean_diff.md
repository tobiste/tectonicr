# Circular Mean Difference

The circular mean difference is based on the sample circular distance

## Usage

``` r
circular_mean_difference(x, w = NULL, axial = TRUE, na.rm = TRUE)

circular_mean_difference_alt(x, w = NULL, axial = TRUE, na.rm = TRUE)
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

numeric

## References

Mardia, K.V., and Jupp, P.E (1999). Directional Statistics, Wiley Series
in Probability and Statistics. John Wiley & Sons, Inc., Hoboken, NJ,
USA. [doi:10.1002/9780470316979](https://doi.org/10.1002/9780470316979)

## See also

[`sample_circular_distance()`](https://tobiste.github.io/tectonicr/reference/sample_dispersion.md)

## Examples

``` r
data("san_andreas")
circular_mean_difference(san_andreas$azi)
#> [1] 0.6951742
circular_mean_difference(san_andreas$azi, weighting(san_andreas$unc))
#> [1] 0.6662931

circular_mean_difference_alt(san_andreas$azi)
#> [1] 27.92053
circular_mean_difference_alt(san_andreas$azi, weighting(san_andreas$unc))
#> [1] 26.6433
```
