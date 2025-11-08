# Sample Circular Median and Deviation

Sample median direction for a vector of circular data

## Usage

``` r
circular_sample_median(x, axial = TRUE, na.rm = TRUE)

circular_sample_median_deviation(x, axial = TRUE, na.rm = TRUE)
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

numeric

## References

N.I. Fisher (1993) Statistical Analysis of Circular Data, Cambridge
University Press.

## Examples

``` r
set.seed(1)
x <- rvm(n = 100, mean = 0, kappa = 10)
circular_sample_median(x)
#> [1] 175.2845
circular_sample_median_deviation(x)
#> [1] 12.74011

data("san_andreas")
circular_sample_median(san_andreas$azi)
#> [1] 9
circular_sample_median_deviation(san_andreas$azi)
#> [1] 19.77709
```
