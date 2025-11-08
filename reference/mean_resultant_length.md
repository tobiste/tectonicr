# Mean Resultant Length

Measure of spread around the circle. It should be noted that: If R=0,
then the data is completely spread around the circle. If R=1, the data
is completely concentrated on one point.

## Usage

``` r
mean_resultant_length(x, w = NULL, na.rm = TRUE)
```

## Arguments

- x:

  numeric vector. Values in degrees, for which the mean, median or
  standard deviation are required.

- w:

  (optional) Weights. A vector of positive numbers, of the same length
  as `x`.

- na.rm:

  logical value indicating whether `NA` values in `x` should be stripped
  before the computation proceeds.

## Value

numeric.

## References

Mardia, K.V. (1972). Statistics of Directional Data: Probability and
Mathematical Statistics. London: Academic Press.

## Examples

``` r
# Example data from Davis (1986), pp. 316
finland_stria <- c(
  23, 27, 53, 58, 64, 83, 85, 88, 93, 99, 100, 105, 113,
  113, 114, 117, 121, 123, 125, 126, 126, 126, 127, 127, 128, 128, 129, 132,
  132, 132, 134, 135, 137, 144, 145, 145, 146, 153, 155, 155, 155, 157, 163,
  165, 171, 172, 179, 181, 186, 190, 212
)
mean_resultant_length(finland_stria, w = NULL, na.rm = FALSE) # 0.800
#> [1] 0.8003694
```
