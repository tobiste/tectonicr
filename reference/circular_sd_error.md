# Standard Error of Mean Direction of Circular Data

Measure of the chance variation expected from sample to sample in
estimates of the mean direction (after Mardia 1972). It is a parametric
estimate of the the circular standard error of the mean direction by the
particular form of the standard error for the von Mises distribution.
The approximated standard error of the mean direction is computed by the
mean resultant length and the MLE concentration parameter \\\kappa\\.

## Usage

``` r
circular_sd_error(x, w = NULL, axial = TRUE, na.rm = TRUE)
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

- Batschelet, E. (1971). Recent statistical methods for orientation
  data. "Animal Orientation, Symposium 1970 on Wallops Island". Amer.
  Inst. Biol. Sciences, Washington.

- Mardia, K.V. (1972). Statistics of Directional Data: Probability and
  Mathematical Statistics. London: Academic Press.

- N.I. Fisher (1993) Statistical Analysis of Circular Data, Cambridge
  University Press.

- Davis (1986) Statistics and data analysis in geology. 2nd ed., John
  Wiley & Sons.

## See also

[`mean_resultant_length()`](https://tobiste.github.io/tectonicr/reference/mean_resultant_length.md),
[`circular_mean()`](https://tobiste.github.io/tectonicr/reference/circle_stats.md)

## Examples

``` r
# Example data from Davis (1986), pp. 316
finland_stria <- c(
  23, 27, 53, 58, 64, 83, 85, 88, 93, 99, 100, 105, 113,
  113, 114, 117, 121, 123, 125, 126, 126, 126, 127, 127, 128, 128, 129, 132,
  132, 132, 134, 135, 137, 144, 145, 145, 146, 153, 155, 155, 155, 157, 163,
  165, 171, 172, 179, 181, 186, 190, 212
)
circular_sd_error(finland_stria, axial = FALSE)
#> [1] 0.09244732

data(san_andreas)
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
circular_sd_error(sa.por$azi.PoR, w = weighting(san_andreas$unc))
#> [1] 0.02387728
```
