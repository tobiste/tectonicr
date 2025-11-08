# Confidence Interval around the Mean Direction of Circular Data after Batschelet (1971)

Probabilistic limit on the location of the true or population mean
direction, assuming that the estimation errors are normally distributed.

## Usage

``` r
confidence_angle(x, conf.level = 0.95, w = NULL, axial = TRUE, na.rm = TRUE)

confidence_interval(x, conf.level = 0.95, w = NULL, axial = TRUE, na.rm = TRUE)
```

## Arguments

- x:

  numeric vector. Values in degrees.

- conf.level:

  Level of confidence: \\(1 - \alpha \\)/100\\. (`0.95` by default).

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

Angle in degrees

## Details

The confidence angle gives the interval, i.e. plus and minus the
confidence angle, around the mean direction of a particular sample, that
contains the true mean direction under a given level of confidence.

## References

- Batschelet, E. (1971). Recent statistical methods for orientation
  data. "Animal Orientation, Symposium 1970 on Wallops Island". Amer.
  Inst. Biol. Sciences, Washington.

- Mardia, K.V. (1972). Statistics of Directional Data: Probability and
  Mathematical Statistics. London: Academic Press. (p. 146)

- Davis (1986) Statistics and data analysis in geology. 2nd ed., John
  Wiley & Sons.

- Jammalamadaka, S. Rao and Sengupta, A. (2001). Topics in Circular
  Statistics, Sections 3.3.3 and 3.4.1, World Scientific Press,
  Singapore.

## See also

[`mean_resultant_length()`](https://tobiste.github.io/tectonicr/reference/mean_resultant_length.md),
[`circular_sd_error()`](https://tobiste.github.io/tectonicr/reference/circular_sd_error.md)

## Examples

``` r
# Example data from Davis (1986), pp. 316
finland_stria <- c(
  23, 27, 53, 58, 64, 83, 85, 88, 93, 99, 100, 105, 113,
  113, 114, 117, 121, 123, 125, 126, 126, 126, 127, 127, 128, 128, 129, 132,
  132, 132, 134, 135, 137, 144, 145, 145, 146, 153, 155, 155, 155, 157, 163,
  165, 171, 172, 179, 181, 186, 190, 212
)
confidence_angle(finland_stria, axial = FALSE)
#> [1] 10.43928
confidence_interval(finland_stria, axial = FALSE)
#> $mu
#> [1] 129.1903
#> 
#> $conf.angle
#> [1] 10.43928
#> 
#> $conf.interval
#> [1] 118.7510 139.6296
#> 

data(san_andreas)
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
confidence_angle(sa.por$azi.PoR, w = weighting(san_andreas$unc))
#> [1] 5.364684
confidence_interval(sa.por$azi.PoR, w = weighting(san_andreas$unc))
#> $mu
#> [1] 140.8843
#> 
#> $conf.angle
#> [1] 5.364684
#> 
#> $conf.interval
#> [1] 135.5196 146.2490
#> 
```
