# Circular Summary Statistics

Circular mean, standard deviation, variance, quasi-quantiles, mode, 95%
confidence angle, standardized skewness and kurtosis

## Usage

``` r
circular_summary(
  x,
  w = NULL,
  axial = TRUE,
  mode = FALSE,
  kappa = NULL,
  fisher.CI = FALSE,
  conf.level = 0.95,
  na.rm = FALSE
)
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

- mode:

  logical. Whether the circular mode should be calculated or not.

- kappa:

  numeric. von Mises distribution concentration parameter used for the
  circular mode. Will be estimated using
  [`est.kappa()`](https://tobiste.github.io/tectonicr/reference/est.kappa.md)
  if not provided.

- fisher.CI:

  logical. Whether Fisher's or the default Mardia/Batchelet's confidence
  interval should be calculated.

- conf.level:

  numeric. Level of confidence: \\(1 - \alpha \\)/100\\. (`0.95` by
  default).

- na.rm:

  logical value indicating whether `NA` values in `x` should be stripped
  before the computation proceeds.

## Value

named vector

## See also

[`circular_mean()`](https://tobiste.github.io/tectonicr/reference/circle_stats.md),
[`circular_sd()`](https://tobiste.github.io/tectonicr/reference/circle_stats.md),
[`circular_var()`](https://tobiste.github.io/tectonicr/reference/circle_stats.md),
[`circular_quantiles()`](https://tobiste.github.io/tectonicr/reference/circle_stats.md),
[`confidence_angle()`](https://tobiste.github.io/tectonicr/reference/confidence.md),
[`second_central_moment()`](https://tobiste.github.io/tectonicr/reference/second_central_moment.md),
[`circular_mode()`](https://tobiste.github.io/tectonicr/reference/circular_mode.md)

## Examples

``` r
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
circular_summary(sa.por$azi.PoR)
#>            n         mean           sd          var          25% quasi-median 
#> 1126.0000000  140.8777709   23.4183183    0.2840287  124.8495460  136.8274701 
#>          75%       median           CI     skewness     kurtosis            R 
#>  150.2163335  138.9450918    5.4487325   -0.1758456    1.4538046    0.7159713 
circular_summary(sa.por$azi.PoR, w = weighting(san_andreas$unc))
#>            n         mean           sd          var          25% quasi-median 
#> 1126.0000000  140.8843069   22.3029729    0.2614357  124.8483102  136.8328725 
#>          75%       median           CI     skewness     kurtosis            R 
#>  150.2174282  138.9450918    5.3646843   -0.2712938    1.6993452    0.7385643 
```
