# Bootstrapped Estimates for Circular Dispersion

Calculates bootstrapped estimates of the circular dispersion, its
standard error and its confidence interval.

## Usage

``` r
circular_dispersion_boot(
  x,
  y = NULL,
  w = NULL,
  w.y = NULL,
  R = 1000,
  conf.level = 0.95,
  ...
)
```

## Arguments

- x, y:

  vectors of numeric values in degrees. `length(y)` is either `1` or
  `length(x)`

- w, w.y:

  (optional) Weights. A vector of positive numbers and of the same
  length as `x`. `w.y` is the (optional) weight of `y`.

- R:

  positive integer. The number of bootstrap replicates (1000 by
  default).

- conf.level:

  Level of confidence: \\(1 - \alpha \\)/100\\. (`0.95` by default).

- ...:

  optional arguments passed to
  [`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html)

## Value

list containing:

- `MLE`:

  the maximum likelihood estimate of the circular dispersion

- `sde`:

  standard error of MLE

- `CI`:

  lower and upper limit of the confidence interval of MLE

## See also

[`circular_dispersion()`](https://tobiste.github.io/tectonicr/reference/dispersion.md)

## Examples

``` r
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
circular_dispersion(sa.por$azi.PoR, y = 135, w = weighting(san_andreas$unc))
#> [1] 0.1384805
circular_dispersion_boot(sa.por$azi.PoR, y = 135, w = weighting(san_andreas$unc), R = 1000)
#> $MLE
#> [1] 0.2610566
#> 
#> $sde
#> [1] 0.01132889
#> 
#> $CI
#> [1] 0.2375632 0.2833176
#> 
```
