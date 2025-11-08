# Confidence Interval around the Mean Direction of Circular Data after Fisher (1993)

For large samples (`n >=25`) i performs are parametric estimate based on
[`sample_circular_dispersion()`](https://tobiste.github.io/tectonicr/reference/sample_dispersion.md).
For smaller size samples, it returns a bootstrap estimate.

## Usage

``` r
confidence_interval_fisher(
  x,
  conf.level = 0.95,
  w = NULL,
  axial = TRUE,
  na.rm = TRUE,
  boot = FALSE,
  R = 1000L,
  quiet = FALSE
)
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

- boot:

  logical. Force bootstrap estimation

- R:

  integer. number of bootstrap replicates

- quiet:

  logical. Prints the used estimation (parametric or bootstrap).

## Value

list

## References

N.I. Fisher (1993) Statistical Analysis of Circular Data, Cambridge
University Press.

## Examples

``` r
# Example data from Davis (1986), pp. 316
finland_stria <- c(
  23, 27, 53, 58, 64, 83, 85, 88, 93, 99, 100, 105, 113,
  113, 114, 117, 121, 123, 125, 126, 126, 126, 127, 127, 128, 128, 129, 132,
  132, 132, 134, 135, 137, 144, 145, 145, 146, 153, 155, 155, 155, 157, 163,
  165, 171, 172, 179, 181, 186, 190, 212
)
confidence_interval_fisher(finland_stria, axial = FALSE)
#> Parametric estimate
#> $mu
#> [1] 129.1903
#> 
#> $conf.angle
#> [1] 7.043583
#> 
#> $conf.interval
#> [1] 122.1467 136.2339
#> 
confidence_interval_fisher(finland_stria, axial = FALSE, boot = TRUE)
#> Bootstrap estimate based on 1000 replicates
#> $mu
#> [1] 129.2045
#> 
#> $conf.angle
#> [1] 5.787037
#> 
#> $conf.interval
#> [1] 123.7251 135.1838
#> 

data(san_andreas)
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
confidence_interval_fisher(sa.por$azi.PoR, w = weighting(san_andreas$unc))
#> Parametric estimate
#> $mu
#> [1] 140.8843
#> 
#> $conf.angle
#> [1] 1.711388
#> 
#> $conf.interval
#> [1] 139.1729 142.5957
#> 
confidence_interval_fisher(sa.por$azi.PoR, w = weighting(san_andreas$unc), boot = TRUE)
#> Bootstrap estimate based on 1000 replicates
#> $mu
#> [1] 140.8711
#> 
#> $conf.angle
#> [1] 0.721629
#> 
#> $conf.interval
#> [1] 140.1712 141.5932
#> 
```
