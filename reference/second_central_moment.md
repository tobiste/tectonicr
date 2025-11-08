# Second Central Momentum

Measures the skewness (a measure of the asymmetry of the probability
distribution) and the kurtosis (measure of the "tailedness" of the
probability distribution). Standardized versions are the skewness and
kurtosis normalized by the mean resultant length (Mardia 1972).

## Usage

``` r
second_central_moment(x, w = NULL, axial = TRUE, na.rm = FALSE)
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

list containing

- `skewness`:

  second central sine momentum, i.e. the skewness

- `std_skewness`:

  standardized skewness

- `kurtosis`:

  second central cosine momentum, i.e. the kurtosis

- `std_kurtosis`:

  standardized kurtosis

## Details

Negative values of skewness indicate skewed data in counterclockwise
direction.

Large kurtosis values indicate tailed, values close to `0` indicate
packed data.

## Examples

``` r
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
sa.por <- PoR_shmax(san_andreas, PoR, "right")
second_central_moment(sa.por$azi.PoR)
#> $skewness
#> [1] -0.02661796
#> 
#> $std_skewness
#> [1] -0.1758456
#> 
#> $kurtosis
#> [1] 0.3800558
#> 
#> $std_kurtosis
#> [1] 1.453805
#> 
second_central_moment(sa.por$azi.PoR, w = weighting(san_andreas$unc))
#> $skewness
#> [1] -0.04106606
#> 
#> $std_skewness
#> [1] -0.2712938
#> 
#> $kurtosis
#> [1] 0.3998642
#> 
#> $std_kurtosis
#> [1] 1.699345
#> 
```
