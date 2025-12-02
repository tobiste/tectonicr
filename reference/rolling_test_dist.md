# Apply Rolling Functions using Circular Statistics

**\[superseded\]** A generic function for applying a function to rolling
margins of an array along an additional value.

## Usage

``` r
distroll_circstats(
  x,
  distance,
  FUN,
  width = NULL,
  min_n = 2,
  align = c("right", "center", "left"),
  w = NULL,
  sort = TRUE,
  ...
)

distroll_confidence(
  x,
  distance,
  w = NULL,
  width = NULL,
  min_n = 2,
  align = c("right", "center", "left"),
  sort = TRUE,
  ...
)

distroll_dispersion(
  x,
  y,
  w = NULL,
  w.y = NULL,
  distance,
  width = NULL,
  min_n = 2,
  align = c("right", "center", "left"),
  sort = TRUE,
  ...
)

distroll_dispersion_sde(
  x,
  y,
  w = NULL,
  w.y = NULL,
  distance,
  width = NULL,
  min_n = 2,
  align = c("right", "center", "left"),
  sort = TRUE,
  ...
)
```

## Arguments

- x, y:

  vectors of numeric values in degrees. `length(y)` is either 1 or
  `length(x)`

- distance:

  numeric. the independent variable along the values in `x` are sorted,
  e.g. the plate boundary distances

- FUN:

  the function to be applied

- width:

  numeric. the range across `distance` on which `FUN` should be applied
  on `x`. If `NULL`, then width is a number that separates the distances
  in 10 equal groups.

- min_n:

  integer. The minimum values that should be considered in `FUN` (2 by
  default), otherwise `NA`.

- align:

  specifies whether the index of the result should be left- or
  right-aligned or centered (default) compared to the rolling window of
  observations. This argument is only used if width represents widths.

- w:

  numeric. the weighting for `x`

- sort:

  logical. Should the values be sorted after `distance` prior to
  applying the function (`TRUE` by default).

- ...:

  optional arguments to `FUN`

- w.y:

  numeric. the weighting for `y`

## Value

two-column vectors of (sorted) `x` and the rolled statistics along
`distance`.

## Note

`distroll_circstats()` and friends are complete, and for new code it is
recommended switching to
[`distance_binned_stats()`](https://tobiste.github.io/tectonicr/reference/distance_binned_stats.md),
which is fasrter, easier to use, more featureful, and still under active
development.

## Examples

``` r
data("plates")
plate_boundary <- subset(plates, plates$pair == "na-pa")
data("san_andreas")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
san_andreas$distance <- distance_from_pb(
  x = san_andreas,
  PoR = PoR,
  pb = plate_boundary,
  tangential = TRUE
)
dat <- san_andreas |> cbind(PoR_shmax(san_andreas, PoR, "right"))

distroll_circstats(dat$azi.PoR,
  distance = dat$distance,
  w = 1 / dat$unc, FUN = circular_mean
) |> head()
#> Warning: `distroll_circstats()` was deprecated in tectonicr 0.4.4.9009.
#> ℹ Please use `distance_binned_stats()` instead.
#>       distance         x   n
#> [1,] -9.226505  68.11323   9
#> [2,] -7.604482 127.98953  19
#> [3,] -5.982458 155.84338  14
#> [4,] -4.360435 154.02577  39
#> [5,] -2.738412 138.55072 125
#> [6,] -1.116388 138.29357 759
distroll_confidence(dat$azi.PoR, distance = dat$distance, w = 1 / dat$unc) |> head()
#>       distance          x   n
#> [1,] -9.226505 180.000000   9
#> [2,] -7.604482  82.942137  19
#> [3,] -5.982458 180.000000  14
#> [4,] -4.360435  55.581771  39
#> [5,] -2.738412  13.587015 125
#> [6,] -1.116388   5.914822 759
distroll_dispersion(dat$azi.PoR,
  y = 135,
  distance = dat$distance, w = 1 / dat$unc
) |> head()
#>       distance          x   n
#> [1,] -9.226505 0.56974161   9
#> [2,] -7.604482 0.31555788  19
#> [3,] -5.982458 0.39396588  14
#> [4,] -4.360435 0.32627346  39
#> [5,] -2.738412 0.06582746 125
#> [6,] -1.116388 0.09488957 759
distroll_dispersion_sde(dat$azi.PoR,
  y = 135,
  distance = dat$distance, w = 1 / dat$unc, R = 100
) |> head()
#>       distance          x   n
#> [1,] -9.226505 0.20502401   9
#> [2,] -7.604482 0.16781628  19
#> [3,] -5.982458 0.20699430  14
#> [4,] -4.360435 0.12168181  39
#> [5,] -2.738412 0.02105609 125
#> [6,] -1.116388 0.01025607 759

# New functions
distance_binned_stats(
  dat$azi.PoR,
  distance = dat$distance, width.breaks = 1, unc = dat$unc, prd = 135
) |> head()
#> # A tibble: 6 × 19
#>   bins      n distance_min distance_median distance_max  mean    sd    var    lq
#>   <fct> <int>        <dbl>           <dbl>        <dbl> <dbl> <dbl>  <dbl> <dbl>
#> 1 [-9.…     2        -9.23           -9.09        -8.95   NA   NA   NA      NA  
#> 2 (-8.…     8        -8.26           -7.91        -7.56  102.  44.7  0.703  71.4
#> 3 (-7.…     9        -7.47           -6.84        -6.60  134.  36.2  0.550  92.2
#> 4 (-6.…    11        -6.44           -6.18        -5.87  150.  29.4  0.410 126. 
#> 5 (-5.…    10        -5.32           -5.06        -4.64  149.  37.9  0.583  65.5
#> 6 (-4.…    23        -4.41           -3.91        -3.57  159.  33.7  0.500  80.8
#> # ℹ 10 more variables: quasimedian <dbl>, uq <dbl>, median <dbl>, mode <dbl>,
#> #   CI <dbl>, skewness <dbl>, kurtosis <dbl>, nchisq <dbl>, dispersion <dbl>,
#> #   dispersion_sde <dbl>
```
