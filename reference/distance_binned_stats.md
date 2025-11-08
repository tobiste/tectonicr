# Distance Binned Summary Statistics

**\[experimental\]** Circular summary statistics over intervals of
distances.

## Usage

``` r
distance_binned_stats(
  azi,
  distance,
  n.breaks = 10,
  width.breaks = NULL,
  unc = NULL,
  prd = NULL,
  prd.error = NULL,
  kappa = 2,
  R = 1000,
  conf.level = 0.95,
  ...
)
```

## Arguments

- azi:

  numeric. Azimuth values in degrees.

- distance:

  numeric. the independent variable along the values in `azi` are
  sorted, e.g. the plate boundary distances

- n.breaks:

  numeric. number (greater than or equal to 2) giving the number of
  equal-sized intervals into which `distance` is to be cut. Default
  is 10. Will be ignored if `width.breaks` is specified.

- width.breaks:

  numeric. The width of the intervals into which `distance` is to be
  cut.

- unc:

  (optional) Uncertainties of `azi` (in degrees) acting as inverse
  weighting factors for statistics.

- prd:

  (optional) numeric. A predicted orientation in degrees.

- prd.error:

  (optional) numeric. The uncertainty of the predicted orientation in
  degrees.

- kappa:

  numeric. Concentration parameter applied for the circular mode.

- R:

  integer. Number of bootstrap iterates for estimating the error of the
  dispersion.

- conf.level:

  The level of confidence for confidence interval and bootstrapped
  standard error of dispersion.

- ...:

  optional arguments passed to
  [`ggplot2::cut_number()`](https://ggplot2.tidyverse.org/reference/cut_interval.html)
  and
  [`ggplot2::cut_width()`](https://ggplot2.tidyverse.org/reference/cut_interval.html)

## Value

tibble containing the `n` values for `azi`in each bin, min/median/max
distance of the bin, and the summary statistics for `azi`. If `prd` is
specified, the normal Chi-squared statistic, dispersion and its standard
error are returned as well.

## See also

[`circular_summary()`](https://tobiste.github.io/tectonicr/reference/circular_summary.md),
[`circular_dispersion()`](https://tobiste.github.io/tectonicr/reference/dispersion.md),
and
[`circular_dispersion_boot()`](https://tobiste.github.io/tectonicr/reference/circular_dispersion_boot.md)

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

distance_binned_stats(dat$azi.PoR,
  distance = dat$distance, width.breaks = 1,
  unc = dat$unc, prd = 135
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
