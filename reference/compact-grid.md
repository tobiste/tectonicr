# Compact Smoothed Stress Field

Filter smoothed stress field containing a range of search radii or
kernel half widths to find shortest wavelength (R) with the least
circular sd. or dispersion (or any statistic) for each coordinate,
respectively.

## Usage

``` r
compact_grid(x, type = c("stress", "dispersion"))

compact_grid2(x, ..., FUN = min)
```

## Arguments

- x:

  output of
  [`stress2grid()`](https://tobiste.github.io/tectonicr/reference/stress2grid.md),
  [`PoR_stress2grid()`](https://tobiste.github.io/tectonicr/reference/PoR_stress2grid.md),
  [`stress2grid_stats()`](https://tobiste.github.io/tectonicr/reference/stress2grid.md),
  or
  [`kernel_dispersion()`](https://tobiste.github.io/tectonicr/reference/kernel_dispersion.md)

- type:

  character. Type of the grid `x`. Either `"stress"` (when input is
  [`stress2grid()`](https://tobiste.github.io/tectonicr/reference/stress2grid.md)
  or
  [`PoR_stress2grid()`](https://tobiste.github.io/tectonicr/reference/PoR_stress2grid.md))
  or `"dispersion"` (when input is
  [`kernel_dispersion()`](https://tobiste.github.io/tectonicr/reference/kernel_dispersion.md)).

- ...:

  `<tidy-select>` One unquoted expression separated by commas. Variable
  names can be used as if they were positions in the data frame.
  Variable must be a column in `x`.

- FUN:

  function is used to aggregate the data using the search radius `R`.
  Default is [`min()`](https://rdrr.io/r/base/Extremes.html).

## Value

`sf` object

## See also

[`stress2grid()`](https://tobiste.github.io/tectonicr/reference/stress2grid.md),
[`PoR_stress2grid()`](https://tobiste.github.io/tectonicr/reference/PoR_stress2grid.md),
[`kernel_dispersion()`](https://tobiste.github.io/tectonicr/reference/kernel_dispersion.md),
[`stress2grid_stats()`](https://tobiste.github.io/tectonicr/reference/stress2grid.md),
[`dplyr::dplyr_tidy_select()`](https://dplyr.tidyverse.org/reference/dplyr_tidy_select.html)

## Examples

``` r
data("san_andreas")
res <- stress2grid(san_andreas)
compact_grid(res) |> head()
#> Simple feature collection with 6 features and 7 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: -108.82 ymin: 24.08 xmax: -108.82 ymax: 36.08
#> Geodetic CRS:  WGS 84
#>     R     lon   lat        azi       sd  N       mdr              geometry
#> 1  50 -108.82 24.08 179.958278 27.62559 22 0.4806675 POINT (-108.82 24.08)
#> 2 150 -108.82 26.08   0.197365 23.97752 37 0.8177757 POINT (-108.82 26.08)
#> 3 200 -108.82 30.08  25.078141 19.52286  3 0.6026966 POINT (-108.82 30.08)
#> 4 150 -108.82 32.08  22.579218 18.54572  4 0.6892362 POINT (-108.82 32.08)
#> 5 100 -108.82 34.08  30.225537 61.46356  6 0.5770263 POINT (-108.82 34.08)
#> 6 150 -108.82 36.08 145.432349 78.41604  3 0.7654557 POINT (-108.82 36.08)

if (FALSE) { # \dontrun{
res2 <- stress2grid_stats(san_andreas)
compact_grid2(res2, var, FUN = min)
} # }
```
