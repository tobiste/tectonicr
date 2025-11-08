# Adaptive Kernel Dispersion

Stress field and wavelength analysis using circular dispersion (or other
statistical estimators for dispersion)

## Usage

``` r
kernel_dispersion(
  x,
  stat = c("dispersion", "nchisq", "rayleigh"),
  grid = NULL,
  lon_range = NULL,
  lat_range = NULL,
  gridsize = 2.5,
  min_data = 3L,
  max_data = Inf,
  min_dist_threshold = 200,
  dist_threshold = 0.1,
  stat_threshold = Inf,
  R_range = seq(100, 2000, 100),
  ...
)

dispersion_grid(...)
```

## Arguments

- x:

  `sf` object containing

  azi

  :   Azimuth in degree

  unc

  :   Uncertainties of azimuth in degree

  prd

  :   Predicted value for azimuth

- stat:

  The measurement of dispersion to be calculated. Either `"dispersion"`
  (default), `"nchisq"`, or `"rayleigh"` for circular dispersion,
  normalized Chi-squared test statistic, or Rayleigh test statistic.

- grid:

  (optional) Point object of class `sf`.

- lon_range, lat_range:

  (optional) numeric vector specifying the minimum and maximum
  longitudes and latitudes (are ignored if `"grid"` is specified).

- gridsize:

  Numeric. Target spacing of the regular grid in decimal degree. Default
  is 2.5. (is ignored if `"grid"` is specified)

- min_data:

  Integer. Minimum number of data per bin. Default is `3`

- max_data:

  integer. The number of nearest observations that should be used for
  prediction, where "nearest" is defined in terms of the space of the
  spatial locations. Default is `Inf`.

- min_dist_threshold:

  Numeric. Maximum distance (in km) of the grid point to the next data
  point. Default is 200

- dist_threshold:

  Numeric. Distance weight to prevent overweight of data nearby (`0` to
  `1`). Default is `0.1`

- stat_threshold:

  numeric. Generates missing values when the kernel `stat` value exceeds
  this threshold. Default is `Inf`.

- R_range:

  Numeric value or vector specifying the (adaptive) kernel half-width(s)
  as search radius (in km). Default is `seq(50, 1000, 50)`

- ...:

  optional arguments to
  [`dist_greatcircle()`](https://tobiste.github.io/tectonicr/reference/dist_greatcircle.md)

## Value

`sf` object containing

- lon,lat:

  longitude and latitude in degree

- stat:

  output of function defined in `stat`

- R:

  The rearch radius in km.

- mdr:

  Mean distance of datapoints per search radius

- N:

  Number of data points in search radius

## Note

`dispersion_grid()` was renamed to `kernel_dispersion()` to create a
more consistent API.

## See also

[`circular_dispersion()`](https://tobiste.github.io/tectonicr/reference/dispersion.md),
[`norm_chisq()`](https://tobiste.github.io/tectonicr/reference/norm_chisq.md),
[`rayleigh_test()`](https://tobiste.github.io/tectonicr/reference/rayleigh_test.md)

## Examples

``` r
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
san_andreas_por <- data2PoR(san_andreas, PoR)
san_andreas_por$prd <- 135
kernel_dispersion(san_andreas_por) |> head()
#> Simple feature collection with 6 features and 6 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: -84.77055 ymin: 52.59628 xmax: -84.77055 ymax: 52.59628
#> Geodetic CRS:  unnamed
#>         lon      lat      stat   R   N       mdr                   geometry
#> 1 -84.77055 52.59628        NA 100   0        NA POINT (-84.77055 52.59628)
#> 2 -84.77055 52.59628        NA 200   1        NA POINT (-84.77055 52.59628)
#> 3 -84.77055 52.59628        NA 300   2        NA POINT (-84.77055 52.59628)
#> 4 -84.77055 52.59628 0.7982977 400   9 0.8735345 POINT (-84.77055 52.59628)
#> 5 -84.77055 52.59628 0.8069604 500  85 0.8866905 POINT (-84.77055 52.59628)
#> 6 -84.77055 52.59628 0.7687554 600 298 0.8673966 POINT (-84.77055 52.59628)
```
