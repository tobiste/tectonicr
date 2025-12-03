# Spatial Interpolation of SHmax

Stress field interpolation and wavelength analysis using a kernel
(weighted) mean/median and standard deviation/IQR of stress data.
Parameters can be adjusted to have inverse-distance-weighting (IDW) or
nearest-neighbor interpolations (NN).

## Usage

``` r
stress2grid(
  x,
  stat = c("mean", "median", "tensor"),
  grid = NULL,
  lon_range = NULL,
  lat_range = NULL,
  gridsize = 2,
  min_data = 3L,
  max_data = Inf,
  max_sd = Inf,
  threshold = deprecated(),
  min_dist_threshold = 200,
  arte_thres = deprecated(),
  method_weighting = FALSE,
  quality_weighting = TRUE,
  dist_weighting = c("inverse", "linear", "none"),
  idp = 1,
  qp = 1,
  mp = 1,
  dist_threshold = 0.1,
  R_range = seq(50, 1000, 50),
  ...
)

stress2grid_stats(
  x,
  grid = NULL,
  lon_range = NULL,
  lat_range = NULL,
  gridsize = 2,
  min_data = 4L,
  max_data = Inf,
  threshold = deprecated(),
  min_dist_threshold = 200,
  arte_thres = deprecated(),
  method_weighting = FALSE,
  quality_weighting = TRUE,
  dist_weighting = c("inverse", "linear", "none"),
  idp = 1,
  qp = 1,
  mp = 1,
  dist_threshold = 0.1,
  R_range = seq(50, 1000, 50),
  mode = FALSE,
  kappa = 10,
  ...
)
```

## Arguments

- x:

  `sf` object containing

  azi

  :   SHmax in degree

  unc

  :   (optional) Uncertainties of SHmax in degree

  type

  :   (optional) Methods used for the determination of the direction of
      SHmax

- stat:

  whether the direction of interpolated SHmax is based on the circular
  mean and standard deviation (`"mean"`, the default), the
  quasi-circular median and quasi-interquartile range (`"median"`), or
  the orientation tensor based principal direction and dispersion
  ("tensor").

- grid:

  (optional) Point object of class `sf`.

- lon_range, lat_range:

  (optional) numeric vector specifying the minimum and maximum
  longitudes and latitudes (ignored if `grid` is specified).

- gridsize:

  numeric. Target spacing of the regular grid in decimal degree. Default
  is `2.5`. (is ignored if `grid` is specified)

- min_data:

  integer. If the number of observations within distance `R_range` is
  less than `min_data`, a missing value `NA` will be generated. Default
  is `3` for `stress2grid()` and `4` for `stress2grid_stats()`.

- max_data:

  integer. The number of nearest observations that should be used for
  prediction, where "nearest" is defined in terms of the space of the
  spatial locations. Default is `Inf`.

- max_sd:

  numeric. Threshold for deviation of direction in degrees; if exceeds,
  missing values will be generated.

- threshold:

  **\[deprecated\]** is no longer supported; use `max_sd` instead.

- min_dist_threshold:

  numeric. Distance threshold for smallest distance of the prediction
  location to the next observation location. Default is `200` km.

- arte_thres:

  **\[deprecated\]** is no longer supported; use `min_dist_threshold`
  instead.

- method_weighting:

  logical. If a method weighting should be applied: Default is `FALSE`.
  If `FALSE`, overwrites `mp`.

- quality_weighting:

  logical. If a quality weighting should be applied: Default is `TRUE`.
  If `FALSE`, overwrites `qp`.

- dist_weighting:

  Distance weighting method which should be used. One of `"none"`,
  `"linear"`, or `"inverse"` (the default).

- idp, qp, mp:

  numeric. The weighting power of inverse distance, quality and method
  (the higher the value, the more weight). Default is `1`. When set to
  `0`, no weighting is applied. Only effective when
  `dist_weighting=="inverse"`.

- dist_threshold:

  numeric. Distance weight to prevent overweight of data nearby (0 to
  1). Default is `0.1`

- R_range:

  numeric value or vector specifying the kernel half-width(s) search
  radii, i.e. the maximum distance from the prediction location to be
  used for prediction (in km). Default is `seq(50, 1000, 50)`. If
  combined with `max_data`, both criteria apply.

- ...:

  (optional) arguments to
  [`dist_greatcircle()`](https://tobiste.github.io/tectonicr/reference/dist_greatcircle.md)

- mode:

  logical. Should the circular mode be included in the statistical
  summary (slow)?

- kappa:

  numeric. von Mises distribution concentration parameter used for the
  circular mode. Will be estimated using
  [`est.kappa()`](https://tobiste.github.io/tectonicr/reference/est.kappa.md)
  if not provided.

## Value

`sf` object containing

- lon,lat:

  longitude and latitude in degrees

- azi:

  Circular mean od median SHmax in degree

- sd:

  Circular standard deviation or Quasi-IQR on the Circle of SHmax in
  degrees

- R:

  Search radius in km

- mdr:

  Mean distance between grid point and datapoints per search radius

- N:

  Number of data points in search radius

When `stress2grid_stats()`, `azi` and `sd` are replaced by the output of
[`circular_summary()`](https://tobiste.github.io/tectonicr/reference/circular_summary.md).

## Details

`stress2grid()` is originally based on the MATLAB script "stress2grid"
by Ziegler and Heidbach (2019):
<https://github.com/MorZieg/Stress2Grid>. The tectonicr version has been
significantly modified to provide better performance and more
flexibility.

`stress2grid_stats()` is based on `stress2grid()` but calculates
circular summary statistics (see
[`circular_summary()`](https://tobiste.github.io/tectonicr/reference/circular_summary.md)).

## References

Ziegler, M. and Heidbach, O. (2019). Matlab Script Stress2Grid v1.1. GFZ
Data Services.
[doi:10.5880/wsm.2019.002](https://doi.org/10.5880/wsm.2019.002)

## See also

[`dist_greatcircle()`](https://tobiste.github.io/tectonicr/reference/dist_greatcircle.md),
[`PoR_stress2grid()`](https://tobiste.github.io/tectonicr/reference/PoR_stress2grid.md),
[`compact_grid()`](https://tobiste.github.io/tectonicr/reference/compact-grid.md),
[`circular_mean()`](https://tobiste.github.io/tectonicr/reference/circle_stats.md),
[`circular_median()`](https://tobiste.github.io/tectonicr/reference/circle_stats.md),
[`circular_sd()`](https://tobiste.github.io/tectonicr/reference/circle_stats.md),
[`circular_summary()`](https://tobiste.github.io/tectonicr/reference/circular_summary.md)

## Examples

``` r
data("san_andreas")

# Inverse Distance Weighting interpolation:
stress2grid(san_andreas, stat = "median") |> head()
#> Simple feature collection with 6 features and 7 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: -112.82 ymin: 24.08 xmax: -112.82 ymax: 24.08
#> Geodetic CRS:  WGS 84
#>       lon   lat      azi sd   R   N       mdr              geometry
#> 1 -112.82 24.08 148.7449  7 150   4 0.8224433 POINT (-112.82 24.08)
#> 2 -112.82 24.08 148.7449  7 200   4 0.6168325 POINT (-112.82 24.08)
#> 3 -112.82 24.08 147.0000  7 250   7 0.6785179 POINT (-112.82 24.08)
#> 4 -112.82 24.08 163.0000  1 300  17 0.7787910 POINT (-112.82 24.08)
#> 5 -112.82 24.08 163.0000  0 350  73 0.8790139 POINT (-112.82 24.08)
#> 6 -112.82 24.08 165.0000  0 400 127 0.8402209 POINT (-112.82 24.08)

stress2grid(san_andreas, stat = "tensor") |> head()
#> Simple feature collection with 6 features and 7 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: -112.82 ymin: 24.08 xmax: -112.82 ymax: 24.08
#> Geodetic CRS:  WGS 84
#>       lon   lat        azi         sd   R   N       mdr              geometry
#> 1 -112.82 24.08 -15.084872 0.15607975 150   4 0.8224433 POINT (-112.82 24.08)
#> 2 -112.82 24.08 -15.084872 0.15607975 200   4 0.6168325 POINT (-112.82 24.08)
#> 3 -112.82 24.08 -22.098607 0.10914104 250   7 0.6785179 POINT (-112.82 24.08)
#> 4 -112.82 24.08 -15.989161 0.06416128 300  17 0.7787910 POINT (-112.82 24.08)
#> 5 -112.82 24.08  -9.946535 0.07734405 350  73 0.8790139 POINT (-112.82 24.08)
#> 6 -112.82 24.08  -6.041035 0.07223013 400 127 0.8402209 POINT (-112.82 24.08)

# Nearest Neighbor interpolation:
stress2grid(san_andreas, stat = "median", max_data = 5) |> head()
#> Simple feature collection with 6 features and 7 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: -120.82 ymin: 30.08 xmax: -120.82 ymax: 30.08
#> Geodetic CRS:  WGS 84
#>       lon   lat      azi sd   R N       mdr              geometry
#> 1 -120.82 30.08 35.00000  0 350 4 0.6995191 POINT (-120.82 30.08)
#> 2 -120.82 30.08 35.53572  0 400 5 0.6828824 POINT (-120.82 30.08)
#> 3 -120.82 30.08 35.53572  0 450 5 0.6070066 POINT (-120.82 30.08)
#> 4 -120.82 30.08 35.53572  0 500 5 0.5463060 POINT (-120.82 30.08)
#> 5 -120.82 30.08 35.53572  0 550 5 0.4966418 POINT (-120.82 30.08)
#> 6 -120.82 30.08 35.53572  0 600 5 0.4552550 POINT (-120.82 30.08)

if (FALSE) { # \dontrun{
stress2grid_stats(san_andreas) |> head()
} # }
```
