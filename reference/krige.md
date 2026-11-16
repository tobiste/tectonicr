# Spatial kriging of SHmax axes after van den Boogaart & Schaeben (2002)

Ordinary kriging of the axial orientation tensors, with a distance-only
trace variogram (isotropy of measurements). The kriged tensor is
decomposed with the eigenvalue method of
[`ot_eigen2d()`](https://tobiste.github.io/tectonicr/reference/ort-eigen.md).
The per-observation uncertainty, given as a standard deviation in
degrees, is used as measurement error on the diagonal of the kriging
matrix.

## Usage

``` r
circular_krige(
  x,
  grid = NULL,
  lon_range = NULL,
  lat_range = NULL,
  gridsize = 2,
  R_range = 300,
  max_data = 50,
  min_data = 3L,
  min_dist_threshold = 200,
  vgm_model = c("spherical", "exponential", "gaussian"),
  nugget = 0,
  range_km = 300,
  sill = NULL,
  use_uncertainty = TRUE
)

PoR_circular_krige(
  x,
  PoR,
  grid = NULL,
  PoR_grid = TRUE,
  lon_range = NULL,
  lat_range = NULL,
  gridsize = 2.5,
  remove_PoR = FALSE,
  ...
)

PoR_circular_krige(
  x,
  PoR,
  grid = NULL,
  PoR_grid = TRUE,
  lon_range = NULL,
  lat_range = NULL,
  gridsize = 2.5,
  remove_PoR = FALSE,
  ...
)
```

## Arguments

- x:

  `sf` object containing

  azi

  :   angles in degrees (e.g. orientation of \\\sigma\_\text{Hmax}\\)

  unc

  :   (optional) Uncertainties of `azi` in degree

  type

  :   (optional) Methods used for the determination of the direction of
      \\\sigma\_\text{Hmax}\\

- grid:

  (optional) Point object of class `sf`.

- lon_range, lat_range:

  (optional) numeric vector specifying the minimum and maximum
  longitudes and latitudes (ignored if `grid` is specified).

- gridsize:

  numeric. Target spacing of the regular grid in decimal degree. Default
  is `2.5`. (is ignored if `grid` is specified)

- R_range:

  numeric value or vector of search radii in km. Only observations
  within the radius are used.

- max_data:

  integer. The number of nearest observations used for each prediction.
  Default is `50`.

- min_data:

  integer. If the number of observations within distance `R_range` is
  less than `min_data`, a missing value `NA` will be generated. Default
  is `3` for
  [`circular_idw()`](https://tobiste.github.io/tectonicr/reference/idw.md)
  and `4` for
  [`circular_idw_stats()`](https://tobiste.github.io/tectonicr/reference/idw.md).

- min_dist_threshold:

  numeric. Distance threshold for smallest distance of the prediction
  location to the next observation location. Default is `200` km.

- vgm_model:

  Variogram model: `"spherical"` (default), `"exponential"` or
  `"gaussian"`.

- nugget:

  numeric in \\\[0, 1)\\. Nugget of the trace variogram.

- range_km:

  numeric. Range of the variogram model in km.

- sill:

  numeric in \\(\text{nugget}, 1\]\\. Sill of the trace variogram. If
  `NULL` (default), it is set to \\1 - \bar{R}^2\\, where \\\bar{R}\\ is
  the mean resultant length of the doubled angles of `x$azi`.

- use_uncertainty:

  logical. If `TRUE` (default), `unc` (standard deviation in degrees) is
  used as measurement error.

- PoR:

  Pole of Rotation. `data.frame` or object of class `"euler.pole"`
  containing the geographical coordinates of the Euler pole

- PoR_grid:

  logical. Whether the grid should be generated based on the coordinate
  range in the PoR (`TRUE`, the default) CRS or the geographical CRS
  (`FALSE`). Is ignored if `grid` is specified.

- remove_PoR:

  logical. Whether PoR azimuths and coordinates will be removed from
  final output or not (the default.)

- ...:

  Arguments passed to `circular_krige()`

## Value

`sf` object containing

- lon,lat:

  longitude and latitude in degrees

- azi:

  kriged axis direction in degrees (0 to 180), the principal direction
  from
  [`ot_eigen2d()`](https://tobiste.github.io/tectonicr/reference/ort-eigen.md)

- sd:

  circular standard deviation of the axis in degrees, computed from
  `conf` as `rad2deg(sqrt(-2 log(conf))) / 2`

- R:

  search radius in km

- mdr:

  mean distance between grid point and data / R

- N:

  number of observations used

- conf:

  confidence in \\\[0, 1\]\\, \\2 \lambda_1 - 1\\, i.e. the axial
  strength relative to the trace of the kriged tensor

- kvar:

  kriging variance on the trace-variogram scale

- eig1:

  largest eigenvalue of the kriged tensor (axial strength), in \\\[0.5,
  1\]\\. The smaller eigenvalue, \\1 = \lambda_1\\, is the axial
  dispersion

Attribute `"vgm"` records the variogram parameters used.

## Details

The method is the simplified ordinary-kriging form derived by van den
Boogaart & Schaeben (2002) for axes in the embedding \\u u^T\\. Their
isotropy-of-measurements assumption gives the trace variogram used here.
The kriged tensor is \\\sum_i \lambda_i u_i u_i^T\\, computed by
[`ot_eigen2d()`](https://tobiste.github.io/tectonicr/reference/ort-eigen.md)
with the kriging weights as weights. The principal direction is the
eigenvector of the largest eigenvalue, and the eigenvalues sum to one.
Their values give the strength and dispersion of the prediction.

Fit `nugget`, `sill`, and `range_km` to the empirical variogram from
[`variogram_axial()`](https://tobiste.github.io/tectonicr/reference/variogram_axial.md).

## Note

Only axial (180° periodic) data are supported.

## References

van den Boogaart, K.G. and Schaeben, H. (2002). Kriging of regionalized
directions, axes, and orientations I. Mathematical Geology, 34(5),
479-503.

## See also

[`variogram_axial()`](https://tobiste.github.io/tectonicr/reference/variogram_axial.md)
for computing the empirical variogram of axial data;
[`ortensor2d()`](https://tobiste.github.io/tectonicr/reference/ortensor2d.md)
and
[`ot_eigen2d()`](https://tobiste.github.io/tectonicr/reference/ort-eigen.md)
for the orientation tensor and its eigenvalue decomposition.

Other smoothing:
[`PoR_idw`](https://tobiste.github.io/tectonicr/reference/PoR_idw.md),
[`idw`](https://tobiste.github.io/tectonicr/reference/idw.md),
[`smooth`](https://tobiste.github.io/tectonicr/reference/smooth.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data("san_andreas")
v <- variogram_axial(san_andreas)
res <- circular_krige(san_andreas, R_range = 300,
                           vgm_model = "spherical",
                           nugget = 0.25, range_km = 1000)
head(res)
attributes(res)
} # }
```
