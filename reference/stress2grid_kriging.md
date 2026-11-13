# Spatial kriging of SHmax axes after van den Boogaart & Schaeben (2002)

Ordinary kriging of the doubled-angle embedding of SHmax axes, with a
distance-only trace variogram (isotropy of measurements). The
per-observation uncertainty, given as a standard deviation in degrees,
is used as measurement error on the diagonal of the kriging matrix.

## Usage

``` r
stress2grid_kriging(
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
  use_uncertainty = TRUE,
  axial = TRUE
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
  [`stress2grid()`](https://tobiste.github.io/tectonicr/reference/stress2grid.md)
  and `4` for
  [`stress2grid_stats()`](https://tobiste.github.io/tectonicr/reference/stress2grid.md).

- min_dist_threshold:

  numeric. Distance threshold for smallest distance of the prediction
  location to the next observation location. Default is `200` km.

- vgm_model:

  Variogram model: `"spherical"` (default), `"exponential"` or
  `"gaussian"`. The sill is fixed at 1 (see details).

- nugget:

  numeric in \\\[0, 1)\\. Nugget of the trace variogram.

- range_km:

  numeric. Range of the variogram model in km.

- use_uncertainty:

  logical. If `TRUE` (default), `unc` (standard deviation in degrees) is
  used as measurement error.

- axial:

  Must be `TRUE`. Only axial data are supported.

## Value

`sf` object containing

- lon,lat:

  longitude and latitude in degrees

- azi:

  kriged axis direction in degrees (0 to 180), from the dominant
  eigenvector of the orientation tensor

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

  confidence in \\\[0, 1\]\\: length of the kriged doubled-angle vector,
  capped at 1. Equals \\2 \times \text{eig}\_1 - 1\\

- kvar:

  kriging variance on the trace-variogram scale

- eig1:

  larger eigenvalue of the kriged orientation tensor, in \\\[0.5, 1\]\\.
  The smaller eigenvalue is `1 - eig1`

Attribute `"vgm"` records the variogram parameters used.

## Details

The method is the simplified ordinary-kriging form derived by van den
Boogaart & Schaeben (2002) for axes in the embedding u u^T. Their
isotropy-of-measurements assumption gives the trace variogram used here.
The matrix-valued estimate is projected back onto the axes by taking the
dominant eigenvector, which equals the direction of the predicted
doubled-angle vector. The length of that vector is the confidence
measure described in the paper.

Fit `nugget` and `range_km` to the empirical variogram from
[`variogram_axial()`](https://tobiste.github.io/tectonicr/reference/variogram_axial.md).
The variogram has sill 1, because the trace of the covariance of unit
vectors is 1 (eq. 36 of the paper).

## References

van den Boogaart, K.G. and Schaeben, H. (2002). Kriging of regionalized
directions, axes, and orientations I. Mathematical Geology, 34(5),
479-503.

## Examples

``` r
if (FALSE) { # \dontrun{
data("san_andreas")
v <- variogram_axial(san_andreas)
res <- stress2grid_kriging(san_andreas, R_range = 300,
                           vgm_model = "spherical",
                           nugget = 0.1, range_km = 200)
head(res)
} # }
```
