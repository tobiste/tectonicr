# Smoothed SHmax on a grid after Hansen & Mount (1990) / Müller et al. (2003)

Spatial smoothing of SHmax orientations with a distance-weighted, robust
smoother. `lambda` balances smoothness against fit to the data. Two
kernels are available: `"fsr"` (fixed search radius) and `"fnn"` (fixed
number of nearest neighbours).

## Usage

``` r
circular_smooth(
  x,
  grid = NULL,
  lon_range = NULL,
  lat_range = NULL,
  gridsize = 2,
  method = c("fsr", "fnn"),
  R_range = 300,
  lambda = 0.5,
  n_iter = 20,
  tol = 1e-06,
  robust_scale_deg = 10,
  min_data = 3L,
  min_sq = 0,
  min_dist_threshold = Inf,
  quality_weighting = FALSE,
  qp = 1
)

PoR_circular_smooth(
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

- method:

  Kernel type: `"fsr"` (fixed search radius) or `"fnn"` (fixed number of
  nearest neighbours).

- R_range:

  numeric. For `method = "fsr"`: search radii in km. For
  `method = "fnn"`: number of nearest neighbours. Can also be a vector
  of different radii.

- lambda:

  numeric. Smoothness parameter. `lambda > 1` favours smoothness;
  `lambda < 1` favours fidelity to the data.

- n_iter:

  integer. Maximum number of fitting iterations.

- tol:

  numeric. Convergence tolerance on the doubled angle (radians).

- robust_scale_deg:

  numeric. Minimum robustness scale in degrees. The robust cutoff is 6
  times the larger of this value and the median residual, so
  observations further than about `6 * robust_scale_deg` from the fit
  are downweighted. Protects real local rotations (e.g. along faults)
  from being removed as outliers. Default is `10`.

- min_data:

  integer. If the number of observations within distance `R_range` is
  less than `min_data`, a missing value `NA` will be generated. Default
  is `3` for
  [`circular_idw()`](https://tobiste.github.io/tectonicr/reference/idw.md)
  and `4` for
  [`circular_idw_stats()`](https://tobiste.github.io/tectonicr/reference/idw.md).

- min_sq:

  numeric. Minimum sum of quality x robustness weights of the
  contributing observations; grid nodes below are dropped.

- min_dist_threshold:

  numeric. Distance threshold for smallest distance of the prediction
  location to the next observation location. Default is `200` km.

- quality_weighting:

  logical. If a quality weighting should be applied: Default is `TRUE`.
  If `FALSE`, overwrites `qp`.

- qp:

  numeric. The weighting power of the data quality (the higher the
  value, the more weight). Default is `1`. When set to `0`, no weighting
  is applied.

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

  Arguments passed to `circular_smooth()`

## Value

`sf` object with columns

- lon, lat:

  longitude and latitude in degrees

- azi:

  smoothed direction in degrees (0 to 180)

- sd:

  circular standard deviation in degrees

- R:

  search radius in km (for FNN, the effective radius at each node)

- mdr:

  mean distance / R

- N:

  number of contributing observations

- sq:

  sum of quality x robustness weights

Attribute `"fit"` is a list (one element per `R_range`) with the
observation-level fit: `lon`, `lat`, `azi_obs`, `azi_fit`, `q` (quality,
scaled so the best observation is 1) and `delta` (robustness weight; low
values mark downweighted outliers).

## Details

Based on Hansen & Mount (1990) and Müller, Wehrle, Hettel, Sperner &
Fuchs (2003). The predicted direction at a grid node is the weighted
circular mean of the fitted observation axes, computed with the same
routine as `stat = "mean"`.

## Note

Only axial (180° periodic) data are supported.

## References

Hansen, K. M., & Mount, V. S. (1990). Smoothing and extrapolation of
crustal stress orientation measurements. Journal of Geophysical
Research: Solid Earth, 95(B2), 1155–1165.
https://doi.org/10.1029/JB095iB02p01155

Müller, B., Wehrle, V., Hettel, S., Sperner, B., & Fuchs, K. (2003). A
new method for smoothing orientated data and its application to stress
data. Geological Society, London, Special Publications, 209(1), 107–126.
https://doi.org/10.1144/GSL.SP.2003.209.01.11

## See also

[`compact_grid()`](https://tobiste.github.io/tectonicr/reference/compact-grid.md)

Other smoothing:
[`PoR_idw`](https://tobiste.github.io/tectonicr/reference/PoR_idw.md),
[`idw`](https://tobiste.github.io/tectonicr/reference/idw.md),
[`krige`](https://tobiste.github.io/tectonicr/reference/krige.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data("san_andreas")
res <- circular_smooth(san_andreas, method = "fsr",
                          R_range = 300, lambda = 0.5)
head(res)
attributes(res)
} # }
```
