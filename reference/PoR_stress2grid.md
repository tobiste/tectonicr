# Spatial Interpolation of SHmax in PoR Coordinate Reference System

The data is transformed into the PoR system before the interpolation.
The interpolation grid is returned in geographical coordinates and
azimuths.

## Usage

``` r
PoR_stress2grid(
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

PoR_stress2grid_stats(
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

  :   \\\sigma\_\text{Hmax}\\ in degree

  unc

  :   Uncertainties of \\\sigma\_\text{Hmax}\\ in degree

  type

  :   Methods used for the determination of the orientation of
      \\\sigma\_\text{Hmax}\\

- PoR:

  Pole of Rotation. `data.frame` or object of class `"euler.pole"`
  containing the geographical coordinates of the Euler pole

- grid:

  (optional) Point object of class `sf`.

- PoR_grid:

  logical. Whether the grid should be generated based on the coordinate
  range in the PoR (`TRUE`, the default) CRS or the geographical CRS
  (`FALSE`). Is ignored if `grid` is specified.

- lon_range, lat_range:

  (optional) numeric vector specifying the minimum and maximum
  longitudes and latitudes (are ignored if `grid` is specified).

- gridsize:

  Numeric. Target spacing of the regular grid in decimal degree. Default
  is `2.5` (is ignored if `grid` is specified)

- remove_PoR:

  logical. Whether PoR azimuths and coordinates will be removed from
  final output or not (the default.)

- ...:

  Arguments passed to
  [`stress2grid()`](https://tobiste.github.io/tectonicr/reference/stress2grid.md)

## Value

`sf` object containing

- lon,lat:

  longitude and latitude in geographical CRS (in degrees)

- lon.PoR,lat.PoR:

  longitude and latitude in PoR CRS (in degrees). Only if
  `remove_PoR=TRUE`

- azi:

  geographical mean \\\sigma\_\text{Hmax}\\ in degree

- azi.PoR:

  PoR mean \\\sigma\_\text{Hmax}\\ in degree. Only if `remove_PoR=TRUE`

- sd:

  Standard deviation of \\\sigma\_\text{Hmax}\\ in degrees

- R:

  Search radius in km

- mdr:

  Mean distance of datapoints per search radius

- N:

  Number of data points in search radius

## Details

Stress field and wavelength analysis in PoR system and back-transformed

## See also

[`stress2grid()`](https://tobiste.github.io/tectonicr/reference/stress2grid.md),
[`compact_grid()`](https://tobiste.github.io/tectonicr/reference/compact-grid.md)

## Examples

``` r
data("san_andreas")
data("nuvel1")
PoR <- subset(nuvel1, nuvel1$plate.rot == "na")
PoR_stress2grid(san_andreas, PoR) |> head()
#> Simple feature collection with 6 features and 10 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: -125.0802 ymin: 34.07892 xmax: -125.0802 ymax: 34.07892
#> Geodetic CRS:  WGS 84
#>     lon.PoR  lat.PoR azi.PoR sd   R N mdr                   geometry      lat
#> 1 -84.77055 52.59628      NA  0  50 0  NA POINT (-125.0802 34.07892) 34.07892
#> 2 -84.77055 52.59628      NA  0 100 0  NA POINT (-125.0802 34.07892) 34.07892
#> 3 -84.77055 52.59628      NA  0 150 0  NA POINT (-125.0802 34.07892) 34.07892
#> 4 -84.77055 52.59628      NA  0 200 1  NA POINT (-125.0802 34.07892) 34.07892
#> 5 -84.77055 52.59628      NA  0 250 1  NA POINT (-125.0802 34.07892) 34.07892
#> 6 -84.77055 52.59628      NA  0 300 2  NA POINT (-125.0802 34.07892) 34.07892
#>         lon azi
#> 1 -125.0802  NA
#> 2 -125.0802  NA
#> 3 -125.0802  NA
#> 4 -125.0802  NA
#> 5 -125.0802  NA
#> 6 -125.0802  NA

if (FALSE) { # \dontrun{
PoR_stress2grid_stats(san_andreas, PoR, mode = TRUE) |> head()
} # }
```
