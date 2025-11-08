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

  :   SHmax in degree

  unc

  :   Uncertainties of SHmax in degree

  type

  :   Methods used for the determination of the orientation of SHmax

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
  is 2.5 (is ignored if `grid` is specified)

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

  geographical mean SHmax in degree

- azi.PoR:

  PoR mean SHmax in degree. Only if `remove_PoR=TRUE`

- sd:

  Standard deviation of SHmax in degrees

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
#>     lon.PoR  lat.PoR  azi.PoR       sd   R   N       mdr
#> 1 -84.77055 52.59628 161.9672 56.71303 400   9 0.8735345
#> 2 -84.77055 52.59628 163.2300 51.03117 450  46 0.9219984
#> 3 -84.77055 52.59628 158.1162 50.46203 500  85 0.8866905
#> 4 -84.77055 52.59628 158.1101 50.39348 550 187 0.8907467
#> 5 -84.77055 52.59628 152.7236 51.47597 600 298 0.8673966
#> 6 -84.77055 52.59628 148.7594 50.96027 650 385 0.8371344
#>                     geometry      lat       lon      azi
#> 1 POINT (-125.0802 34.07892) 34.07892 -125.0802 34.47049
#> 2 POINT (-125.0802 34.07892) 34.07892 -125.0802 35.73332
#> 3 POINT (-125.0802 34.07892) 34.07892 -125.0802 30.61948
#> 4 POINT (-125.0802 34.07892) 34.07892 -125.0802 30.61337
#> 5 POINT (-125.0802 34.07892) 34.07892 -125.0802 25.22692
#> 6 POINT (-125.0802 34.07892) 34.07892 -125.0802 21.26265

if (FALSE) { # \dontrun{
PoR_stress2grid_stats(san_andreas, PoR, mode = TRUE) |> head()
} # }
```
