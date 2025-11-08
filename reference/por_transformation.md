# Conversion between spherical PoR to geographical coordinate system

Transformation from spherical PoR to geographical coordinate system and
vice versa

## Usage

``` r
geographical_to_PoR(x, PoR)

PoR_to_geographical(x, PoR)
```

## Arguments

- x:

  Can be either a `"data.frame"` containing `lat` and `lon` coordinates
  of a point in the geographical CRS or the `lat.PoR`, `lon.PoR`) of the
  point in the PoR CRS, a two-column matrix containing the lat and lon
  coordinates, a `sf` object, or a `raster` object.

- PoR:

  Pole of Rotation. `"data.frame"` containing the geographical
  coordinates of the Euler pole

## Value

object of same type of `x` with the transformed coordinates. If `x` was
a `data.frame`, transformed coordinates are named `lat.PoR` and
`lon.PoR` for PoR CRS, or `lat` and `lon` for geographical CRS).

## Examples

``` r
data("nuvel1")
por <- subset(nuvel1, nuvel1$plate.rot == "na") # North America relative to Pacific plate
data("san_andreas")
san_andreas.por <- geographical_to_PoR(san_andreas, por)
head(san_andreas.por)
#> Simple feature collection with 6 features and 9 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: -85.74431 ymin: 59.05548 xmax: -78.16683 ymax: 63.39714
#> Geodetic CRS:  unnamed
#>         id   lat     lon azi unc type depth quality regime
#> 1 wsm00892 38.14 -118.84  50  25  FMS     7       C      S
#> 2 wsm00893 35.97 -114.71  54  25  FMS     5       C      S
#> 3 wsm00894 37.93 -118.17  24  25  FMS     5       C      S
#> 4 wsm00896 38.63 -118.21  41  25  FMS    17       C      N
#> 5 wsm00897 39.08 -115.62  30  25  FMS     5       C      N
#> 6 wsm00903 38.58 -112.58  27  25  FMS     7       C      N
#>                     geometry
#> 1 POINT (-85.46364 59.05548)
#> 2 POINT (-78.16683 60.50445)
#> 3 POINT (-84.55317 59.38025)
#> 4 POINT (-85.74431 59.73632)
#> 5 POINT (-84.31053 61.68018)
#> 6 POINT (-80.60932 63.39714)
head(PoR_to_geographical(san_andreas.por, por))
#> Simple feature collection with 6 features and 9 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: -118.84 ymin: 35.97 xmax: -112.58 ymax: 39.08
#> Geodetic CRS:  WGS 84
#>         id   lat     lon azi unc type depth quality regime
#> 1 wsm00892 38.14 -118.84  50  25  FMS     7       C      S
#> 2 wsm00893 35.97 -114.71  54  25  FMS     5       C      S
#> 3 wsm00894 37.93 -118.17  24  25  FMS     5       C      S
#> 4 wsm00896 38.63 -118.21  41  25  FMS    17       C      N
#> 5 wsm00897 39.08 -115.62  30  25  FMS     5       C      N
#> 6 wsm00903 38.58 -112.58  27  25  FMS     7       C      N
#>                geometry
#> 1 POINT (-118.84 38.14)
#> 2 POINT (-114.71 35.97)
#> 3 POINT (-118.17 37.93)
#> 4 POINT (-118.21 38.63)
#> 5 POINT (-115.62 39.08)
#> 6 POINT (-112.58 38.58)
```
