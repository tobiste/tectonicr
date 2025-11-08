# Transforms coordinates and azimuths into PoR coordinates system

Convenience function to add PoR coordinates and PoR azimuths to data

## Usage

``` r
data2PoR(x, PoR)
```

## Arguments

- x:

  `sf` object or a `data.frame` containing the coordinates of the
  point(s) (`lat`, `lon` columns). `x` must contain the direction of
  \\\sigma\_{Hmax}\\ as column `azi`, its standard deviation (column
  `unc`) is optional).

- PoR:

  `data.frame` or object of class `euler.pole` containing the
  geographical coordinates of the Eule pole.

## Value

`sf` object in PoR CRS with additional columns `lon.PoR`, `lat.PoR`, and
`azi.PoR`

## Examples

``` r
por <- subset(nuvel1, nuvel1$plate.rot == "na")
data2PoR(san_andreas, por)
#> Simple feature collection with 1126 features and 12 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: -91.02055 ymin: 51.34628 xmax: -53.77881 ymax: 66.80636
#> Geodetic CRS:  unnamed
#> First 10 features:
#>          id   lat     lon azi unc type depth quality regime   lon.PoR  lat.PoR
#> 1  wsm00892 38.14 -118.84  50  25  FMS     7       C      S -85.46364 59.05548
#> 2  wsm00893 35.97 -114.71  54  25  FMS     5       C      S -78.16683 60.50445
#> 3  wsm00894 37.93 -118.17  24  25  FMS     5       C      S -84.55317 59.38025
#> 4  wsm00896 38.63 -118.21  41  25  FMS    17       C      N -85.74431 59.73632
#> 5  wsm00897 39.08 -115.62  30  25  FMS     5       C      N -84.31053 61.68018
#> 6  wsm00903 38.58 -112.58  27  25  FMS     7       C      N -80.60932 63.39714
#> 7  wsm00904 38.71 -112.04  29  25  FMS     9       C      N -80.32836 63.82044
#> 8  wsm00905 37.30 -114.94  58  25  FMS     5       C     NS -80.59498 61.13839
#> 9  wsm00906 37.61 -113.30  32  25  FMS     2       C     TS -79.55181 62.38235
#> 10 wsm00916 37.69 -115.05  47  25  FMS     3       C      S -81.36361 61.28913
#>       azi.PoR                   geometry
#> 1  173.240166 POINT (-85.46364 59.05548)
#> 2    1.058195 POINT (-78.16683 60.50445)
#> 3  147.609558 POINT (-84.55317 59.38025)
#> 4  163.607388 POINT (-85.74431 59.73632)
#> 5  152.233017 POINT (-84.31053 61.68018)
#> 6  150.611386 POINT (-80.60932 63.39714)
#> 7  152.525837 POINT (-80.32836 63.82044)
#> 8    3.075318 POINT (-80.59498 61.13839)
#> 9  156.996990 POINT (-79.55181 62.38235)
#> 10 171.468619 POINT (-81.36361 61.28913)
```
