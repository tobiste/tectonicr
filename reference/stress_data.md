# Example crustal stress dataset

Subsets of the World Stress Map (WSM) compilation of information on the
crustal present-day stress field (Version 1.1. 2019).

## Usage

``` r
data('san_andreas')

data('tibet')

data('iceland')
```

## Format

A `sf` object / `data.frame` with 10 columns. Each row represents a
different in-situ stress measurement:

- id:

  Measurement identifier

- lat:

  Latitude in degrees

- lon:

  Longitude in degrees

- azi:

  SHmax azimuth in degrees

- unc:

  Measurement standard deviation (in degrees)

- type:

  Type of measurement

- depth:

  Depth in km

- quality:

  WSM quality rank

- regime:

  Stress regime

An object of class `sf` (inherits from `data.frame`) with 1126 rows and
10 columns.

An object of class `sf` (inherits from `data.frame`) with 1165 rows and
10 columns.

An object of class `sf` (inherits from `data.frame`) with 490 rows and
10 columns.

## Source

<https://www.world-stress-map.org/>

## Details

- `'san_andreas"`:

  contains 407 stress data adjacent to the San Andreas Fault to be
  tested against a tangentially displaced plate boundary.

- `"tibet"`:

  contains 947 stress data from the Himalaya and Tibetan plateau to be
  tested against an inward-moving displaced plate boundary.

- `'iceland`:

  contains 201 stress data from Iceland to be tested against a
  outward-moving displaced plate boundary.

## References

Heidbach, O., Barth, A., Müller, B., Reinecker, J., Stephansson, O.,
Tingay, M., & Zang, A. (2016). WSM quality ranking scheme, database
description and analysis guidelines for stress indicator. WSM Technical
Report; 16-01. GFZ German Research Centre for Geosciences.
[doi:10.2312/WSM.2016.001](https://doi.org/10.2312/WSM.2016.001)

## See also

[`download_WSM()`](https://tobiste.github.io/tectonicr/reference/import_WSM.md)
for description of columns and stress regime acronyms

## Examples

``` r
data("san_andreas")
head(san_andreas)
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

data("tibet")
head(tibet)
#> Simple feature collection with 6 features and 9 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 73.7 ymin: 35.4 xmax: 77.4 ymax: 42.63
#> Geodetic CRS:  WGS 84
#>         id   lat   lon azi unc type depth quality regime            geometry
#> 1 wsm01766 35.40 73.70 100  40  FMS     0       D      T   POINT (73.7 35.4)
#> 2 wsm01860 42.10 75.98  62  25  FMS    10       C      T  POINT (75.98 42.1)
#> 3 wsm01861 42.00 77.40  18  25  FMS     0       C     NS     POINT (77.4 42)
#> 4 wsm01869 42.57 75.08   5  25  FMS    10       C      T POINT (75.08 42.57)
#> 5 wsm01877 40.90 73.93 177  25  FMS    10       C      T  POINT (73.93 40.9)
#> 6 wsm01880 42.63 75.33 167  25  FMS     5       C      T POINT (75.33 42.63)

data("iceland")
head(iceland)
#> Simple feature collection with 6 features and 9 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: -16.99 ymin: 64.38 xmax: -15.16 ymax: 66.33
#> Geodetic CRS:  WGS 84
#>         id   lat    lon azi unc type depth quality regime             geometry
#> 1 wsm00187 64.38 -15.16 167  25  FMS  26.0       C      N POINT (-15.16 64.38)
#> 2 wsm00188 66.33 -16.29 172  25  FMS  33.0       C      S POINT (-16.29 66.33)
#> 3 wsm00190 66.20 -16.70 174  25  FMS   0.0       C      S   POINT (-16.7 66.2)
#> 4 wsm00191 64.90 -16.94 164  25  FMS  10.0       C      T  POINT (-16.94 64.9)
#> 5 wsm00192 64.57 -16.95 140  90  FMS  13.8       E      N POINT (-16.95 64.57)
#> 6 wsm00193 64.66 -16.99 129  90  FMS  14.5       E      N POINT (-16.99 64.66)
```
