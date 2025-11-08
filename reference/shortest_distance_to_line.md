# Shortest distance between pairs of geometries

The shortest Great Circle distance between pairs of geometries

## Usage

``` r
shortest_distance_to_line(x, line, ellipsoidal = FALSE)
```

## Arguments

- x, line:

  objects of class `sfg`, `sfc` or `sf`

- ellipsoidal:

  Logical. Whether the distance is calculated using spherical distances
  ([`sf::st_distance()`](https://r-spatial.github.io/sf/reference/geos_measures.html))
  or ellipsoidal distances (`lwgeom::st_geod_distance()`).

## Value

numeric. Shortest distance in meters

## Examples

``` r
plate_boundary <- subset(plates, plates$pair == "na-pa")
shortest_distance_to_line(san_andreas, plate_boundary) |>
  head()
#> [1] 270330.7 278345.4 299623.0 346884.7 547464.8 622915.5
```
