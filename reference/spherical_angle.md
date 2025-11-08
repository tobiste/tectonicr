# Angle along great circle on spherical surface

Smallest angle between two points on the surface of a sphere, measured
along the surface of the sphere

## Usage

``` r
orthodrome(lat1, lon1, lat2, lon2)

haversine(lat1, lon1, lat2, lon2)

vincenty(lat1, lon1, lat2, lon2)
```

## Arguments

- lat1, lat2:

  numeric vector. latitudes of point 1 and 2 (in radians)

- lon1, lon2:

  numeric vector. longitudes of point 1 and 2 (in radians)

## Value

numeric. Angle in radians

## Details

- `"orthodrome"`:

  based on the spherical law of cosines

- `"haversine"`:

  uses haversine formula that is optimized for 64-bit floating-point
  numbers

- `"vincenty"`:

  uses Vincenty formula for an ellipsoid with equal major and minor axes

## References

- Imboden, C. & Imboden, D. (1972). Formel fuer Orthodrome und Loxodrome
  bei der Berechnung von Richtung und Distanz zwischen Beringungs- und
  Wiederfundort. *Die Vogelwarte* **26**, 336-346.

- Sinnott, Roger W. (1984). Virtues of the Haversine. *Sky and
  telescope* **68**(2), 158. Vincenty, T. (1975). Direct and inverse
  solutions of geodesics on the ellipsoid with application of nested
  equations. *Survey Review*, **23**(176), 88\<U+2013\>93.
  [doi:10.1179/sre.1975.23.176.88](https://doi.org/10.1179/sre.1975.23.176.88)
  .

- [http://www.movable-type.co.uk/scripts/latlong.html](http://www.movable-type.co.uk/scripts/latlong.md)

- <http://www.edwilliams.org/avform147.htm>

## Examples

``` r
berlin <- c(52.52, 13.41) |> deg2rad()
calgary <- c(51.04, -114.072) |> deg2rad()
orthodrome(berlin[1], berlin[2], calgary[1], calgary[2]) # 1.176406
#> [1] 1.176404
haversine(berlin[1], berlin[2], calgary[1], calgary[2]) # 1.176406
#> [1] 1.176404
vincenty(berlin[1], berlin[2], calgary[1], calgary[2]) # 1.176406
#> [1] 1.176404
```
