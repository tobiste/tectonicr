# Distance between points

Returns the great circle distance between a location and all grid point
in km

## Usage

``` r
dist_greatcircle(
  lat1,
  lon1,
  lat2,
  lon2,
  r = earth_radius(),
  method = c("haversine", "orthodrome", "vincenty", "euclidean")
)
```

## Arguments

- lat1, lon1:

  numeric vector. coordinate of point(s) 1 (degrees).

- lat2, lon2:

  numeric vector. coordinates of point(s) 2 (degrees).

- r:

  numeric. radius of the sphere (default = 6371.0087714 km, i.e. the
  radius of the Earth)

- method:

  Character. Formula for calculating great circle distance, one of:

  `"haversine"`

  :   great circle distance based on the haversine formula that is
      optimized for 64-bit floating-point numbers (the default)

  `"orthodrome"`

  :   great circle distance based on the spherical law of cosines

  `"vincenty"`

  :   distance based on the Vincenty formula for an ellipsoid with equal
      major and minor axes

  "euclidean"

  :   Euclidean distance (not great circle distance!)

## Value

numeric vector with length equal to `length(lat1)` or `length(lat2)`

## See also

[`orthodrome()`](https://tobiste.github.io/tectonicr/reference/spherical_angle.md),
[`haversine()`](https://tobiste.github.io/tectonicr/reference/spherical_angle.md),
[`vincenty()`](https://tobiste.github.io/tectonicr/reference/spherical_angle.md)

## Examples

``` r
# Haversine: (4149.157, 2296.583) km
dist_greatcircle(lat1 = 20, lon1 = 12, lat2 = c(50, 30), lon2 = c(40, 32))
#> [1] 4149.157 2296.583

# Orthodrome: (4149.157, 2296.583) km
dist_greatcircle(
  lat1 = 20, lon1 = 12, lat2 = c(50, 30), lon2 = c(40, 32),
  method = "orthodrome"
)
#> [1] 4149.157 2296.583

# Vincenty: (4149.157, 2296.583) km
dist_greatcircle(
  lat1 = 20, lon1 = 12, lat2 = c(50, 30), lon2 = c(40, 32),
  method = "vincenty"
)
#> [1] 4149.157 2296.583

# Euclidean (4076.220, 2284.169) km
dist_greatcircle(
  lat1 = 20, lon1 = 12, lat2 = c(50, 30), lon2 = c(40, 32),
  method = "euclidean"
)
#> [1] 4076.220 2284.169
```
