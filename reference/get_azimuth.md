# Azimuth Between two Points

Calculate initial bearing (or forward azimuth/direction) to go from
point `a` to point `b` following great circle arc on a sphere.

## Usage

``` r
get_azimuth(lat_a, lon_a, lat_b, lon_b)
```

## Arguments

- lat_a, lat_b:

  Numeric. Latitudes of a and b (in degrees).

- lon_a, lon_b:

  Numeric. Longitudes of a and b (in degrees).

## Value

numeric. Azimuth in degrees

## Details

`get_azimuth()` is based on the spherical law of tangents. This formula
is for the initial bearing (sometimes referred to as forward azimuth)
which if followed in a straight line along a great circle arc will lead
from the start point `a` to the end point `b`. \$\$\theta = \arctan2
(\sin \Delta\lambda \cos\psi_2, \cos\psi_1 \sin\psi_1-\sin\psi_1
\cos\psi_2 \cos\Delta\lambda)\$\$ where \\\psi_1, \lambda_1\\ is the
start point, \\\psi_2\\, \\\lambda_2\\ the end point (\\\Delta\lambda\\
is the difference in longitude).

## References

[http://www.movable-type.co.uk/scripts/latlong.html](http://www.movable-type.co.uk/scripts/latlong.md)

## Examples

``` r
berlin <- c(52.517, 13.4) # Berlin
tokyo <- c(35.7, 139.767) # Tokyo
get_azimuth(berlin[1], berlin[2], tokyo[1], tokyo[2]) # 41.57361
#> [1] 41.57361
```
