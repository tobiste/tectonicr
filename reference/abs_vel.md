# Absolute Plate Velocity

Calculates the absolute angular velocity of plate motion

## Usage

``` r
abs_vel(w, alpha, r = earth_radius())
```

## Arguments

- w:

  Angular velocity or rate or angle of rotation

- alpha:

  Angular distance to Euler pole or small circle around Euler pole

- r:

  Radius. Default is WGS84 Earth's radius (6371.009 km)

## Value

numeric (unit of velocity: km/Myr)

## See also

[`earth_radius()`](https://tobiste.github.io/tectonicr/reference/earth_radius.md)

## Examples

``` r
abs_vel(0.21, 0)
#> [1] 0
abs_vel(0.21, 45)
#> [1] 16.51163
abs_vel(0.21, 90)
#> [1] 23.35097
```
