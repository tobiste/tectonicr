# Euler pole object

Creates an object of the orientation of the Euler pole axis

## Usage

``` r
euler_pole(x, y, z = NA, geo = TRUE, angle = NA)
```

## Arguments

- x:

  latitude or x coordinate of Euler pole axis

- y:

  longitude or y

- z:

  z coordinate

- geo:

  logical, `TRUE` (the default) if Euler pole axis is given in
  geographical coordinates (latitude, longitude). `FALSE` if given in
  Cartesian coordinates (`x`, `y`, `z`)

- angle:

  (optional) Angle of rotation in degrees (CCW rotation if angle is
  positive)

## Value

An object of class `"euler.pole"` containing the Euler pole axis in both
geographical and Cartesian coordinates and the angle of rotation in
radians.

## Examples

``` r
euler_pole(90, 0, angle = 45)
#>   lat lon x y z     angle
#> 1  90   0 0 0 1 0.7853982
euler_pole(0, 0, 1, geo = FALSE)
#>   lat lon x y z angle
#> 1  90   0 0 0 1    NA
```
