# Coordinate Transformations

Converts vector between Cartesian and spherical coordinate systems

## Usage

``` r
cartesian_to_spherical(n)

spherical_to_cartesian(p)

spherical_to_geographical(p)
```

## Arguments

- n:

  Cartesian coordinates (x, y, z) as three-column vector

- p:

  Spherical coordinates (colatitude, azimuth) as two-column vector

## Value

Functions return a (2- or 3-dimensional) vector representing a point in
the requested coordinate system.

## See also

[`cartesian_to_geographical()`](https://tobiste.github.io/tectonicr/reference/coordinates.md)
and
[`geographical_to_cartesian()`](https://tobiste.github.io/tectonicr/reference/coordinates.md)
for conversions to geographical coordinates

## Examples

``` r
n <- c(1, -2, 3)
cartesian_to_spherical(n) # 36.699, -63.435
#> [1]  36.69923 -63.43495
p <- c(50, 10)
spherical_to_cartesian(p) # 0.75, 0.13, 0.64
#> [1] 0.7544065 0.1330222 0.6427876
```
