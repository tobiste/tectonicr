# Coordinate Transformations

Converts vector between Cartesian and geographical coordinate systems

## Usage

``` r
cartesian_to_geographical(n)

geographical_to_cartesian(p)

geographical_to_spherical(p)
```

## Arguments

- n:

  Cartesian coordinates (x, y, z) as vector

- p:

  Geographical coordinates (latitude, longitude) as vector

## Value

Functions return a (2- or 3-dimensional) vector representing a point in
the requested coordinate system.

## See also

[`cartesian_to_spherical()`](https://tobiste.github.io/tectonicr/reference/coordinates2.md)
and
[`spherical_to_cartesian()`](https://tobiste.github.io/tectonicr/reference/coordinates2.md)
for conversions to spherical coordinates

## Examples

``` r
n <- c(1, -2, 3)
cartesian_to_geographical(n)
#> [1]  53.30077 -63.43495
p <- c(50, 10)
geographical_to_cartesian(p)
#> [1] 0.6330222 0.1116189 0.7660444
```
