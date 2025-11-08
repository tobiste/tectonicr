# Angle Between Two Vectors

Calculates the angle between two vectors

## Usage

``` r
angle_vectors(x, y)
```

## Arguments

- x, y:

  Vectors in Cartesian coordinates. Can be vectors of three numbers or a
  matrix of 3 columns (x, y, z)

## Value

numeric. angle in degrees

## Examples

``` r
u <- c(1, -2, 3)
v <- c(-2, 1, 1)
angle_vectors(u, v) # 96.26395
#> [1] 96.26395
```
