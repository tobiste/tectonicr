# Normalize Angle Between Two Directions

Normalizes the angle between two directions to the acute angle in
between, i.e. angles between 0 and 90\\^\circ\\

## Usage

``` r
deviation_norm(x, y = NULL)
```

## Arguments

- x, y:

  Minuend and subtrahend. Both numeric vectors of angles in degrees. If
  `y` is missing, it treats `x` as difference. If not, length of
  subtrahend `y` is either `1` or equal to `length(x)`.

## Value

numeric vector, acute angles between two directions, i.e. values between
0 and 90\\^\circ\\

## Author

Tobias Stephan

## Examples

``` r
deviation_norm(175, 5)
#> [1] 10
deviation_norm(c(175, 95, 0), c(5, 85, NA))
#> [1] 10 10 NA
deviation_norm(c(-5, 85, 95, 175, 185, 265, 275, 355, 365))
#> [1]  5 85 85  5  5 85 85  5  5
```
