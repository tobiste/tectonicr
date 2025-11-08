# Relative rotation between two rotations

Calculates the relative rotation between two rotations, i.e. the
difference from rotation 1 to rotation 2.

## Usage

``` r
relative_rotation(r1, r2)
```

## Arguments

- r1, r2:

  Objects of class `"euler.pole"`. First rotation is `r1`, followed
  rotation `r2`.

## Value

`list`. Euler axes (geographical coordinates) and Euler angles (in
degrees)

## References

Schaeben, H., Kroner, U. and Stephan, T. (2021). Euler Poles of Tectonic
Plates. In B. S. Daza Sagar, Q. Cheng, J. McKinley and F. Agterberg
(Eds.), *Encyclopedia of Mathematical Geosciences. Encyclopedia of Earth
Sciences Series* (pp. 1–7). Springer Nature Switzerland AG 2021. doi:
10.1007/978-3-030-26050-7_435-1.

## See also

[`euler_pole()`](https://tobiste.github.io/tectonicr/reference/euler_pole.md)
for class `"euler.pole"`

## Examples

``` r
a <- euler_pole(90, 0, angle = 45)
b <- euler_pole(0, 0, 1, geo = FALSE, angle = -15)
relative_rotation(a, b) # axis: -90, -180; angle: 60
#> $axis
#> [1]  -90 -180
#> 
#> $angle
#> [1] 60
#> 
relative_rotation(b, a) # axis: 90, 0; angle: 60
#> $axis
#> [1] 90  0
#> 
#> $angle
#> [1] 60
#> 
```
