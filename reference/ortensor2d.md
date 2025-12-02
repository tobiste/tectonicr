# Orientation Tensor

2D orientation tensor, which characterizes data distribution using the
Eigenvalue (Watson 1966, Scheidegger 1965).

## Usage

``` r
ortensor2d(x, w = NULL, norm = FALSE)
```

## Arguments

- x:

  numeric vector. Values in degrees.

- w:

  (optional) Weights. A vector of positive numbers and of the same
  length as `x`.

- norm:

  logical. Whether the tensor should be normalized.

## Value

2x2 matrix

## Note

\$\$E = x \cdot x^{T}\$\$ where \\x\\ is the Cartesian vector of the
orientations.

Orientation tensor is also called "inertia tensor".

## References

Watson, G. S. (1966). The Statistics of Orientation Data. The Journal of
Geology, 74(5), 786–797.

Scheidegger, A. E. (1964). The tectonic stress and tectonic motion
direction in Europe and Western Asia as calculated from earthquake fault
plane solutions. Bulletin of the Seismological Society of America,
54(5A), 1519–1528. doi:10.1785/BSSA05405A1519

## See also

[`ot_eigen2d()`](https://tobiste.github.io/tectonicr/reference/ort-eigen.md)

## Examples

``` r
test <- rvm(100, mean = 0, k = 10)
ortensor2d(test)
#>            [,1]       [,2]
#> [1,] 0.90068830 0.01410628
#> [2,] 0.01410628 0.09931170
```
