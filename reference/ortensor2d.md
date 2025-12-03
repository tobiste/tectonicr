# Orientation Tensor

2D orientation tensor characterizes distribution of axial angles using
the Eigenvalue method (Watson 1966, Scheidegger 1965).

## Usage

``` r
ortensor2d(x, w = NULL, norm = FALSE)
```

## Arguments

- x:

  numeric. Axial angular data (in degrees).

- w:

  (optional) Weights. A vector of positive numbers and of the same
  length as `x`.

- norm:

  logical. Whether the tensor should be normalized.

## Value

2x2 matrix

## Details

The moment of inertia can be minimized by calculating the Cartesian
coordinates of the orientation data, and calculating their covariance
matrix. This yields \$\$I = x \cdot x^\intercal\$\$ where \\x\\ is the
Cartesian vector of the orientations. Orientation tensor \\T\\ and the
inertia tensor \\I\\ are related by \$\$I = E - T\$\$ where \\E\\
denotes the unit matrix, so that \$\$T = \frac{1}{n} \sum\_{i=i}^{n} x_i
\cdot x_i^\intercal\$\$

## References

Watson, G. S. (1966). The Statistics of Orientation Data. The Journal of
Geology, 74(5), 786–797.

Scheidegger, A. E. (1964). The tectonic stress and tectonic motion
direction in Europe and Western Asia as calculated from earthquake fault
plane solutions. Bulletin of the Seismological Society of America,
54(5A), 1519–1528. doi:10.1785/BSSA05405A1519

Bachmann, F., Hielscher, R., Jupp, P. E., Pantleon, W., Schaeben, H., &
Wegert, E. (2010). Inferential statistics of electron backscatter
diffraction data from within individual crystalline grains. Journal of
Applied Crystallography, 43(6), 1338–1355.
https://doi.org/10.1107/S002188981003027X

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
