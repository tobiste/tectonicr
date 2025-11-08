# Quantile-Quantile Linearised Plot for Circular Distributions

Uniformly distributed orientations should yield a straight line through
the origin. Systematic departures from linearity will indicate preferred
orientation.

## Usage

``` r
circular_qqplot(
  x,
  axial = TRUE,
  xlab = paste("i/(n+1)"),
  ylab = NULL,
  main = "Circular Quantile-Quantile Plot",
  add_line = TRUE,
  col = "#B63679FF",
  ...
)
```

## Arguments

- x:

  numeric. Angles in degrees

- axial:

  Logical. Whether data are uniaxial (`axial=FALSE`)

- xlab, ylab, main:

  plot labels.

- add_line:

  logical. Whether to connect the points by straight lines?

- col:

  color for the dots.

- ...:

  graphical parameters

## Value

plot

## References

Borradaile, G. J. (2003). Statistics of earth science data: their
distribution in time, space, and orientation (Vol. 351, p. 329). Berlin:
Springer.

## Examples

``` r
# von Mises distribution
x_vm <- rvm(100, mean = 0, kappa = 2)
circular_qqplot(x_vm, pch = 20)


x_norm <- rnorm(100, mean = 0, sd = 25)
circular_qqplot(x_norm, pch = 20)


# uniform (random) data
x_unif <- runif(100, 0, 360)
circular_qqplot(x_unif, pch = 20)
```
