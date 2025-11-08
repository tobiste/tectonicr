# von Mises Quantile-Quantile Plot

Produces a Q-Q plot of the data against a specified von Mises
distribution to graphically assess the goodness of fit of the model.

## Usage

``` r
vm_qqplot(
  x,
  w = NULL,
  axial = TRUE,
  mean = NULL,
  kappa = NULL,
  xlab = "von Mises quantile function",
  ylab = "Empirical quantile function",
  main = "von Mises Q-Q Plot",
  col = "#B63679FF",
  add_line = TRUE,
  ...
)
```

## Arguments

- x:

  numeric. Angles in degrees

- w:

  numeric. optional weightings for `x` to estimate `mean` and `kappa`.

- axial:

  Logical. Whether data are uniaxial (`axial=FALSE`)

- mean:

  numeric. Circular mean of the von Mises distribution. If `NULL`, it
  will be estimated from `x`.

- kappa:

  numeric. Concentration parameter of the von Mises distribution. If
  `NULL`, it will be estimated from `x`.

- xlab, ylab, main:

  plot labels.

- col:

  color for the dots.

- add_line:

  logical. Whether to connect the points by straight lines?

- ...:

  graphical parameters

## Value

plot

## Examples

``` r
# von Mises distribution
x_vm <- rvm(100, mean = 0, kappa = 4)
vm_qqplot(x_vm, axial = FALSE, pch = 20)


# uniform distribution
x_unif <- runif(100, 0, 360)
vm_qqplot(x_unif, axial = FALSE, pch = 20)
```
