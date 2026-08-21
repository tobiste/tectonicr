# Circular Kernel Density Plot

Plots multiples of a von Mises density or wrapped Cauchy distribution in
a circular plot

## Usage

``` r
plot_density(
  x,
  bw = NULL,
  kernel = c("vonmises", "wrappedcauchy"),
  axial = TRUE,
  n = 512L,
  norm.density = TRUE,
  kappa = NULL,
  rho = NULL,
  ...,
  fill = FALSE,
  scale = 0,
  shrink = 1,
  add = TRUE,
  main = NULL,
  labels = TRUE,
  at = seq(0, 360 - 45, 45),
  cborder = TRUE,
  grid = FALSE
)
```

## Arguments

- x:

  numeric. Data to be plotted, i.e. vector containing angles (in
  degrees).

- bw, kappa, rho:

  numeric. Smoothing bandwidth expressed as the concentration parameter
  \\\kappa\\ for the von Mises distribution or \\\rho\\ for the wrapped
  Cauchy distribution. Small and large values for the von Mises and
  wrapped Cauchy distribution, respectively, gives smooth density lines.
  If not specified, parameter will be estimated using
  [`est.kappa()`](https://tobiste.github.io/tectonicr/reference/estimate-kappa.md)
  for the von Mises distribution, or set to \\p \exp(-1)\\ for the
  wrapped Cauchy distribution (where \\p = 2\\ when `axial=TRUE` and 1
  otherwise).

- kernel:

  character. The smoothing kernel to be used; one of `"vonmises"` (the
  default) or `"wrappedcauchy"` for the von Mises or the Wrapped Cauchy
  distribution.

- axial:

  Logical. Whether data are uniaxial (`axial=FALSE`) or biaxial (`TRUE`,
  the default).

- n:

  integer. the number of equally spaced points at which the density is
  to be estimated.

- norm.density:

  logical. Normalize the density?

- ...:

  Further graphical parameters may also be supplied as arguments.

- fill:

  logical. Whether to fill the density curve or draw just a line (the
  default)

- scale:

  numeric. radius of plotted circle. Default is `1.1`.

- shrink:

  numeric. parameter that controls the size of the plotted function.
  Default is `1`.

- add:

  logical. Add to existing plot? (`TRUE` by default).

- main:

  Character string specifying the title of the plot.

- labels:

  Either a logical value indicating whether to plot labels next to the
  tick marks, or a vector of labels for the tick marks.

- at:

  Optional vector of angles at which tick marks should be plotted. Set
  `at=numeric(0)` to suppress tick marks.

- cborder:

  logical. Border of rose plot.

- grid:

  logical. Whether a grid should be added.

## Value

plot or calculated densities as numeric vector

## See also

[`dvm()`](https://tobiste.github.io/tectonicr/reference/vonmises.md) and
[`dwcauchy()`](https://tobiste.github.io/tectonicr/reference/wcauchy.md)

Other rose-plot:
[`plot_points()`](https://tobiste.github.io/tectonicr/reference/plot_points.md),
[`rose()`](https://tobiste.github.io/tectonicr/reference/rose.md),
[`rose_geom`](https://tobiste.github.io/tectonicr/reference/rose_geom.md),
[`rose_stats()`](https://tobiste.github.io/tectonicr/reference/rose_stats.md)

## Examples

``` r
# Filled von Mises kernel density curve inside the plot
plot_density(san_andreas$azi,
  kappa = 100,
  fill = TRUE, col = "#51127C80", border = "#51127CFF",
  grid = TRUE,
  add = FALSE
)

# Superimpose a wrapped Cauchy kernel distribution curve
plot_density(san_andreas$azi,
  rho = 0.9, kernel = "wrappedcauchy",
  fill = FALSE, col = "#FB8861FF",
  add = TRUE
)


# Superimpose a von Mises kernel density curve on a rose diagram:
rose(san_andreas$azi, grid = TRUE)
plot_density(san_andreas$azi,
  bw = 100, col = "#51127CFF",
  add = TRUE, lwd = 3
)


# Corona plot (censity curve outside of a rose diagram plot):
rose(san_andreas$azi, dots = TRUE, stack = TRUE, dot_cex = 0.5, dot_pch = 21)
plot_density(san_andreas$azi,
  bw = 100,
  scale = 1.1, shrink = 3, xpd = NA,
  col = "#51127CFF"
)
```
