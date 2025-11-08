# Direction Lines and Fans in Circular Diagram

Direction Lines and Fans in Circular Diagram

## Usage

``` r
rose_line(x, radius = 1, axial = TRUE, add = TRUE, ...)

rose_fan(x, d, radius = 1, axial = TRUE, add = TRUE, ...)
```

## Arguments

- x:

  angles in degrees

- radius:

  of the plotted circle

- axial:

  Logical. Whether `x` are uniaxial (`axial=FALSE`) or biaxial (`TRUE`,
  the default).

- add:

  logical. Add to existing plot?

- ...:

  optional arguments passed to
  [`graphics::segments()`](https://rdrr.io/r/graphics/segments.html) or
  [`graphics::polygon()`](https://rdrr.io/r/graphics/polygon.html)

- d:

  width of a fan (in degrees)

## Value

No return value, called for side effects

## See also

Other rose-plot:
[`plot_density()`](https://tobiste.github.io/tectonicr/reference/plot_density.md),
[`plot_points()`](https://tobiste.github.io/tectonicr/reference/plot_points.md),
[`rose()`](https://tobiste.github.io/tectonicr/reference/rose.md),
[`rose_stats()`](https://tobiste.github.io/tectonicr/reference/rose_stats.md)

## Examples

``` r
angles <- c(0, 10, 45)
radius <- c(.7, 1, .2)
lwd <- c(2, 1, .75)
col <- c(1, 2, 3)
rose_line(angles, radius = radius, axial = FALSE, add = FALSE, lwd = lwd, col = col)
```
