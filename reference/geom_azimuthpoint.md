# Azimuth + point visualization

`geom_azimuthpoint()` draws line segments (spokes) like
[`geom_azimuth()`](https://tobiste.github.io/tectonicr/reference/geom_azimuth.md),
but also places a point (marker) at the spoke's center `(x, y)`.

Aesthetic rules:

- `linewidth`, `linetype` affect the spoke only

- `shape` affects the point only

- `colour`, `alpha` affect both spoke and point

- `size` sets the size of the point only

## Usage

``` r
geom_azimuthpoint(
  mapping = NULL,
  data = NULL,
  stat = "identity",
  center = TRUE,
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  size = 2,
  ...
)
```

## Arguments

- mapping:

  Set of aesthetic mappings created by
  [`ggplot2::aes()`](https://ggplot2.tidyverse.org/reference/aes.html).

- data:

  A data frame. If `NULL`, the default, the data is inherited from the
  plot data as specified in the call to
  [`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html).

- stat:

  The statistical transformation to use on the data. Defaults to
  `"identity"`.

- center:

  Logical; if `TRUE` spokes are centered on (x, y) using
  [`PositionCenterSpoke`](https://tobiste.github.io/tectonicr/reference/PositionCenterSpoke.md).
  If `FALSE`, behaves like
  [`ggplot2::geom_spoke()`](https://ggplot2.tidyverse.org/reference/geom_spoke.html)
  (line starts at (x, y)).

- na.rm:

  If `FALSE`, the default, missing values are removed with a warning. If
  `TRUE`, missing values are silently removed.

- show.legend:

  Logical. Should this layer be included in the legends?

- inherit.aes:

  If `FALSE`, overrides the default aesthetics, rather than combining
  with them.

- size:

  Size of the point marker (default = 2).

- ...:

  Other arguments passed on to
  [`ggplot2::geom_spoke()`](https://ggplot2.tidyverse.org/reference/geom_spoke.html)
  and
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html).
  These may include `arrow`, `fill`, etc.

## Value

A list of ggplot2 layers (spokes + points).

## Aesthetics

`geom_azimuthpoint()` understands the following aesthetics (required
aesthetics in **bold**):

- **x**

- **y**

- angle (in degrees, transformed internally; spoke only)

- radius (spoke only)

- colour (shared)

- alpha (shared)

- linewidth (spoke only)

- linetype (spoke only)

- shape (point only)

- size (point only, or via argument)

- fill (point only, for shapes that accept fill)

## See also

[`geom_azimuth()`](https://tobiste.github.io/tectonicr/reference/geom_azimuth.md),
[`ggplot2::geom_spoke()`](https://ggplot2.tidyverse.org/reference/geom_spoke.html),
[`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html)

## Examples

``` r
set.seed(20250411)
df <- data.frame(
  x = runif(5), y = runif(5),
  angle_deg = rvm(5, mean = 90, kappa = 10),
  radius = runif(5, 0.5, 2),
  group = rep(1:2, length.out = 5)
)

if (require("ggplot2")) {
ggplot(df, aes(x, y)) +
  geom_azimuthpoint(aes(angle = angle_deg, radius = radius,
                    colour = factor(group), shape = factor(group)),
                linewidth = 1.1, linetype = "dashed",
                size = 3, alpha = 0.8)
}

```
