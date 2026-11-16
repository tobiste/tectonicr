# Azimuth visualization

`geom_azimuth()` visualizes axial-directional vector fields using a geom
to produce a new graphical layer, which allows aesthetic options. This
layer can be overlaid on a map to improve visualisation of mapped data.
The geom draws line segments (spokes) centered at (x, y) with a given
orientation (`angle` in degrees) and length (`radius`). By default the
spoke is centered using
[`PositionCenterSpoke`](https://tobiste.github.io/tectonicr/reference/PositionCenterSpoke.md),
so that the given coordinates mark the middle of the line. The azimuths
are given as angles in degrees increasing clockwise from North.

## Usage

``` r
geom_azimuth(
  mapping = NULL,
  data = NULL,
  stat = "azimuth",
  center = TRUE,
  radius = NULL,
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  ...
)
```

## Arguments

- mapping:

  Set of aesthetic mappings created by
  [`aes()`](https://ggplot2.tidyverse.org/reference/aes.html). If
  specified and `inherit.aes = TRUE` (the default), it is combined with
  the default mapping at the top level of the plot. You must supply
  `mapping` if there is no plot mapping.

- data:

  The data to be displayed in this layer. There are three options:

  If `NULL`, the default, the data is inherited from the plot data as
  specified in the call to
  [`ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html).

  A `data.frame`, or other object, will override the plot data. All
  objects will be fortified to produce a data frame. See
  [`fortify()`](https://ggplot2.tidyverse.org/reference/fortify.html)
  for which variables will be created.

  A `function` will be called with a single argument, the plot data. The
  return value must be a `data.frame`, and will be used as the layer
  data. A `function` can be created from a `formula` (e.g.
  `~ head(.x, 10)`).

- stat:

  The statistical transformation to use on the data for this layer. When
  using a `geom_*()` function to construct a layer, the `stat` argument
  can be used to override the default coupling between geoms and stats.
  The `stat` argument accepts the following:

  - A `Stat` ggproto subclass, for example `StatCount`.

  - A string naming the stat. To give the stat as a string, strip the
    function name of the `stat_` prefix. For example, to use
    [`stat_count()`](https://ggplot2.tidyverse.org/reference/geom_bar.html),
    give the stat as `"count"`.

  - For more information and other ways to specify the stat, see the
    [layer
    stat](https://ggplot2.tidyverse.org/reference/layer_stats.html)
    documentation.

- center:

  Logical; if `TRUE` (the default) spokes are centered on (x, y) using
  [`PositionCenterSpoke`](https://tobiste.github.io/tectonicr/reference/PositionCenterSpoke.md) -
  useful for axial data. If `FALSE`, behaves like
  [`ggplot2::geom_spoke()`](https://ggplot2.tidyverse.org/reference/geom_spoke.html)
  (line starts at (x, y)) - useful for directional data (especially when
  in combination with [`arrow()`](https://rdrr.io/r/grid/arrow.html)).

- radius:

  Length of spoke

- na.rm:

  If `FALSE`, the default, missing values are removed with a warning. If
  `TRUE`, missing values are silently removed.

- show.legend:

  logical. Should this layer be included in the legends? `NA`, the
  default, includes if any aesthetics are mapped. `FALSE` never
  includes, and `TRUE` always includes. It can also be a named logical
  vector to finely select the aesthetics to display. To include legend
  keys for all levels, even when no data exists, use `TRUE`. If `NA`,
  all levels are shown in legend, but unobserved levels are omitted.

- inherit.aes:

  If `FALSE`, overrides the default aesthetics, rather than combining
  with them. This is most useful for helper functions that define both
  data and aesthetics and shouldn't inherit behaviour from the default
  plot specification, e.g.
  [`annotation_borders()`](https://ggplot2.tidyverse.org/reference/annotation_borders.html).

- ...:

  Other arguments passed on to
  [`ggplot2::layer()`](https://ggplot2.tidyverse.org/reference/layer.html).
  These are often aesthetics (e.g. `colour`, `linetype`, `linewidth`,
  `alpha`).

## Value

A ggplot2 layer that adds axis-like spokes.

## Aesthetics

`geom_azimuth()` understands the following aesthetics (required
aesthetics in **bold**):

- **x**

- **y**

- angle (in degrees, measured clockwise from North; transformed
  internally)

- radius

- colour

- alpha

- linewidth

- linetype

## See also

[`ggplot2::geom_spoke()`](https://ggplot2.tidyverse.org/reference/geom_spoke.html),
[`geom_azimuthpoint()`](https://tobiste.github.io/tectonicr/reference/geom_azimuthpoint.md)

## Examples

``` r
set.seed(20250411)
df <- data.frame(
  x = runif(5), y = runif(5),
  angle_deg = rvm(5, mean = 90, kappa = 10),
  radius = runif(5, 0.1, 2)
)

if (require("ggplot2")) {
  ggplot(df, aes(x, y)) +
    geom_azimuth(aes(angle = angle_deg), radius = .1, linewidth = 1.2, colour = "blue")
  if (require("grid")) {
    ggplot(df, aes(x, y, radius = radius)) +
      geom_azimuth(aes(angle = angle_deg), center = FALSE, colour = "red", arrow = grid::arrow())
  }
}
#> Loading required package: ggplot2
#> Loading required package: grid
```
