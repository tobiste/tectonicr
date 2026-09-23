# Circular plot

Polar diagram where angles increase clockwise.

## Usage

``` r
circular_plot(
  main = NULL,
  sub = NULL,
  labels = TRUE,
  at = c(grid.lines, grid.lines + 180),
  cborder = TRUE,
  radius = 1,
  scale = 1.1,
  grid = FALSE,
  grid.lines = seq(0, 135, 45),
  grid.circles = seq(0.2, radius, 0.2),
  origin.text = "N",
  ...
)
```

## Arguments

- main:

  a main title for the plot, see also
  [`title`](https://rdrr.io/r/graphics/title.html).

- sub:

  a subtitle for the plot.

- labels:

  Either a logical value indicating whether to plot labels next to the
  tick marks, or a vector of labels for the tick marks.

- at:

  Optional vector of angles at which tick marks should be plotted. Set
  `at=numeric(0)` to suppress tick marks.

- cborder:

  logical. Whether to plot the primitive of the circular plot (`TRUE` by
  default) or not.

- radius:

  numeric. The radius of the circular plot. `1` by default.

- scale:

  numeric. Increases the radius of the border to have some in-plot
  clearance. `1.1` by default

- grid:

  logical. Whether to add a grid to the circular plot. Default is
  `FALSE`.

- grid.lines, grid.circles:

  numeric. Adds a sequence of straight grid lines and circles based on
  angles (in degrees) and radii, respectively. Ignored when `grid=FALSE`

- origin.text:

  character. String to be drawn at the top margin of the plot (`"N"` by
  default)

- ...:

  optional parameters passed to
  [`graphics::symbols()`](https://rdrr.io/r/graphics/symbols.html) for
  plotting border of circular plot. Ignored if `cborder=FALSE`.

## Value

none

## Examples

``` r
circular_plot(grid = TRUE)
```
