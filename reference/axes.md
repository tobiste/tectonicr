# Plot axes

Show direction axes in a map

## Usage

``` r
axes(
  x,
  y,
  angle,
  radius = 0.5,
  arrow.code = 1,
  arrow.length = 0,
  add = FALSE,
  ...
)
```

## Arguments

- x, y:

  coordinates of points

- angle:

  Azimuth in degrees

- radius:

  length of axis

- arrow.code:

  integer. Kind of arrow head. The default is `1`, i.e. no arrow head.
  See [`graphics::arrows()`](https://rdrr.io/r/graphics/arrows.html) for
  details

- arrow.length:

  numeric Length of the edges of the arrow head (in inches). (Ignored if
  `arrow.code = 1`)

- add:

  logical. add to existing plot?

- ...:

  optional arguments passed to
  [`graphics::arrows()`](https://rdrr.io/r/graphics/arrows.html)

## Value

No return value, called for side effects

## Examples

``` r
data("san_andreas")
axes(san_andreas$lon, san_andreas$lat, san_andreas$azi, add = FALSE)
```
