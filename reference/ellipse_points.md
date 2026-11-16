# Vertices of an oriented ellipse

Computes the outline of an ellipse defined by its centre, the length of
its semi-major axis, its aspect ratio, and the azimuth of its long axis.
The coordinate system is assumed to be map-like: `x` points east and `y`
points north.

## Usage

``` r
ellipse_points(x0, y0, a, R, angle, n = 100)
```

## Arguments

- x0, y0:

  Numeric. Coordinates of the ellipse centre.

- a:

  Numeric. Length of the semi-major axis (the full long axis is
  `2 * a`).

- R:

  Numeric. Aspect ratio of the ellipse, \\R = a / b\\, where `b` is the
  semi-minor axis. Must be \\\ge 1\\.

- angle:

  Numeric. Azimuth of the long axis in degrees, measured clockwise from
  geographic north (i.e. from the positive `y` axis). Because an axis
  has no direction, `angle` and `angle + 180` describe the same ellipse.

- n:

  Integer. Number of vertices used to draw the outline.

## Value

A `data.frame` with `n` rows and columns `x` and `y`. The first and last
vertices are the same, so the outline is closed.

## Details

The long-axis unit vector is \\(\sin\phi, \cos\phi)\\ and the short-axis
unit vector is \\(\cos\phi, -\sin\phi)\\, where \\\phi\\ is `angle` in
radians. The points are \$\$\mathbf{p}(t) = \mathbf{c} + a\cos
t\\\mathbf{u} + b\sin t\\\mathbf{v}, \quad t \in \[0, 2\pi\].\$\$

## Examples

``` r
# Long axis N-S: y-range is +/- 1, x-range is +/- 0.5
p <- ellipse_points(0, 0, a = 1, R = 2, angle = 0)
range(p$x); range(p$y)
#> [1] -0.4999371  0.4999371
#> [1] -0.9994965  1.0000000

# Long axis trending 045
plot(ellipse_points(0, 0, a = 1, R = 3, angle = 45), type = "l", asp = 1)

```
