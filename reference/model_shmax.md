# Theoretical Direction of Maximum Horizontal Stress in the geographical reference system.

Models the direction of maximum horizontal stress \\\sigma\_{Hmax}\\
along great circles, small circles, and loxodromes at a given point or
points according to the relative plate motion in the geographical
coordinate reference system.

## Usage

``` r
model_shmax(df, euler)
```

## Arguments

- df:

  `data.frame` containing the coordinates of the point(s) (`lat`,
  `lon`).

- euler:

  `"data.frame"` or object of class `"euler.pole"` containing the
  geographical coordinates of the Euler pole

## Value

`data.frame`

- gc:

  Azimuth of modeled \\\sigma\_{Hmax}\\ following great circles

- sc:

  Small circles

- ld.cw:

  Clockwise loxodromes

- ld.ccw:

  Counter-clockwise loxodromes

## Details

\\\sigma\_{Hmax}\\ following *great circles* is the (initial) bearing
between the given point and the pole of relative plate motion.
\\\sigma\_{Hmax}\\ along *small circles*, clockwise, and
counter-clockwise *loxodromes* is 90\\^{\circ}\\, +45\\^{\circ}\\, and
135\\^{\circ}\\ (-45\\^{\circ}\\) to this great circle bearing,
respectively.

## References

Stephan, T., Enkelmann, E., and Kroner, U. "Analyzing the horizontal
orientation of the crustal stress adjacent to plate boundaries". *Sci
Rep* 13. 15590 (2023).
[doi:10.1038/s41598-023-42433-2](https://doi.org/10.1038/s41598-023-42433-2)
.

## See also

[`deviation_shmax()`](https://tobiste.github.io/tectonicr/reference/deviation_shmax.md)
to compute the deviation of the modeled direction from the observed
direction of \\\sigma\_{Hmax}\\.
[`PoR_shmax()`](https://tobiste.github.io/tectonicr/reference/PoR_azi.md)
to calculate the azimuth of \\\sigma\_{Hmax}\\ in the pole of rotation
reference system.

## Author

Tobias Stephan

## Examples

``` r
data("nuvel1")
# North America relative to Pacific plate:
euler <- subset(nuvel1, nuvel1$plate.rot == "na")

# the point where we mant to model the SHmax direction:
point <- data.frame(lat = 45, lon = 20)

model_shmax(point, euler)
#>         sc   ld.ccw       gc    ld.cw
#> 1 42.45436 87.45436 132.4544 177.4544
```
